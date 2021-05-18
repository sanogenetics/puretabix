import errno
import io
import itertools
import logging
import os
import struct
import zlib

try:
    # this supports e.g. partial funcs, lambdas
    import multiprocessing_on_dill as multiprocessing
except ImportError:
    import multiprocessing

logger = logging.getLogger(__name__)

headerpattern = "<BBBBIBBHBBHH"
headersize = struct.calcsize(headerpattern)
tailpattern = "<II"
tailsize = struct.calcsize(tailpattern)


def check_is_header(header):
    if header[0] != 31:
        return False
    if header[1] != 139:
        return False
    if header[2] != 8:
        return False
    if header[3] != 4:
        return False
    if header[8] != 66:
        return False
    if header[9] != 67:
        return False
    if header[10] != 2:
        return False
    return True


def get_cdata(header, infile):
    blocksize = header[11] - header[7] - 19
    cdata = infile.read(blocksize)
    assert len(cdata) == blocksize, f"Unable to read up to {blocksize} of cdata"
    return cdata


def get_tail(infile, decompressed=None):
    # read isize and crc check
    tailbytes = infile.read(tailsize)
    if len(tailbytes) != tailsize:
        raise ValueError(f"Unable to read {tailsize} bytes for tail")
    tail_crc, tail_isize = struct.unpack(tailpattern, tailbytes)
    # if we were given the decompressed data, check it matches expectation
    if decompressed:
        # check decompressed size is expected
        assert len(decompressed) == tail_isize
        # check crc check is expected
        assert zlib.crc32(decompressed) == tail_crc
    return tail_crc, tail_isize


def get_cdata_decompressed(header, infile):
    cdata = get_cdata(header, infile)
    # now do the actual decompression
    decompressor = zlib.decompressobj(
        wbits=-15
    )  # we've alread read the header, so ignore it
    decompressed = decompressor.decompress(cdata)
    assert (
        not decompressor.unconsumed_tail
    ), f"unconsumed tail of {len(decompressor.unconsumed_tail)}"
    assert (
        not decompressor.unused_data
    ), f"unused data present of {len(decompressor.unused_data)}"
    return cdata, decompressed


def blocks_in_file_slice(infile, start, end):
    assert start >= 0
    assert end >= 0
    assert start < end, "must span at least one byte!"

    infile.seek(start, os.SEEK_SET)

    # read the bgzip header
    # see https://samtools.github.io/hts-specs/SAMv1.pdf

    # use a buffer so we can look at the whole header, but
    # step through the file
    buffer = b""

    while start < end:

        if len(buffer) < headersize:
            bytesread = infile.read(headersize - len(buffer))
            buffer = buffer + bytesread
        if len(buffer) < headersize:
            logger.warning(f"Unable to read up to {headersize} at {start}")
            break

        header = struct.unpack(headerpattern, buffer)
        # this is a valid location for a block
        if check_is_header(header):
            logger.debug(f"Matched header at {start}")
            buffer = b""

            cdata = get_cdata(header, infile)

            # need to handle the tail
            tail = get_tail(infile)

            block_end = start + headersize + len(cdata) + tailsize
            logger.debug(f"Read block from {start} to {block_end}")

            yield start, header, cdata, tail

            # prepare for the next block
            start = block_end

        else:
            # move ahead a byte
            buffer = buffer[1:]
            start += 1


def blocks_decompressed_in_file_slice(infile, start, end):
    for start, header, cdata, tail in blocks_in_file_slice(infile, start, end):
        # now do the actual decompression
        decompressor = zlib.decompressobj(
            wbits=-15
        )  # we've alread read the header, so ignore it
        decompressed = decompressor.decompress(cdata)

        # check it with the tail
        tail_crc, tail_isize = tail
        # check decompressed size is expected
        assert len(decompressed) == tail_isize
        # check crc check is expected
        assert zlib.crc32(decompressed) == tail_crc

        yield start, header, cdata, decompressed, tail


def get_bgzip_lines_ranged(infilename, start, end):
    with open(infilename, "rb") as infile:
        line_last_old = None if start else b""
        for (
            start,
            header,
            cdata,
            decompressed,
            tail,
        ) in blocks_decompressed_in_file_slice(infile, start, end):
            lines = decompressed.split(b"\n")
            line_first, lines, line_last = lines[0], lines[1:-1], lines[-1]

            if line_last_old is not None:
                line_first = line_last_old + line_first
                yield line_first

            for line in lines:
                yield line

            line_last_old = line_last

        # after getting all the lines in the block containing the "end" location
        # we also need to get the first partial line in the next block too
        block_peek_start = infile.tell()
        block = tuple(
            itertools.islice(
                blocks_decompressed_in_file_slice(
                    infile, block_peek_start, block_peek_start + 1
                ),
                1,
            )
        )

        if block:
            start, header, cdata, decompressed, tail = block[0]
            yield line_last_old + decompressed.split(b"\n", 1)[0]
        elif line_last_old:
            # normally the last block in the file is empty
            yield line_last_old


def _get_lines_to_pipe(
    infilename, start, end, preproc, preprocargs, preprockwargs, q, linebatch
):
    lines = []
    for line in get_bgzip_lines_ranged(infilename, start, end):
        # turn bytes into string, easier for preproc to handle
        line = line.decode()
        # optionally pre-process the line *within* the subprocess
        # this can filter & trim in parallel
        if preproc:
            line = preproc(line, *preprocargs, **preprockwargs)
        # only pass back "real" lines
        if line:
            lines.append(line)
            # batch is full, send it and start a new one
            if len(lines) >= linebatch:
                q.put(lines)
                lines = []
    # send any leftover lines smaller than a batch
    q.put(lines)
    # send a sentinel
    q.put(None)


def get_bgzip_lines_parallel(
    infilename,
    nprocs=multiprocessing.cpu_count(),
    linebatch=2048,  # number of lines in a batch
    queuesize=2,  # number of batches queued per subprocess
    preproc=None,
    preprocargs=[],
    preprockwargs={},
):
    # this is the number of bytes in each chunk of the file
    # each chunk will be handled by a separate child process
    procsize = os.stat(infilename).st_size / nprocs

    # this is the queue the lines will be returned from subprocs via
    q = multiprocessing.Queue(nprocs * queuesize)

    # create the subprocesses
    subprocs = []
    for i in range(nprocs):
        start = int(i * procsize)
        end = int((i + 1) * procsize - 1)
        subproc = multiprocessing.Process(
            target=_get_lines_to_pipe,
            args=(
                infilename,
                start,
                end,
                preproc,
                preprocargs,
                preprockwargs,
                q,
                linebatch,
            ),
        )
        subprocs.append(subproc)

    # start all the subprocesses
    for subproc in subprocs:
        subproc.start()

    # to keep track of the number of active subprocs
    active = nprocs
    while active > 0:
        lines = q.get()
        # we recieved a sentinel value to say that this chunk is complete
        if lines is None:
            active -= 1
        else:
            for line in lines:
                yield line

    # make sure all the subprocesses terminate tidily
    # at this point we have recieved all the sentinels, so we should be fine
    for subproc in subprocs:
        subproc.join(1)


class BlockGZipWriter(io.BufferedWriter):
    # buffer size is 64kb which is 1 block
    def __init__(self, raw, buffer_size=65536):
        super().__init__(raw, buffer_size)

    def _flush_unlocked(self):
        if self.closed:
            raise ValueError("flush on closed file")
        while self._write_buf:
            block = self.make_block(self._write_buf)
            try:
                n = self.raw.write(block)
            except BlockingIOError:
                raise RuntimeError(
                    "self.raw should implement RawIOBase: it "
                    "should not raise BlockingIOError"
                )
            if n is None:
                raise BlockingIOError(
                    errno.EAGAIN, "write could not complete without blocking", 0
                )
            if n > len(block) or n < 0:
                raise OSError("write() returned incorrect number of bytes")
            del self._write_buf[:n]

    @staticmethod
    def compress_content(content: bytes):
        # make a new compressor each time
        compressor = zlib.compressobj(wbits=-15)
        compressed = compressor.compress(content)
        compressed = compressed + compressor.flush()

        return compressed

    @staticmethod
    def generate_header(compressed: bytes) -> bytes:
        patternhead = "<BBBBIBBHBBHH"
        header = [0] * 12
        header[0] = 31  # ID1
        header[1] = 139  # ID2
        header[2] = 8  # compression method
        header[3] = 4  # flags bit2 FEXTRA
        header[4] = 0  # MTIME
        header[5] = 0  # eXtra FLags
        header[6] = 255  # OS 255 is default unspecified
        header[7] = 6  # XLEN
        header[8] = 66
        header[9] = 67
        header[10] = 2
        header[11] = len(compressed) + 6 + 19
        headerbytes = struct.pack(patternhead, *header)
        return headerbytes

    @staticmethod
    def generate_tail(content: bytes) -> bytes:
        patterntail = "<II"
        tail_crc = zlib.crc32(content)
        tail_isize = len(content)
        tail = [tail_crc, tail_isize]
        tailbytes = struct.pack(patterntail, *tail)
        return tailbytes

    @classmethod
    def make_block(cls, content: bytes) -> bytes:
        compressed = cls.compress_content(content)
        headerbytes = cls.generate_header(compressed)
        tailbytes = cls.generate_tail(content)

        return headerbytes + compressed + tailbytes
