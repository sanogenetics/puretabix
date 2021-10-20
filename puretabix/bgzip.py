import io
import itertools
import logging
import os
import struct
import zlib

from .mp import from_multiprocess_generator

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


def check_is_gzip(fileobj):
    """is bytes readable file a valid gzip file?"""
    fileobj.seek(0)
    bytes_data = fileobj.read(3)
    header = struct.unpack("<BBB", bytes_data)
    return header[0] == 31 and header[1] == 139 and header[2] == 8


def check_is_block_gzip(fileobj):
    """is bytes readable file is a valid block gzip file?"""
    if not check_is_gzip(fileobj):
        return False
    # NOTE assumes there is only one extra header
    # not sure if this is required by block gzip spec or not
    fileobj.seek(12)
    bytes_data = fileobj.read(4)
    header = struct.unpack("<ccH", bytes_data)
    return header[0] == b"B" and header[1] == b"C" and header[2] == 2


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


def get_header(infile: io.IOBase):
    bytesread = infile.read(headersize)
    if len(bytesread) < headersize:
        raise EOFError
    header = struct.unpack(headerpattern, bytesread)
    assert check_is_header(header)
    return header


def get_cdata(header, infile: io.IOBase):
    blocksize = header[11] - header[7] - 19
    cdata = infile.read(blocksize)
    assert len(cdata) == blocksize, f"Unable to read up to {blocksize} of cdata"
    return cdata


def get_tail(infile: io.IOBase, decompressed=None):
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


def get_block(infile):
    header = get_header(infile)
    cdata, decompressed = get_cdata_decompressed(header, infile)
    tail = get_tail(infile, decompressed)
    return header, cdata, decompressed, tail


def get_lines(infile):
    start = infile.tell()
    header, cdata, decompressed, tail = get_block(infile)
    decompressedlines = decompressed.split(b"\n")
    if start == 0:
        # first block has no partial start line
        firstline = b""
        lines = decompressedlines[0:-1]
        lastline = decompressedlines[-1]
    else:
        firstline = decompressedlines[0]
        lines = decompressedlines[1:-1]
        lastline = decompressedlines[-1]
    return firstline, lines, lastline


def generate_offset_lines(infile):
    # go to the start
    infile.seek(0)
    start_previous = None
    start = 0
    lastline = ""
    offset = 0
    while True:
        start_previous = start
        start = infile.tell()
        header, cdata, decompressed, tail = get_block(infile)
        if not decompressed:
            # empty block is last block
            if lastline:
                yield start_previous, offset, start_previous, offset + len(
                    lastline
                ), lastline

        decompressedlines = decompressed.split(b"\n")
        if start == 0:
            # first block has no partial start line
            firstline = b""
            offset = 0
            lines = decompressedlines[0:-1]
            lastline = decompressedlines[-1]
        else:
            firstline = decompressedlines[0]
            # this should be added to the last line from the previous block
            yield start_previous, offset, start, len(firstline), lastline + firstline

            offset = len(firstline) + 1
            lines = decompressedlines[1:-1]
            lastline = decompressedlines[-1]

        for line in lines:
            yield start, offset, start, offset + len(line), line
            offset += len(line) + 1


def get_file_ranged_blocks(infile, start, end):
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


def get_file_ranged_decompressed(infile, start, end):
    for start, header, cdata, tail in get_file_ranged_blocks(infile, start, end):
        # now do the actual decompression
        decompressor = zlib.decompressobj(
            wbits=-15
        )  # we've alread read the header, so ignore it
        decompressed = decompressor.decompress(cdata)

        assert not decompressor.unconsumed_tail
        assert not decompressor.unused_data

        # check it with the tail
        tail_crc, tail_isize = tail
        # check decompressed size is expected
        assert len(decompressed) == tail_isize
        # check crc check is expected
        assert zlib.crc32(decompressed) == tail_crc

        yield start, header, cdata, decompressed, tail


def get_file_ranged_lines(infile, start, end):
    line_last_old = None if start else b""
    for (
        start,
        header,
        cdata,
        decompressed,
        tail,
    ) in get_file_ranged_decompressed(infile, start, end):
        lines = decompressed.split(b"\n")
        line_first, lines, line_last = lines[0], lines[1:-1], lines[-1]

        if line_last_old is not None:
            line_first = line_last_old + line_first
            yield line_first.decode("utf-8")

        for line in lines:
            yield line.decode("utf-8")

        line_last_old = line_last

    # after getting all the lines in the block containing the "end" location
    # we also need to get the first partial line in the next block too
    block_peek_start = infile.tell()
    block = tuple(
        itertools.islice(
            get_file_ranged_decompressed(
                infile, block_peek_start, block_peek_start + 1
            ),
            1,
        )
    )

    if block:
        start, header, cdata, decompressed, tail = block[0]
        line_last_old = line_last_old + decompressed.split(b"\n", 1)[0]

    if line_last_old:
        yield line_last_old.decode("utf-8")


def get_filename_ranged_lines(infilename, start, end):
    with open(infilename, "rb") as infile:
        for line in get_file_ranged_lines(infile, start, end):
            yield line


def get_filename_parallel_lines(infilename, nprocs=multiprocessing.cpu_count()):
    # this is the number of bytes in each chunk of the file
    # each chunk will be handled by a separate child process
    procsize = os.stat(infilename).st_size / nprocs
    genfunckwargs = []
    for i in range(nprocs):
        start = int(i * procsize)
        end = int((i + 1) * procsize - 1)
        genfunckwargs.append({"infilename": infilename, "start": start, "end": end})

    # dealing with individual lines, so use large batches
    for _, item in from_multiprocess_generator(
        get_filename_ranged_lines, genfunckwargs, batchsize=1024
    ):
        # don't yield empty lines
        if item:
            yield item


class BlockGZipWriter(io.BufferedIOBase):
    # buffer size is 64kb which is 1 block
    # 65280 is what bgzip uses, for some reason?
    def __init__(self, raw: io.IOBase, block_size=65536 - 256):
        self.raw = raw
        assert self.raw.writable()
        self.block_size = block_size
        self.block_buffer = b""

    def write(self, b):
        self.block_buffer = self.block_buffer + b
        while len(self.block_buffer) > self.block_size:
            content = self.block_buffer[: self.block_size]
            self.block_buffer = self.block_buffer[len(content) :]
            block = self.make_block(content)
            self.raw.write(block)

    def flush(self):
        while len(self.block_buffer):
            content = self.block_buffer[: self.block_size]
            self.block_buffer = self.block_buffer[len(content) :]
            block = self.make_block(content)
            self.raw.write(block)

        self.raw.flush()

    def close(self):
        self.flush()
        # add an empty block at the end
        self.raw.write(self.make_block(b""))
        self.raw.close()

    @staticmethod
    def compress_content(content: bytes):
        # make a new compressor each time
        compressor = zlib.compressobj(wbits=-15)
        compressed = compressor.compress(content)
        compressed = compressed + compressor.flush()

        return compressed

    @staticmethod
    def generate_header(compressed: bytes) -> bytes:
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
        headerbytes = struct.pack(headerpattern, *header)
        return headerbytes

    @staticmethod
    def generate_tail(content: bytes) -> bytes:
        tail_crc = zlib.crc32(content)
        tail_isize = len(content)
        tail = [tail_crc, tail_isize]
        tailbytes = struct.pack(tailpattern, *tail)
        return tailbytes

    @classmethod
    def make_block(cls, content: bytes) -> bytes:
        compressed = cls.compress_content(content)
        headerbytes = cls.generate_header(compressed)
        tailbytes = cls.generate_tail(content)

        return headerbytes + compressed + tailbytes
