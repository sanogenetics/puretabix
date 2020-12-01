import collections
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


# jump into the middle of a file and try to scan to the first block
def check_bytes_header(header: bytes):
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


def _ranged_scan_and_read(infile, start, end, buffer=collections.deque()):
    assert start < end, "must span at least one byte!"

    infile.seek(start, os.SEEK_SET)

    # read the bgzip header
    # see https://samtools.github.io/hts-specs/SAMv1.pdf
    pattern = "<BBBBIBBHBBHH"
    patternsize = struct.calcsize(pattern)
    got_content = False
    block_start = start - len(buffer)
    while block_start < end:
        block_start_next = block_start

        if len(buffer) < patternsize:
            bytesread = infile.read(patternsize - len(buffer))
            buffer.extend(bytesread)
        if len(buffer) < patternsize:
            logger.warning(f"Unable to read up to {patternsize} at {block_start}")
            break

        headerbuffer = bytes((buffer[i] for i in range(patternsize)))

        header = struct.unpack(pattern, headerbuffer)
        logger.debug(
            f"Header at {block_start} is {header} or in raw bytes {headerbuffer}"
        )

        if check_bytes_header(header):
            logger.debug(f"Matched header at {block_start} {header}")
            # this is a valid location for a block
            blocksize = header[11] - header[7] - 19
            logger.debug(f"Block at {block_start} has cdata size = {blocksize}")
            cdata = infile.read(blocksize)
            assert (
                len(cdata) == blocksize
            ), f"Unable to read up to {blocksize} of cdata at {block_start}"
            buffer.clear()

            # now do the actual decompression
            decompressor = zlib.decompressobj(
                wbits=-15
            )  # we've alread read the header, so ignore it
            decompressed = decompressor.decompress(cdata)
            assert (
                not decompressor.unconsumed_tail
            ), f"unconsumed tail of {len(decompressor.unconsumed_tail)} at {block_start}"
            assert (
                not decompressor.unused_data
            ), f"unused data present of {len(decompressor.unused_data)} at {block_start}"

            # read isize and crc check
            tailpattern = "<II"
            tailpatternsize = struct.calcsize(tailpattern)
            tailbytes = infile.read(tailpatternsize)
            if len(tailbytes) != tailpatternsize:
                raise ValueError(
                    f"Unable to read {tailpatternsize} bytes for tail at {block_start}"
                )
            tail_crc, tail_isize = struct.unpack(tailpattern, tailbytes)
            # check decompressed size is expected
            assert len(decompressed) == tail_isize
            # check crc check is expected
            assert zlib.crc32(decompressed) == tail_crc

            # last block has no compressed content
            if decompressed:
                got_content = True
                yield decompressed
            else:
                logger.warning(f"Found empty block at {block_start}")
                break

            logger.debug(f"Read block from {block_start} to {infile.tell()}")

            block_start_next += patternsize + blocksize + tailpatternsize
        else:
            # move ahead a byte
            buffer.popleft()
            block_start_next += 1

        block_start = block_start_next

    logger.debug(f"End of scan from {start} to {end} at {block_start}-{infile.tell()}")

    if not got_content:
        logger.warning(f"Got to end of {start}-{end} without finding a block")


def get_bgzip_lines_ranged(infilename, start, end):
    with open(infilename, "rb") as infile:
        line_last_old = None if start else b""
        buffer = collections.deque()
        for decompressed in _ranged_scan_and_read(infile, start, end, buffer):
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
        # pass the buffer along to solve sometimes needing to go back a little bit in the file!
        block_peek_start = infile.tell()
        decompressed = tuple(
            itertools.islice(
                _ranged_scan_and_read(
                    infile, block_peek_start, block_peek_start + 1, buffer
                ),
                1,
            )
        )
        if len(decompressed):
            yield line_last_old + decompressed[0].split(b"\n", 1)[0]
        elif line_last_old:
            # normally the last blockin the file is empty
            logger.debug(f"Enpty block for peak into {block_peek_start}")
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
