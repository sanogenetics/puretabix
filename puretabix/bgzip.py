import io
import logging
import struct
import zlib
from typing import Tuple

logger = logging.getLogger(__name__)

headerpattern = "<BBBBIBBHBBHH"
headersize = struct.calcsize(headerpattern)
tailpattern = "<II"
tailsize = struct.calcsize(tailpattern)


class BlockGZipReader:
    raw: io.IOBase

    def __init__(self, raw: io.IOBase):
        assert raw.seekable()
        self.raw = raw
        assert self.check_is_block_gzip()

    def seek(self, offset: int) -> int:
        return self.raw.seek(offset)

    def tell(self) -> int:
        return self.raw.tell()

    def check_is_gzip(self) -> bool:
        """
        returns a boolean for if it has a gzip header
        will seek to start of file
        """
        self.raw.seek(0)
        bytes_data = self.raw.read(3)
        header = struct.unpack("<BBB", bytes_data)
        return header[0] == 31 and header[1] == 139 and header[2] == 8

    def check_is_block_gzip(self) -> bool:
        """
        returns a boolean for if it has a block gzip header
        also checks if it has a gzip header
        will seek to start of file
        """
        if not self.check_is_gzip():
            return False
        # NOTE assumes there is only one extra header
        # not sure if this is required by block gzip spec or not
        self.raw.seek(12)
        bytes_data = self.raw.read(4)
        header = struct.unpack("<ccH", bytes_data)
        return header[0] == b"B" and header[1] == b"C" and header[2] == 2

    @staticmethod
    def check_is_header(header: Tuple) -> bool:
        """
        tests if a series of integers matches a blockgzip header
        """
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

    def get_header(self) -> Tuple:
        """
        reads the next header from the file from the current point, assuming file is currently
        at start of a block
        """
        bytesread = self.raw.read(headersize)
        if len(bytesread) < headersize:
            raise EOFError(f"Expected to read {headersize} read {len(bytesread)}")
        header = struct.unpack(headerpattern, bytesread)
        assert self.check_is_header(header)
        return header

    def get_cdata(self, header) -> bytes:
        """
        reads the compressed data of a block from the current point, given the bytes from the
        header of that block to determine size
        """
        blocksize = header[11] - header[7] - 19
        cdata = self.raw.read(blocksize)
        assert len(cdata) == blocksize, f"Unable to read up to {blocksize} of cdata"
        return cdata

    def get_tail(self, decompressed=None) -> Tuple[int, int]:
        """
        reads the tail of the block from the current point as a tuple of crc and isize
        if the decompressed bytes are provided, will use the crc checksum in the tail
        to validate the bytes match expectation
        """
        # read isize and crc check
        tailbytes = self.raw.read(tailsize)
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

    def get_cdata_decompressed(self, header) -> Tuple[bytes, bytes]:
        """
        reads the compressed data of a block from the current point, given the bytes from the
        header of that block to determine size

        decompresses the data and returns both compressed and decompressed forms
        """
        cdata = self.get_cdata(header)
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

    def get_block(self, header=None) -> Tuple[Tuple, bytes, bytes, Tuple[int, int]]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        if not header:
            header = self.get_header()
        cdata, decompressed = self.get_cdata_decompressed(header)
        tail = self.get_tail(decompressed)
        return header, cdata, decompressed, tail

    def get_block_lines(
        self, header: Tuple = None
    ) -> Tuple[Tuple, bytes, bytes, bytes, Tuple[bytes, ...], bytes, Tuple[int, int]]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        will split by line and return the partial first and last lines, as well as the complete lines in between
        note these are returned as bytes not true string objects

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        start = self.raw.tell()
        header, cdata, decompressed, tail = self.get_block(header)
        # empty block
        if len(decompressed) == 0:
            return header, cdata, decompressed, b"", (), b"", tail

        # line endings can abut block ending so keep them
        decompressedlines = decompressed.splitlines(keepends=True)
        if start == 0:
            # first block has no partial start line
            firstline = b""
            lines = decompressedlines[0:-1]
            lastline = decompressedlines[-1]
        else:
            firstline = decompressedlines[0]
            lines = decompressedlines[1:-1]
            lastline = decompressedlines[-1]
        return header, cdata, decompressed, firstline, tuple(lines), lastline, tail

    def get_block_lines_offset(
        self, header=None
    ) -> Tuple[
        Tuple,
        bytes,
        bytes,
        bytes,
        Tuple[bytes, ...],
        bytes,
        Tuple[int, ...],
        Tuple[int, ...],
        Tuple[int, int],
    ]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        will split by line and return the partial first and last lines, as well as the complete lines in between
        note these are returned as bytes not true string objects

        also returns the offsets within the block of start and end of each line (inclusive, including separator)

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        (
            header,
            cdata,
            decompressed,
            firstline,
            lines,
            lastline,
            tail,
        ) = self.get_block_lines(header)
        offsetstarts = []
        offsetends = []
        offset = len(firstline)
        for line in lines:
            offsetstarts.append(offset)
            offset += len(line)
            offsetends.append(offset)
            offset += 1
        return (
            header,
            cdata,
            decompressed,
            firstline,
            lines,
            lastline,
            tuple(offsetstarts),
            tuple(offsetends),
            tail,
        )

    def scan_block_lines_offset(self, end: int = -1):
        """
        starting from the current position, scan forward through the file for the next block start

        will read the block and leave the file pointing at the start of the next block
        """
        buffer = b""
        blockstart = self.raw.tell()

        # as long as there is more file to read
        while end < 0 or blockstart < end:
            # populate buffer
            if len(buffer) < headersize:
                bytesread = self.raw.read(headersize - len(buffer))
                buffer = buffer + bytesread
            # check not at end
            if len(buffer) < headersize:
                logger.warning(f"Unable to read up to {headersize}")
                raise EOFError()

            header = struct.unpack(headerpattern, buffer)
            # this is a valid location for a block
            if not self.check_is_header(header):
                # move ahead a byte
                buffer = buffer[1:]
                blockstart += 1
            else:
                # this is a valid location for a block
                (
                    header,
                    cdata,
                    decompressed,
                    firstline,
                    lines,
                    lastline,
                    offsetstarts,
                    offsetends,
                    tail,
                ) = self.get_block_lines_offset(header)
                blockend = self.raw.tell()
                return (
                    blockstart,
                    blockend,
                    header,
                    cdata,
                    decompressed,
                    firstline,
                    lines,
                    lastline,
                    offsetstarts,
                    offsetends,
                    tail,
                )
        # reach the end of the file without finding a block
        raise EOFError()

    def generate_lines_offset(self, end: int = -1):
        """
        starting from the current position, scan forward through the file
        generator that yields for each line in the file
        will stop at the end of the block that includes the end point, or
        at an empty block that indicates the end of the file
        """
        partialline = b""
        blockstart_previous = 0
        offsetstart_previous = 0
        more = True
        while more:
            (
                blockstart,
                _,
                _,
                _,
                decompressed,
                firstline,
                lines,
                lastline,
                offsetstarts,
                offsetends,
                _,
            ) = self.scan_block_lines_offset(end=end)

            if not decompressed:
                # empty block is end of file
                more = False
                # process the last partial line first
                blockstarts = (blockstart_previous,)
                blockends = (blockstart_previous,)
                offsetstarts = (offsetstart_previous,)
                offsetends = (offsetstart_previous + len(partialline),)
                lines = (partialline,)
            else:
                blockstarts = (blockstart_previous,) + ((blockstart,) * len(lines))
                blockends = (blockstart,) + ((blockstart,) * len(lines))
                offsetstarts = (offsetstart_previous,) + offsetstarts
                offsetends = (len(partialline),) + offsetends
                # append holdover partial line to initial line to make a new line
                lines = (partialline + firstline,) + lines
                # keep the last partial line for next block
                # last block ended on a line ending no rollover needed
                if lastline.endswith(b"\n"):
                    blockstarts = blockstarts + (blockstart,)
                    blockends = blockends + (blockstart,)
                    offsetstarts = offsetstarts + (offsetends[-1],)
                    offsetends = offsetends + (offsetends[-1] + len(lastline),)
                    lines = lines + (lastline,)
                    partialline = b""
                else:
                    partialline = lastline

            for line, start_block, start_offset, end_block, end_offset in zip(
                lines, blockstarts, offsetstarts, blockends, offsetends
            ):
                if line:
                    yield start_block, start_offset, end_block, end_offset, line

            # prepare to do next block
            blockstart_previous = blockstart
            offsetstart_previous = offsetends[-1] + 1

    def generate_lines(self, end: int = -1):
        for _, _, _, _, line in self.generate_lines_offset(end):
            yield line


class BlockGZipWriter(io.BufferedIOBase):
    # buffer size is 64kb which is 1 block
    # 65280 is what bgzip uses, for some reason?
    def __init__(self, raw: io.RawIOBase, block_size=65536 - 256):
        self.raw = raw
        assert self.raw.writable()
        self.block_size = block_size
        self.block_buffer = b""

    def write(self, data: bytes):
        self.block_buffer = self.block_buffer + data
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
