import gzip
import io
import itertools
import struct
import zlib


def region_to_bins(begin, end, n_levels=5, min_shift=14):
    """
    generator of keys to bins of records which *may* overlap the given region

    n_levels: int, optional
        cluster level, 5 for tabix
    min_shift: int, optional
        minimum shift, 14 for tabix
    """
    t = 0
    s = min_shift + (n_levels << 1) + n_levels
    for level in range(n_levels + 1):
        b = t + (begin >> s)
        e = t + (end >> s)
        n = e - b + 1
        for k in range(b, e + 1):
            yield k
            n += 1
        t += 1 << ((level << 1) + level)
        s -= 3


class TabixIndex:
    def __init__(self, fileobj):
        """
        In-memory representation of a Tabix index. See https://samtools.github.io/hts-specs/tabix.pdf
        for more information.

        Generally these are pretty small files that need to be read entirely, thus downloading them
        locally before processing is recommented e.g. io.BytesIO
        """
        self._fileobject = fileobj
        self.magic = None
        self.n_sequences = None
        self.file_format = None
        self.column_sequence = 0  # column for sequence ids, 1-based
        self.column_begin = 0  # column for region start, 1-based
        self.column_end = 0  # column for region end, 1-based
        self.meta = None
        self.headerlines_count = None
        self.names = None
        self.index_bin = {}
        self.index_interval = {}

        # pre-process the file
        self._parse_index()

    def _parse_index(self):
        # the index file is block-gzipped but small enough we can
        # load it into memory and process like a regular gzip file
        with gzip.GzipFile(fileobj=self._fileobject) as f:
            header_pattern = "<4siiiii4sii"
            header = struct.unpack(
                header_pattern, f.read(struct.calcsize(header_pattern))
            )

            self.magic = header[0]
            if self.magic != b"TBI\01":  # check magic
                raise RuntimeError(f"invalid tabix index magic {self.magic}.")

            self.file_format = header[2]
            # 0 = generic tab-delemited
            # 1 = SAM
            # 2 = VCF
            if self.file_format not in (0, 1, 2):
                raise RuntimeError(f"invalid tabix index format {self.file_format}.")

            # these are 1 based
            # value of 0 states not included in file
            # e.g. VCF has no explicit end column
            self.column_sequence = header[3]  # Column for the sequence name
            self.column_begin = header[4]  # Column for the start of a region
            self.column_end = header[5]  # Column for the end of a region

            # this is the comment marker, usually #
            self.meta = header[6].decode("ascii")

            # number of lines of header at the start of the file
            # this does not include lines marked as comments
            self.headerlines_count = header[7]

            # sequence names are a series of bytes followed by a null byte
            self.names = tuple(
                map(bytes.decode, f.read(header[8]).split(b"\x00")[:-1])
            )  # throw the last empty one away
            if len(self.names) != header[1]:
                raise RuntimeError(
                    f"unexpected number of sequences {header[1]} vs {len(self.names)}"
                )

            # for each sequence
            for name in self.names:
                # each sequence has a bin index and an interval index

                # parse the bin index
                n_bins = struct.unpack("<i", f.read(4))[0]
                bins = {}
                for _ in range(n_bins):
                    # each bin has a key, and a series of chunks
                    bin_key, n_chunks = struct.unpack("<Ii", f.read(8))
                    chunks = [
                        chunk
                        for chunk in struct.iter_unpack("<QQ", f.read(16 * n_chunks))
                    ]

                    assert bin_key not in bins
                    bins[bin_key] = chunks
                if name in self.index_bin:
                    raise RuntimeError(f"duplicate sequence name {name}")
                self.index_bin[name] = bins

                # parse the interval index
                n_intervals = struct.unpack("<i", f.read(4))[0]
                intervals = [
                    i[0] for i in struct.iter_unpack("<Q", f.read(8 * n_intervals))
                ]

                if name in self.index_interval:
                    raise RuntimeError(f"duplicate sequence name {name}")
                self.index_interval[name] = intervals

    def _lookup_linear(self, sequence_name, start):
        """
        For each tiling 16 kb window keep the virtual file offset of the leftmost record (i.e.
        having the smallest start coordinate) that overlaps the window.

        Given a location, get the smallest start location of all records that overlap the 16kb window
        containing the location
        """
        # linear index is in 16kb intervals
        # 16kb = 16 * (2 ** 10) = 2 ** 14 = 1 << 14
        # throw away the first 14 bits to get index position
        i = start >> 14
        # if this sequence_name isn't valid, say that
        if sequence_name not in self.index_interval:
            return None
        # if it would be beyond the index, say that
        if i >= len(self.index_interval[sequence_name]):
            return None
        # its a valid sequnce name and a valid interval window
        return self.index_interval[sequence_name][i]

    def _lookup_bin_chunks(self, sequence_name, start, end):
        """
        Records are assigned to a bin if the entirely fit in the bin.
        So we want all the records in all the bins that overlap with the region of interest.
        These records *might* overlap with the region of interest.
        """
        for chunks_bin_index in reversed(tuple(region_to_bins(start, end))):
            if chunks_bin_index in self.index_bin[sequence_name]:
                for chunk in self.index_bin[sequence_name][chunks_bin_index]:
                    yield chunk

    def lookup_virtual(self, sequence_name, start, end):
        virtual_start = None
        virtual_end = None

        linear_start = self._lookup_linear(sequence_name, start)
        # if this is not in the linear index, cant return anything
        if not linear_start:
            return None, None

        for chunk_start, chunk_end in self._lookup_bin_chunks(
            sequence_name, start, end
        ):
            if chunk_end <= linear_start:
                # if the chunk finished before this section of the linear starts, skip the chunk
                # rare, but does happen sometimes
                continue

            # move the chunk start to where the linear start begins
            chunk_start = min(chunk_start, linear_start)

            if virtual_start is None or chunk_start < virtual_start:
                virtual_start = chunk_start

            if virtual_end is None or chunk_end > virtual_end:
                virtual_end = chunk_end

        # either both or neither must be set
        assert (virtual_start is None) == (virtual_end is None)

        return virtual_start, virtual_end


class TabixIndexedFile:
    def __init__(self, fileobj, index_fileobj):
        self.index = TabixIndex(index_fileobj)
        self.fileobj = fileobj
        # check fileobject is a real block-gzipped file
        if not self.is_block_gzip():
            raise ValueError("fileobj must be a block gzipped file-like object")

    def is_gzip(self):
        """is bytes_data a valid gzip file?"""
        self.fileobj.seek(0)
        bytes_data = self.fileobj.read(3)
        header = struct.unpack("<BBB", bytes_data)
        return header[0] == 31 and header[1] == 139 and header[2] == 8

    def is_block_gzip(self):
        """is bytes_data is a valid block gzip file?"""
        if not self.is_gzip():
            return False
        # NOTE assumes there is only one extra header
        # not sure if this is required by block gzip spec or not
        self.fileobj.seek(12)
        bytes_data = self.fileobj.read(4)
        header = struct.unpack("<ccH", bytes_data)
        return header[0] == b"B" and header[1] == b"C" and header[2] == 2

    def fetch(self, name, start, end=None):
        """
        Returns a text block of lines that are included in the region of interest
        """
        # quick check
        if name not in self.index.names:
            return ""

        # default if only start specified
        if not end:
            end = start

        # use the index to get "virtual" file offsets that include all of the region of interest
        virtual_start, virtual_end = self.index.lookup_virtual(name, start, end)

        # location not indexed, return empty string
        if not virtual_start and not virtual_end:
            return ""

        # the lower 16 bits store the offset of the byte inside the gzip block
        # the rest store the offset of gzip block
        block_start = virtual_start >> 16
        offset_start = virtual_start & 0xFFFF
        block_end = virtual_end >> 16
        offset_end = None  # this is important only in the last block

        value = io.BytesIO()
        while block_start <= block_end:
            # if this is the last block, then we need the ending offset
            if block_start == block_end:
                offset_end = virtual_end & 0xFFFF

            # read this block
            block, block_size = self._read_block(block_start)

            # take the content of interest out of the block
            if offset_end is None:
                # we want everything else in the block
                value.write(block[offset_start:])
            else:
                # we want to stop within this block
                value.write(block[offset_start:offset_end])

            # move to next block
            offset_start = 0
            block_start += block_size

        # turn the bytes into a list of strings
        # TODO ascii vs utf-8?
        lines = io.StringIO(value.getvalue().decode("ascii"))

        # skip header lines defined in index
        if self.index.headerlines_count:
            lines = itertools.islice(lines, self.index.headerlines_count, None)

        # filter out comments
        lines = (line for line in lines if not line.startswith(self.index.meta))

        # filter lines of wrong lengths i.e. cut off around chunk boundries
        lines = (
            line
            for line in lines
            if len(line.split("\t"))
            >= max(
                self.index.column_sequence,
                self.index.column_begin,
                self.index.column_end,
            )
        )

        # filter lines before start
        lines = (
            line
            for line in lines
            if int(line.split("\t")[self.index.column_begin - 1]) >= start
        )

        # filter lines after end
        if self.index.column_end:
            # VCF doesnt need this
            raise NotImplementedError()
        else:
            lines = (
                line
                for line in lines
                if int(line.split("\t")[self.index.column_begin - 1]) <= end
            )

        # return the lines as a text block
        return "".join(lines)

    def _read_block(self, offset):
        # move the underlying file-like object to the right place
        self.fileobj.seek(offset)
        # read the first few bytes of the block from the source
        # TODO validate this block
        block_header = self.fileobj.read(18)
        size = struct.unpack("<H", block_header[16:18])[0] + 1

        # read the appropriate amount of content, excluding what has already been read
        compressed_bytes = block_header + self.fileobj.read(size - 18)
        if len(compressed_bytes) != size:
            raise RuntimeError(
                f"Unexpected number of bytes read. Got {len(compressed_bytes)}, expected {size}"
            )

        # decompress the content, which should be a complete and valid gzip file
        decompressor = zlib.decompressobj(15 + 32)
        uncompressed_bytes = decompressor.decompress(compressed_bytes)
        if decompressor.unconsumed_tail:
            raise RuntimeError(
                f"Had unconsumed tail after decompressing of {len(decompressor.unconsumed_tail)}"
            )
        if decompressor.unused_data:
            raise RuntimeError(
                f"Had unused data after decompressing of {len(decompressor.unused_data)}"
            )

        return uncompressed_bytes, size
