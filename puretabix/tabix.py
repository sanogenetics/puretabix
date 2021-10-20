import gzip
import itertools
import logging
import struct
import zlib

from .bgzip import (
    check_is_block_gzip,
    check_is_gzip,
    generate_offset_lines,
    get_file_ranged_decompressed,
)
from .vcf import LINE_START, VCFAccumulator, get_vcf_fsm

logger = logging.getLogger(__name__)


class TabixIndex:
    def __init__(
        self,
        file_format,
        column_sequence,
        column_begin,
        column_end,
        meta,
        headerlines_count,
        indexes,
    ):
        """
        In-memory representation of a Tabix index. See https://samtools.github.io/hts-specs/tabix.pdf
        for more information.

        Generally these are pretty small files that need to be read entirely, thus downloading them
        locally before processing is recommented e.g. io.BytesIO
        """

        # 0 = generic tab-delemited
        # 1 = SAM
        # 2 = VCF
        self.file_format = file_format
        # column for sequence ids, 1-based
        self.column_sequence = column_sequence
        # column for region start, 1-based
        self.column_begin = column_begin
        # column for region end, 1-based
        self.column_end = column_end
        self.meta = meta
        self.headerlines_count = headerlines_count

        # a dictionary of names to (bin_index, interval_index)
        self.indexes = indexes

    @classmethod
    def from_file(cls, fileobj):
        """
        Generally these are pretty small files that need to be read entirely, thus downloading them
        locally before processing is recommented e.g. io.BytesIO
        """
        # the index file is block-gzipped but small enough we can
        # load it into memory and process like a regular gzip file
        with gzip.GzipFile(fileobj=fileobj) as f:
            header_pattern = "<4siiiii4sii"
            header = struct.unpack(
                header_pattern, f.read(struct.calcsize(header_pattern))
            )

            magic = header[0]
            if magic != b"TBI\01":  # check magic
                raise RuntimeError(f"invalid tabix index magic {magic}.")

            # number of named sequences (e.g. chromosomes)
            n_sequences = header[1]

            file_format = header[2]
            # 0 = generic tab-delemited
            # 1 = SAM
            # 2 = VCF
            if file_format not in (0, 1, 2):
                raise RuntimeError(f"invalid tabix index format {file_format}.")

            # these are 1 based
            # value of 0 states not included in file
            # e.g. VCF has no explicit end column
            column_sequence = header[3]  # Column for the sequence name
            column_begin = header[4]  # Column for the start of a region
            column_end = header[5]  # Column for the end of a region

            # this is the comment marker, usually #
            meta = header[6].decode("ascii")[0]
            assert meta == "#", (header[6], meta)

            # number of lines of header at the start of the file
            # this does not include lines marked as comments
            headerlines_count = header[7]

            # sequence names are a series of bytes followed by a null byte
            names = tuple(
                map(bytes.decode, f.read(header[8]).split(b"\x00")[:-1])
            )  # throw the last empty one away
            if len(names) != n_sequences:
                raise RuntimeError(
                    f"unexpected number of sequences {n_sequences} vs {len(names)}"
                )

            indexes = {}
            # for each sequence
            for name in names:
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

                # parse the interval index
                n_intervals = struct.unpack("<i", f.read(4))[0]
                intervals = [
                    i[0] for i in struct.iter_unpack("<Q", f.read(8 * n_intervals))
                ]

                if name in indexes:
                    raise RuntimeError(f"duplicate sequence name {name}")
                indexes[name] = (bins, intervals)

        return cls(
            file_format,
            column_sequence,
            column_begin,
            column_end,
            meta,
            headerlines_count,
            indexes,
        )

    def __repr__(self):
        return f"TabixIndex({self.file_format}, {self.column_sequence}, {self.column_begin}, {self.column_end}, {self.meta}, {self.headerlines_count}, {self.indexes})"

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
        if sequence_name not in self.indexes:
            return None
        linear_index = self.indexes[sequence_name][1]
        # if it would be beyond the index, say that
        if i >= len(linear_index):
            return None
        # its a valid sequnce name and a valid interval window
        return linear_index[i]

    def _lookup_bin_chunks(self, sequence_name, start, end):
        """
        Records are assigned to a bin if the entirely fit in the bin.
        So we want all the records in all the bins that overlap with the region of interest.
        These records *might* overlap with the region of interest.
        """
        bin_index = self.indexes[sequence_name][0]
        for chunks_bin_index in reversed(tuple(self.region_to_bins(start, end))):
            if chunks_bin_index in bin_index:
                for chunk in bin_index[chunks_bin_index]:
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

    def write_to(self, outfile):
        # header
        outfile.write(
            struct.pack(
                "<4siiiii4si",
                b"TBI\01",  # magic number
                len(self.indexes),  # n sequences
                self.file_format,  # file format 0 generic, 1 sam, 2 vcf
                self.column_sequence,  # column for sequence ids, 1-based
                self.column_begin,  # column for region start, 1-based
                self.column_end,  # column for region end, 1-based
                self.meta.encode("ascii")
                + b"\x00\x00\x00",  # this is a character, but represented as a int
                self.headerlines_count,
            )
        )

        # length of concatenated zero terminated names
        names = tuple(self.indexes.keys())  # ensure consistent order
        names_concat = b"".join((i.encode("ascii") + b"\0" for i in names))
        outfile.write(
            struct.pack(f"<i{len(names_concat)}s", len(names_concat), names_concat)
        )

        for name in names:
            # n_bin
            #   bin
            #   n_chunk
            #     chunk_begin
            #     chunk_end
            # n_intv
            #   ioff
            bin_index = self.indexes[name][0]
            outfile.write(struct.pack("<i", len(bin_index)))
            for bin_i in bin_index.keys():
                # unsigned bin number, n_chunk
                outfile.write(struct.pack("<Ii", bin_i, len(bin_index[bin_i])))
                for chunk_begin, chunk_end in bin_index[bin_i]:
                    outfile.write(struct.pack("<QQ", chunk_begin, chunk_end))

            intv_index = self.indexes[name][1]
            outfile.write(struct.pack("<i", len(intv_index)))
            for ioff in intv_index:
                outfile.write(struct.pack("<Q", ioff))

    @classmethod
    def build_from(cls, vcffile):
        """
        read a vcf file in blockgzip format and create an index object for it
        """

        file_format = 2  # 0 generic, 1 sam, 2 vcf
        column_sequence = 1  # column for sequence ids, 1-based
        column_begin = 2  # column for region start, 1-based
        column_end = 0  # column for region end, 1-based
        meta = "#"

        # dictionary of names to (bin_index, interval_index)
        # bin_index is a dictionary of bin numbers to [(chunk_start, chunk_end),...]
        # interval_index is a list of 16kbase interval start locations
        indexes = {}

        # go through each line in turn
        # parse it into a structured object
        vcf_fsm = get_vcf_fsm()
        accumulator = VCFAccumulator()

        # these are the internal two index types
        bin_index = {}
        interval_index = {}

        try:
            for (
                start_block,
                start_offset,
                end_block,
                end_offset,
                line,
            ) in generate_offset_lines(vcffile):
                if not line:
                    continue

                vcf_fsm.run(line.decode() + "\n", LINE_START, accumulator)
                vcf_line = accumulator.to_vcfline()
                accumulator.reset()
                # TODO add better debugging for unexpected lines

                # skip any comment lines with # or ##
                if vcf_line.comment_raw or vcf_line.comment_key:
                    continue

                # have we started a new chromosome?
                if vcf_line.chrom not in indexes:
                    indexes[vcf_line.chrom] = ({}, [])
                    bin_index = indexes[vcf_line.chrom][0]
                    interval_index = indexes[vcf_line.chrom][1]

                # get the combined number for the block & offset
                start_virtual = start_block << 16 | start_offset
                end_virtual = end_block << 16 | end_offset

                # subtract 1 because its 0 offset
                record_start = vcf_line.pos - 1
                # subtract another 1 because half-open end
                record_end = vcf_line.pos - 1 + len(vcf_line.ref) - 1

                # bin index
                # smallest bin that completely contains the record

                bin_i = cls.region_to_bin(record_start, record_end)

                if bin_i not in bin_index.keys():
                    bin_index[bin_i] = [(start_virtual, end_virtual + 1)]
                else:
                    # extend chunk if directly continuous
                    if bin_index[bin_i][-1][1] == start_virtual:
                        bin_index[bin_i][-1] = (
                            bin_index[bin_i][-1][0],
                            end_virtual + 1,
                        )
                    else:
                        bin_index[bin_i].append((start_virtual, end_virtual + 1))

                # interval index
                # is the lowest virtual offset of all records that overlap interval
                # half-closed half-open records

                # throw away the first 14 bits to get interval index
                # e.g. 0 => 0, 16384 => 1, etc
                interval_i_start = record_start >> 14
                # subtract another 1 because half-open end
                interval_i_end = record_end >> 14

                for linear_i in sorted(set((interval_i_start, interval_i_end))):
                    # add any padding intervals
                    while len(interval_index) <= linear_i:
                        interval_index.append(-1)

                    # if there is no previous value, or the previous value was higher
                    # replace previous value
                    if (
                        interval_index[linear_i] == -1
                        or start_virtual < interval_index[linear_i]
                    ):
                        interval_index[linear_i] = start_virtual

                # pass through again and replace any -1 values with the last "not -1" value
                # these are intervals with no overlapping records
                lastgood = -1
                for i in range(len(interval_index)):
                    if interval_index[i] == -1:
                        interval_index[i] = lastgood
                    else:
                        lastgood = interval_index[i]
        except EOFError:
            # if we reach the end of the file, ignore it because index will be complete
            # TODO look for bgzip empty block instead
            pass

        return cls(
            file_format, column_sequence, column_begin, column_end, meta, 0, indexes
        )

    @staticmethod
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

    @staticmethod
    def region_to_bin(begin, end):
        """
        returns the index of the smallest bin that contains the region
        as a half-closed half-open interval
        """
        if begin >> 14 == end >> 14:
            return ((1 << 15) - 1) // 7 + (begin >> 14)
        elif begin >> 17 == end >> 17:
            return ((1 << 12) - 1) // 7 + (begin >> 17)
        elif begin >> 20 == end >> 20:
            return ((1 << 9) - 1) // 7 + (begin >> 20)
        elif begin >> 23 == end >> 23:
            return ((1 << 6) - 1) // 7 + (begin >> 23)
        elif begin >> 26 == end >> 26:
            return ((1 << 3) - 1) // 7 + (begin >> 26)
        return 0

    @staticmethod
    def bin_start(k):
        # bit_length-1 is int log2
        lvl = (((7 * k) + 1).bit_length() - 1) // 3
        ol = ((2 ** (3 * lvl)) - 1) // 7
        sl = 2 ** (29 - (3 * lvl))
        start = (k - ol) * sl
        return start

    @staticmethod
    def bin_size(k):
        # bit_length-1 is int log2
        lvl = (((7 * k) + 1).bit_length() - 1) // 3
        sl = 2 ** (29 - (3 * lvl))
        return sl


class TabixIndexedFile:
    def __init__(self, fileobj, index):
        self.index = index
        self.fileobj = fileobj
        # check fileobject is a real block-gzipped file
        if not self.is_block_gzip():
            raise ValueError("fileobj must be a block gzipped file-like object")

    @staticmethod
    def from_files(fileobj, index_fileobj):
        return TabixIndexedFile(fileobj, TabixIndex.from_file(index_fileobj))

    def is_gzip(self):
        """is fileobj a valid gzip file?"""
        return check_is_gzip(self.fileobj)

    def is_block_gzip(self):
        """is fileobj is a valid block gzip file?"""
        return check_is_block_gzip(self.fileobj)

    def fetch_bytes_virtual(self, virtual_start, virtual_end):
        # the lower 16 bits store the offset of the byte inside the gzip block
        # the rest store the offset of gzip block
        block_start = virtual_start >> 16
        offset_start = virtual_start & 0xFFFF
        block_end = virtual_end >> 16
        offset_end = virtual_end & 0xFFFF
        return self.fetch_bytes_block_offset(
            block_start, offset_start, block_end, offset_end
        )

    def fetch_bytes_block_offset(
        self, block_start, offset_start, block_end, offset_end
    ):
        value = b""
        for start, _, _, decompressed, _ in get_file_ranged_decompressed(
            self.fileobj, block_start, block_end + 1
        ):
            if start == block_end:
                decompressed = decompressed[:offset_end]
            if start == block_start:
                decompressed = decompressed[offset_start:]
            value = value + decompressed
        return value

    def fetch_bytes(self, name, start, end=None):
        """
        Returns bytes in the region of interest
        """
        # quick check
        if name not in self.index.indexes.keys():
            return b""

        # default if only start specified
        if not end:
            end = start

        # use the index to get "virtual" file offsets that include all of the region of interest
        virtual_start, virtual_end = self.index.lookup_virtual(name, start, end)

        # location not indexed, return empty string
        if not virtual_start and not virtual_end:
            return b""

        return self.fetch_bytes_virtual(virtual_start, virtual_end)

    def fetch(self, name, start, end=None):
        """
        Returns lines in the region of interest
        """
        # default if only start specified
        if not end:
            end = start

        region = self.fetch_bytes(name, start, end)
        # no region, no lines
        if not region:
            return ""

        # turn the bytes into a list of strings
        # TODO ascii vs utf-8?
        lines = region.decode("ascii").splitlines()

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
        if offset is not None:
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
