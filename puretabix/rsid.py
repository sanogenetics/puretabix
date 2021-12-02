import gzip
import logging
import os
import struct
from types import FunctionType
from typing import Any, Dict, Iterable, List, Tuple

from .bgzip import BlockGZipReader
from .mp import MultiprocessGeneratorPool

logger = logging.getLogger(__name__)


def _get_block_rsid_line(
    infile: BlockGZipReader, block_rsids: Iterable[Tuple[int, Iterable[int]]]
):
    lineloc = None
    lines = []
    lastline = None

    for blockstart, rsids in block_rsids:
        if not lineloc or lineloc != blockstart:
            pos = infile.tell()
            if pos != blockstart:
                infile.seek(blockstart)
            _, _, _, _, lines, lastline, _ = infile.get_block_lines()
            lineloc = blockstart

        linemapping = {}
        for line in lines:
            # skip comments, blanks
            if line.startswith(b"#") or not line.strip():
                continue
            _, _, linersid, _ = line.split(b"\t", 3)
            assert linersid[0:2] == b"rs"
            linersidnum = int(linersid[2:])
            linemapping[linersidnum] = line

        # work out which rsids we have not found yet
        missing = []
        for rsid in rsids:
            # rsids is a sorted set
            if rsid in linemapping:
                yield rsid, linemapping[rsid]
                missing = []  # clear the last missing list because we found one
            else:
                missing.append(rsid)

        # any missing rsids at the end of the block _might_ be split into the next block
        if missing:
            # read the next block, and track where we read it from
            nextblockstart = infile.tell()
            (
                _,
                _,
                _,
                nextfirstline,
                nextlines,
                nextlastline,
                _,
            ) = infile.get_block_lines()

            # combine the bits of the blocks to make the line that spans the gap
            line = lastline + nextfirstline
            # skip comments, blanks
            if line.startswith(b"#") or not line.strip():
                continue
            _, _, linersid, _ = line.split(b"\t", 3)
            assert linersid[0:2] == b"rs"
            linersidnum = int(linersid[2:])

            # return it if this is one of the missing
            if linersidnum in missing:
                yield linersidnum, line

            # record that we moved on to the next block so we don't have to seek back and read it again
            lines = nextlines
            lastline = nextlastline
            lineloc = nextblockstart


def _pool_func(
    vcf_f,
    vcf_f_args,
    block_rsids: Iterable[Tuple[int, Iterable[int]]],
    post_process: FunctionType,
):
    with vcf_f(*vcf_f_args) as infile:
        bzipreader = BlockGZipReader(infile)
        for rsid, line in _get_block_rsid_line(bzipreader, block_rsids):
            # now processes it as desired
            result = post_process(line)
            yield rsid, result


class RsidIndex:
    rsidnums: Tuple[int, ...]
    fileposes: Tuple[int, ...]

    def __init__(self, rsidnums: Iterable[int], fileposes: Iterable[int]):
        self.rsidnums = tuple(rsidnums)
        self.fileposes = tuple(fileposes)

    @classmethod
    def from_file(cls, _file):
        # check if it is gzipped
        # peak at the first 3 bytes, then reset
        _file.seek(0)
        gzip_peek = _file.read(3)
        _file.seek(0)
        # do the actual check
        header = struct.unpack("<BBB", gzip_peek)
        gzipped = header[0] == 31 and header[1] == 139 and header[2] == 8

        if gzipped:
            # wrap the file in a decompressor
            # typical compression is ~50%
            index = cls.from_file(gzip.open(_file))
        else:
            rsidnums = []
            fileposes = []
            indexpattern = "<QI"
            bytecount = struct.calcsize(indexpattern)
            while True:
                recordbytes = _file.read(bytecount)
                if len(recordbytes) != bytecount:
                    break
                filepos, rsidnum = struct.unpack(indexpattern, recordbytes)
                rsidnums.append(rsidnum)
                fileposes.append(filepos)
            index = cls(rsidnums, fileposes)
        return index

    @classmethod
    def from_filename(cls, filename):
        with open(filename, "rb") as indexfile:
            index = cls.from_file(indexfile)
        return index

    @classmethod
    def from_vcf_filename(cls, vcffilename):
        with open(vcffilename, mode="rb") as vcffile:
            return cls.from_vcf(vcffile, os.stat(vcffilename).st_size)

    @classmethod
    def from_vcf(cls, vcfile, filesize):
        bgzipin = BlockGZipReader(vcfile)
        bgzipin.seek(0)
        rsidnums = []
        fileposes = []
        # read through the vcf blocks
        while True:
            (
                block_start,
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
            ) = bgzipin.scan_block_lines_offset()

            # stop on an empty block
            if not decompressed:
                break
            # get the first non-comment line in each block
            line_rests = decompressed
            while True:
                # first line in the file is the start of the first line
                if block_start == 0:
                    line, line_rests = line_rests.split(b"\n", 1)
                else:
                    assert line_rests.count(b"\n") > 1, line_rests
                    _, line, line_rests = line_rests.split(b"\n", 2)

                # ignore comments, move to next line
                if not line.startswith(b"#"):
                    break

            # get the rsid at the start of each block
            _, _, rsid, _ = line.split(b"\t", 3)
            rsidnum = int(rsid[2:])
            rsidnums.append(rsidnum)
            fileposes.append(block_start)

        index = cls(rsidnums, fileposes)
        return index

    def write_to(self, outfile):
        """
        Writes the current index to a provided file-like binary object
        Also compresses with GZip on the way
        """
        with gzip.GzipFile(fileobj=outfile) as outfilegzip:
            for offset, rsidnum in zip(self.fileposes, self.rsidnums):
                # write blockoffset, rsidnum to file
                # needs to be 64 bit for rsids over 4294967295
                outfilegzip.write(struct.pack("<QI", offset, rsidnum))

    def get_block_rsids(self, rsids: Iterable[int]) -> Dict[int, List[int]]:
        """
        This creates a dicrionary of the start of a block to
        the rsids are interested in inside that block.
        Thus we can ignore the other blocks.
        """

        indexpointer = 0
        blocksrsids = {}
        for rsid in sorted(rsids):
            # while there is a next block in the index
            # and if the rsid is beyond the current block
            # point to the next index block
            while (
                indexpointer + 1 < len(self.fileposes)
                and rsid >= self.rsidnums[indexpointer + 1]
            ):
                indexpointer += 1

            # get the file position corresponding to the index block
            filepos = self.fileposes[indexpointer]
            # get existing rsids in this index block (if any)
            thisblockrsids = blocksrsids.get(filepos, [])
            # track this rsid along any others in the same index block
            thisblockrsids.append(rsid)
            # update the mappipng with the expanded rsid list
            blocksrsids[filepos] = thisblockrsids

        return blocksrsids

    def gen_rsid_lines(
        self,
        post_process: FunctionType,
        rsids: Iterable[int],
        vcf_f: FunctionType,
        vcf_f_args: Tuple[Any],
        pool: MultiprocessGeneratorPool,
    ):
        """
        generator of (rsid, line) pairs from a VCF file
        The VCF file must be ordered by increasing numeric RSID component
        and must be block-gzip compressed
        and block indexed by numeric RSID component of first complete record

        rsids that are not in the VCF file will have "None" as the line

        index must be a RsidIndex
        infilename must be a string for the filename of the vcf file
        nproces is the number of processors to use,
        """

        # see how many distinct blocks this is
        blocksrsids = self.get_block_rsids(rsids)

        # make separate list of block->rsids for each subprocess
        # each subprocess should get blocks that are close by
        args = tuple(blocksrsids.items())
        blocksize = (len(blocksrsids) // len(pool)) + 1
        fkwargss = []
        for i in range(0, len(blocksrsids), blocksize):
            chunk_block_rsids = args[i : i + blocksize]
            fkwargss.append(
                {
                    "vcf_f": vcf_f,
                    "vcf_f_args": vcf_f_args,
                    "block_rsids": chunk_block_rsids,
                    "post_process": post_process,
                }
            )

        # send argument to pool
        pool.submit(_pool_func, fkwargss)

        # get results
        for i, (_, result) in enumerate(pool.results()):
            yield result
