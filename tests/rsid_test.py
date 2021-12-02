import os
import random
from io import SEEK_END
from tempfile import TemporaryDirectory

from pytest import fixture

from puretabix import BlockGZipWriter, RsidIndex, VCFLine
from puretabix.mp import MultiprocessGeneratorPool

# create a VCF file


def gen_vcf_lines(random_numbers):
    # header
    yield VCFLine.as_comment_key_string("fileformat", "VCFv4.3")
    yield VCFLine.as_comment_key_dict(
        "FORMAT",
        {
            "ID": "GT",
            "Number": "1",
            "Type": "String",
            "Description": "Genotype",
        },
    )
    # ##contig=<ID=1,length=1000000,assembly=test>
    yield VCFLine.as_comment_key_dict(
        "contig",
        {
            "ID": 1,
            "length": 1000000,
            "assembly": "test",
        },  # TODO fix me to be either b37 or b38
    )
    yield VCFLine.as_comment_raw(
        "\t".join(
            (
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "TESTSAMPLE",  # customizable name of sample
            )
        )
    )
    # variants
    for pos in random_numbers:
        ref = "A"
        alt = ["T"]
        _id = f"rs{pos}"
        yield VCFLine(
            comment_raw="",
            comment_key="",
            comment_value_str="",
            comment_value_dict={},
            chrom="1",
            pos=pos,
            _id=[_id],
            ref=ref,
            alt=alt,
            qual_str=".",
            _filter=["PASS"],
            info={},
            sample=[{"GT": "0/0"}],
        )


@fixture
def random_numbers():
    rng = random.Random(42)
    return tuple(sorted(rng.sample(range(1, 1000000), 1000)))


@fixture
def vcf_temp_dir():
    with TemporaryDirectory() as tmpdir:
        yield tmpdir


@fixture
def vcf_generated_file(random_numbers, vcf_temp_dir):
    tmpdir = vcf_temp_dir
    # stored in git as txt for line editability
    # recompress into blockgzip
    # from https://stackoverflow.com/a/61979188 - override close temporarily
    # this is so the blockgzip writer writes the end of the file correctly
    # AND we can get all the bytes
    with open(os.path.join(tmpdir, "vcf.gz"), "wb") as target:
        with BlockGZipWriter(target) as bgzipfile:
            for vcf_line in gen_vcf_lines(random_numbers):
                vcf_line_bytes = str(vcf_line).encode("utf-8")
                print(vcf_line_bytes)
                bgzipfile.write(vcf_line_bytes)
                bgzipfile.write(b"\n")
    # reopen in read mode
    with open(os.path.join(tmpdir, "vcf.gz"), "rb") as target:
        yield target


@fixture
def rsidindex(vcf_generated_file):
    # get file size by going to the end
    vcf_generated_file.seek(0, SEEK_END)
    size = vcf_generated_file.tell()
    vcf_generated_file.seek(0)
    return RsidIndex.from_vcf(vcf_generated_file, size)


@fixture
def pool():
    with MultiprocessGeneratorPool() as pool:
        # with SingleProcessGeneratorPool() as pool:
        yield pool


def openfile(filename):
    return open(filename, "rb")


def passthrough(i):
    return i


class TestBasic:
    def test_lookup(
        self, vcf_temp_dir, vcf_generated_file, rsidindex, random_numbers, pool
    ):
        lines = tuple(
            rsidindex.gen_rsid_lines(
                passthrough,
                random_numbers[:10],
                openfile,
                [os.path.join(vcf_temp_dir, "vcf.gz")],
                pool,
            )
        )
        assert len(lines) == 10
