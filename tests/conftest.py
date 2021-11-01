import gzip
import os

import pytest

from puretabix.bgzip import BlockGZipReader


@pytest.fixture
def vcf_filename():
    pth = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data",
        "CEU.exon.2010_03.genotypes.trimmed.vcf.gz",
    )
    return pth


@pytest.fixture
def vcf(vcf_filename):
    with open(vcf_filename, "rb") as vcf:
        yield vcf


@pytest.fixture
def vcf_tbi(vcf_filename):
    with open(vcf_filename + ".tbi", "rb") as vcf:
        yield vcf


@pytest.fixture
def vcf_gz(vcf_filename):
    with gzip.open(vcf_filename, "rb") as vcf:
        yield vcf


@pytest.fixture
def vcf_bgzreader(vcf_filename):
    with open(vcf_filename, "rb") as vcf:
        yield BlockGZipReader(vcf)
