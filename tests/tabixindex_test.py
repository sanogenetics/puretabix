import os.path

import pytest

import puretabix


@pytest.fixture
def vcf():
    pth = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data",
        "CEU.exon.2010_03.genotypes.trimmed.vcf.gz",
    )
    with open(pth, "rb") as vcf:
        yield vcf


@pytest.fixture
def vcf_tbi():
    pth = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "data",
        "CEU.exon.2010_03.genotypes.trimmed.vcf.gz.tbi",
    )
    with open(pth, "rb") as vcf:
        yield vcf


class TestBasic:
    def test_basic(self, vcf, vcf_tbi):
        indexed = puretabix.TabixIndexedFile(vcf, vcf_tbi)

        fetched = indexed.fetch("1", 1108138)
        assert len(fetched.strip().split("\n")) == 1
        assert "rs61733845" in fetched

        fetched = indexed.fetch("1", 1108138 - 10, 1108138 + 10)
        assert len(fetched.strip().split("\n")) == 1
        assert "rs61733845" in fetched
