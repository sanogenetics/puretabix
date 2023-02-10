import pytest

import puretabix


class TestQuery:
    @pytest.fixture
    def indexed(self, vcf, vcf_tbi):
        return puretabix.TabixIndexedFile.from_files(vcf, vcf_tbi)

    @pytest.fixture
    def indexed_vcf(self, vcf, vcf_tbi):
        return puretabix.TabixIndexedVCFFile.from_files(vcf, vcf_tbi)

    def test_hit(self, indexed):
        fetched = indexed.fetch("1", 1108138)
        assert len(fetched.strip().split("\n")) == 1, fetched
        assert "rs61733845" in fetched, fetched

        fetched = indexed.fetch("1", 1108138 - 10, 1108138 + 10)
        assert len(fetched.strip().split("\n")) == 1, fetched
        assert "rs61733845" in fetched, fetched

    def test_vcf_line(self, indexed_vcf):
        fetched = tuple(indexed_vcf.fetch_vcf_lines("1", 1108138))
        assert len(fetched) == 1, fetched
        assert "rs61733845" in fetched[0]._id, fetched

        fetched = tuple(indexed_vcf.fetch_vcf_lines("1", 1108138 - 10, 1108138 + 10))
        assert len(fetched) == 1, fetched
        assert "rs61733845" in fetched[0]._id, fetched

    def test_beyond_end(self, indexed):
        fetched = indexed.fetch("1", 245804116 + 1)
        assert fetched == "", fetched

        fetched = indexed.fetch("1", 245804116 * 2)
        assert fetched == "", fetched

    def test_before_first(self, indexed):
        fetched = indexed.fetch("1", 100)
        assert fetched == "", fetched

        fetched = indexed.fetch("1", 1105365)
        assert fetched == "", fetched


class TestCreatedQuery(TestQuery):
    # subclass the query tests
    # but instead of loading existing index files, generate the index

    @pytest.fixture
    def indexed(self, vcf):
        idx = puretabix.TabixIndex.build_from(vcf)
        return puretabix.TabixIndexedFile(vcf, idx)

    @pytest.fixture
    def indexed_vcf(self, vcf):
        idx = puretabix.TabixIndex.build_from(vcf)
        return puretabix.TabixIndexedVCFFile(vcf, idx)
