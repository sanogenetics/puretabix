import puretabix


class TestBasic:
    def test_basic(self, vcf, vcf_tbi):
        indexed = puretabix.TabixIndexedFile.from_files(vcf, vcf_tbi)

        fetched = indexed.fetch("1", 1108138)
        assert len(fetched.strip().split("\n")) == 1, fetched
        assert "rs61733845" in fetched, fetched

        fetched = indexed.fetch("1", 1108138 - 10, 1108138 + 10)
        assert len(fetched.strip().split("\n")) == 1, fetched
        assert "rs61733845" in fetched, fetched

    def test_beyond_end(self, vcf, vcf_tbi):
        indexed = puretabix.TabixIndexedFile.from_files(vcf, vcf_tbi)
        fetched = indexed.fetch("1", 245804116 + 1)
        assert fetched == "", fetched

        fetched = indexed.fetch("1", 245804116 * 2)
        assert fetched == "", fetched

    def test_before_first(self, vcf, vcf_tbi):
        indexed = puretabix.TabixIndexedFile.from_files(vcf, vcf_tbi)

        fetched = indexed.fetch("22", 100)
        assert fetched == "", fetched

        fetched = indexed.fetch("22", 16042344)
        assert fetched == "", fetched

    def test_index(self, vcf):
        idx = puretabix.TabixIndex.build_from(vcf)
        print(idx)
        assert idx.indexes
