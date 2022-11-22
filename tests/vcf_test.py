from puretabix.vcf import VCFLine, read_vcf_lines


class TestVCFLineConstructors:
    def test_comment_raw(self):
        line = VCFLine.as_comment_raw("Hello world")
        linestr = str(line)
        assert linestr == "#Hello world"

    def test_comment_key_string(self):
        line = VCFLine.as_comment_key_string("Hello", "world")
        linestr = str(line)
        assert linestr == "##Hello=world"

    def test_comment_key_dict(self):
        line = VCFLine.as_comment_key_dict("foo", {"Hello": "world", "spam": "eggs"})
        linestr = str(line)
        assert linestr == "##foo=<Hello=world,spam=eggs>"


class TestVCFFSM:
    def test_dbsnp(self, vcf_gz):
        lines = tuple(map(bytes.decode, vcf_gz.readlines()))
        lines_parsed = read_vcf_lines(lines)

        for line_in, line_out in zip(lines, lines_parsed):
            line_in = line_in.strip()
            print(line_in)
            print(str(line_out).strip())
            assert line_in == str(line_out), (line_in, line_out)

    def test_dbsnp_gt(self, vcf_gz):
        lines = tuple(map(bytes.decode, vcf_gz.readlines()))
        lines_parsed = read_vcf_lines(lines)

        for line_in, line_out in zip(lines, lines_parsed):
            line_in = line_in.strip()
            print(line_in)
            print(str(line_out).strip())
            if not line_out.is_comment:
                assert line_out.get_genotype()
