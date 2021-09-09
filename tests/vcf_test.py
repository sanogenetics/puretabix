import puretabix
from puretabix.vcf import VCFLine


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
        vcf_fsm = puretabix.vcf.get_vcf_fsm()

        accumulator = puretabix.vcf.VCFAccumulator()
        lines_parsed = []
        for line in lines:
            vcf_fsm.run(line, puretabix.vcf.LINE_START, accumulator)
            vcfline = accumulator.to_vcfline()
            accumulator.reset()
            line_parsed = str(vcfline)
            lines_parsed.append(line_parsed)

        for line_in, line_out in zip(lines, lines_parsed):
            print(line_in, line_out)
            line_in = line_in.strip()
            line_out = str(line_out)
            assert line_in == line_out, (line_in, line_out)
