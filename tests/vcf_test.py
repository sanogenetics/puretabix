from puretabix.vcf import VCFLine, read_vcf_lines


class TestVCFLineConstructors:
    def test_comment_raw(self):
        line = VCFLine.as_comment_raw("Hello world")
        linestr = str(line)
        assert linestr == "#Hello world"
        linerepr = repr(line)
        assert linerepr == "VCFLine('Hello world','','',{},'',0,(),'',(),'',(),{},())"

    def test_comment_key_string(self):
        line = VCFLine.as_comment_key_string("Hello", "world")
        linestr = str(line)
        assert linestr == "##Hello=world"
        linerepr = repr(line)
        assert linerepr == "VCFLine('','Hello','world',{},'',0,(),'',(),'',(),{},())"

    def test_comment_key_dict(self):
        line = VCFLine.as_comment_key_dict("foo", {"Hello": "world", "spam": "eggs"})
        linestr = str(line)
        assert linestr == "##foo=<Hello=world,spam=eggs>"
        linerepr = repr(line)
        assert (
            linerepr
            == "VCFLine('','foo','',{'Hello': 'world', 'spam': 'eggs'},'',0,(),'',(),'',(),{},())"
        )

    def test_data(self):
        line = VCFLine.as_data(
            "chr1",
            123,
            ["rs1"],
            "A",
            ["C"],
            ".",
            ["PASS"],
            {"foo": "bar"},
            [{"GT": "1/1"}, {"GT": "1/1"}],
        )
        linestr = str(line)
        assert linestr == "chr1\t123\trs1\tA\tC\t.\tPASS\tfoo=b,a,r\tGT\t1/1\t1/1"
        linerepr = repr(line)
        assert (
            linerepr
            == "VCFLine('','','',{},'chr1',123,('rs1',),'A',('C',),'.',('PASS',),{'foo': 'bar'},({'GT': '1/1'}, {'GT': '1/1'}))"
        )


class TestVCFFSM:
    def test_dbsnp(self, vcf_gz):
        lines = tuple(map(bytes.decode, vcf_gz.readlines()))
        lines_parsed = read_vcf_lines(lines)

        for line_in, line_out in zip(lines, lines_parsed):
            line_in = line_in.strip()
            assert line_in == str(line_out), (line_in, line_out)

    def test_dbsnp_gt(self, vcf_gz):
        lines = tuple(map(bytes.decode, vcf_gz.readlines()))
        lines_parsed = read_vcf_lines(lines)

        for line_in, line_out in zip(lines, lines_parsed):
            line_in = line_in.strip()
            if not line_out.is_comment:
                assert line_out.get_genotype()

    def test_dbsnp_header(self, vcf_gz):
        lines = tuple(map(bytes.decode, vcf_gz.readlines()))
        lines_parsed = read_vcf_lines(lines, header_only=True)

        for line_in, line_out in zip(lines, lines_parsed):
            line_in = line_in.strip()
            assert line_out.is_comment
            assert str(line_out).startswith("#")
