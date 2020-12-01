from puretabix import get_bgzip_lines_parallel


class TestBlockGZip:
    def test_get_lines(self, vcf_filename, vcf_gz):
        lines = tuple(sorted(map(bytes.decode, vcf_gz.readlines())))
        lines_parsed = tuple(sorted(get_bgzip_lines_parallel(vcf_filename)))

        for line_in, line_out in zip(lines, lines_parsed):
            print(line_in, line_out)
            line_in = line_in.strip()
            line_out = str(line_out)
            assert line_in == line_out, (line_in, line_out)
