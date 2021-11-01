import os.path
import tempfile

from puretabix.bgzip import BlockGZipWriter


class TestBlockGZip:
    def test_get_lines(self, vcf_bgzreader, vcf_gz):
        lines = tuple(sorted(map(bytes.decode, vcf_gz.readlines())))
        print(f"read {len(lines)} linesfor testing")

        vcf_bgzreader.seek(0)
        try:
            lines_parsed = tuple(sorted(vcf_bgzreader.generate_lines()))
        except EOFError:
            print("caught EOF")
        print(f"read {len(lines_parsed)} lines")

        assert len(lines) == len(lines_parsed)

        for line_in, line_out in zip(lines, lines_parsed):
            print(line_in, line_out)
            line_in = line_in.strip()
            line_out = line_out.decode()
            assert line_in == line_out, (line_in, line_out)

    def test_write_bgzip(self, vcf_bgzreader, vcf_gz):
        lines = tuple(sorted(map(bytes.decode, vcf_gz.readlines())))
        with tempfile.TemporaryDirectory() as tmpdir:
            bgzfilename = os.path.join(tmpdir, "out.bzg")
            with BlockGZipWriter(open(bgzfilename, "wb")) as bgzipwriter:
                for line in lines:
                    bgzipwriter.write(line.encode())

            vcf_bgzreader.seek(0)
            lines_parsed = tuple(sorted(vcf_bgzreader.generate_lines()))

            for line_in, line_out in zip(lines, lines_parsed):
                print(line_in, line_out)
                line_in = line_in.strip()
                line_out = line_out.decode()
                assert line_in == line_out, (line_in, line_out)
