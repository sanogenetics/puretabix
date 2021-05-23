import os.path
import tempfile

from puretabix.bgzip import (
    BlockGZipWriter,
    get_filename_parallel_lines,
    get_filename_ranged_lines,
)


class TestBlockGZip:
    def test_get_filename_ranged_lines(self, vcf_filename):
        for line in get_filename_ranged_lines(vcf_filename, 0, 100):
            print(line)

    def test_get_lines(self, vcf_filename, vcf_gz):
        lines = tuple(sorted(map(bytes.decode, vcf_gz.readlines())))
        lines_parsed = tuple(sorted(get_filename_parallel_lines(vcf_filename)))

        for line_in, line_out in zip(lines, lines_parsed):
            print(line_in, line_out)
            line_in = line_in.strip()
            line_out = str(line_out)
            assert line_in == line_out, (line_in, line_out)

    def test_write_bgzip(self, vcf_gz):
        lines = tuple(sorted(map(bytes.decode, vcf_gz.readlines())))
        with tempfile.TemporaryDirectory() as tmpdir:
            bgzfilename = os.path.join(tmpdir, "out.bzg")
            with BlockGZipWriter(open(bgzfilename, "wb")) as bgzipwriter:
                for line in lines:
                    bgzipwriter.write(line.encode())

            lines_parsed = tuple(sorted(get_filename_parallel_lines(bgzfilename)))

            for line_in, line_out in zip(lines, lines_parsed):
                print(line_in, line_out)
                line_in = line_in.strip()
                line_out = str(line_out)
                assert line_in == line_out, (line_in, line_out)
