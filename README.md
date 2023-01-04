Pure Tabix
==========

[![Build Status](https://circleci.com/gh/sanogenetics/puretabix.svg?style=svg)](https://app.circleci.com/pipelines/github/sanogenetics/puretabix)
[![PyPI version](https://badge.fury.io/py/puretabix.svg)](https://badge.fury.io/py/puretabix)

This is a pure-python Tabix index parser. Useful as an alternative to [PySAM](https://pypi.org/project/pysam) and [PyTabix](https://pypi.org/project/pytabix)
for rapid read access by position to Tabix indexed block gzipped files such as VCFs and other common bioinfomatics formats.

See https://samtools.github.io/hts-specs/tabix.pdf and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176 for information
about Tabix and the detailed file format specification.

```py
from puretabix import TabixIndexedFile

tabix_indexed_file = TabixIndexedFile.from_files(open('somefile.vcf.gz', 'rb'), open('somefile.vcf.gz.tbi', 'rb'))
tabix_indexed_file.fetch("1", 1000, 5000)
```

Documentation is supported via Python built-in module [PyDoc](https://docs.python.org/3/library/pydoc.html): `python3 -m pydoc -b puretabix`

VCF
---

Included in this package is tooling for reading and writing VCF lines.

To read a file:

```python
from puretabix.vcf import read_vcf_lines

with open("source.vcf") as input:
    for vcfline in read_vcf_lines(input):
        if vcfline.is_comment:
            # its a comment or meta-information
            pass
        else:
            # access the parsed information
            if "PASS" not in vcfline._filter:
                print(f"{vcfline.chrom} {vcfline.pos} {vcfline.get_genotype()}")
```

To write some lines:

```python
from puretabix.vcf import VCFLine

with open("output.vcf") as output:
    output.write(str(VCFLine.as_comment_key_dict("fileformat", "VCFv4.2")))
    output.write("\n")
    output.write(
        str(
            VCFLine.as_comment_raw(
                "\t".join(
                    (
                        "CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                        "SAMPLE",
                    )
                )
            )
        )
    )
    output.write("\n")
    output.write(
        str(
            VCFLine.as_data(
                "chr1",
                123,
                ("rs123",),
                "A",
                ("C",),
                ".",
                ("PASS",),
                {},
                ({"GT": "1/0"},),
            )
        )
    )
    output.write("\n")
```
VCF with index
--------------

If there is a tabix index for a block gzipped VCF file, that index can be used for fast random access

```python
import puretabix

with open("input.vcf.gz", "rb") as vcf:
    with open("input.vcf.gz.tbi", "rb") as vcf_tbi:
        indexed = puretabix.TabixIndexedVCFFile.from_files(vcf, vcf_tbi)
        vcfline = tuple(indexed.fetch_vcf_lines("chr1", 1108138))
        assert vcfline.chrom == "chr1"
        assert vcfline.pos == 1108138
        print(f"gt = {vcfline.get_genotype()}")
```

development
-----------

TL;DR: `pip install -e '.[dev]' && pre-commit install`

```sh
pip install -e '.[dev]'  # Install using pip including development extras
pre-commit install  # Enable pre-commit hooks
pre-commit run --all-files  # Run pre-commit hooks without committing
# Note pre-commit is configured to use:
# - seed-isort-config to better categorise third party imports
# - isort to sort imports
# - black to format code
pip-compile  # Freeze dependencies
pytest  # Run tests
coverage run --source=puretabix -m pytest && coverage report -m  # Run tests, print coverage
mypy .  # Type checking
pipdeptree  # Print dependencies
scalene --outfile tests/perf_test.txt --profile-all --cpu-sampling-rate 0.0001 tests/perf_test.py  # performance measurements
```

Global git ignores per https://help.github.com/en/github/using-git/ignoring-files#configuring-ignored-files-for-all-repositories-on-your-computer

For release to PyPI see https://packaging.python.org/tutorials/packaging-projects/

```sh
git checkout master
git pull
git add setup.py CHANGES.txt
git commit -m"prepare for x.x.x"
git push
git tag x.x.x
git push origin x.x.x
python3 setup.py sdist && python3 -m twine upload dist/*
```

acknowledgements
----------------

Inspired by @yangmqglobe code in https://github.com/cggh/scikit-allel/pull/297
