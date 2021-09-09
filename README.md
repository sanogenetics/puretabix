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

tabix_indexed_file = TabixIndexedFile(open('somefile.vcf.gz', 'rb'), open('somefile.vcf.gz.tbi', 'rb'))
tabix_indexed_file.fetch("1", 1000, 5000)
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
python3 setup.py sdist && python3 -m twine upload dist/*
```

acknowledgements
----------------

Inspired by @yangmqglobe code in https://github.com/cggh/scikit-allel/pull/297
