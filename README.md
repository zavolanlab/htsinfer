# HTSinfer

[![license][badge-license]][badge-url-license]
[![CI][badge-ci]][badge-url-ci]
[![coverage][badge-coverage]][badge-url-coverage]

HTSinfer infers metadata from High Throughput Sequencing (HTS) data.

## Usage

```sh
htsinfer [-h] -f1 FILE [-f2 FILE] [-n INT] [--verbose] [--debug] [--version]
```

## Parameters

```console
required arguments:
  -f1 FILE, --file-1 FILE
      file path to read/first mate library

optional arguments:
  -f2 FILE, --file-2 FILE
      file path to second mate library
  -n INT, --max-records INT
      maximum number of records to process, starting with
      first record; set to 0 to process entire file(s)
  --verbose, -v
      print logging messages to STDERR
  --debug
      print debugging messages to STDERR
  --version
      show version information and exit
  -h, --help
      show this help message and exit
```

## Installation

In order to use the HTSinfer, clone the repository and install the
dependencies with [Conda][conda]:

```sh
git clone https://github.com/zavolanlab/htsinfer
cd htsinfer
conda env create --file environment.yml
conda env update --file environment-dev.yml  # optional: install development/testing dependencies
conda activate htsinfer
python setup.py install
```

If you have installed the development/testing dependencies, you can verify
that HTSinfer was installed correctly by executing the tests shipped with
the package:

```sh
python -m pytest
```

Run HTSinfer via the CLI script `htsinfer` as described in the [Usage](#Usage)
section.

## Contributing

This project lives off your contributions, be it in the form of bug reports,
feature requests, discussions, or fixes and other code changes. Please refer
to the [contributing guidelines](CONTRIBUTING.md) if you are interested to
contribute. Please mind the [code of conduct](CODE_OF_CONDUCT.md) for all
interactions with the community.

## Contact

For questions or suggestions regarding the code, please use the
[issue tracker][issue-tracker]. For any other inquiries, please contact us
by email: <zavolab-biozentrum@unibas.ch>

(c) 2020 [Zavolan lab, Biozentrum, University of Basel][contact]

[badge-ci]: <https://travis-ci.com/zavolanlab/htsinfer.svg?branch=master>
[badge-coverage]: <https://img.shields.io/coveralls/github/zavolanlab/htsinfer/master>
[badge-license]: <https://img.shields.io/badge/license-Apache%202.0-orange.svg?style=flat&color=important>
[badge-url-ci]: <https://travis-ci.com/zavolanlab/htsinfer>
[badge-url-coverage]: <https://coveralls.io/github/zavolanlab/htsinfer>
[badge-url-license]: <http://www.apache.org/licenses/LICENSE-2.0>
[conda]: <https://docs.conda.io/en/latest/miniconda.html>
[contact]: <https://zavolan.biozentrum.unibas.ch/>
[issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
