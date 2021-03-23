# HTSinfer

[![license][badge-license]][badge-url-license]
[![CI][badge-ci]][badge-url-ci]
[![coverage][badge-coverage]][badge-url-coverage]

HTSinfer infers metadata from High Throughput Sequencing (HTS) data.

## Usage

```sh
htsinfer [--output-directory PATH] [--temporary-directory PATH]
         [--cleanup-regime {default,keep_all,keep_none,keep_results}]
         [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}] [-h] [--version]
         FASTQ_PATH [FASTQ_PATH]
```

## Parameters

```console
positional arguments:
  FASTQ_PATH            either one or two file paths to the read library to be
                        evaluated.

optional arguments:
  --output-directory PATH
                        path to directory where output is written to (default:
                        current working directory)
  --temporary-directory PATH
                        path to directory where temporary output is written to
                        (default: default temporary directory)
  --cleanup-regime {default,keep_all,keep_none,keep_results}
                        determine which data to keep after each run; in default
                        mode, both temporary data and results are kept when
                        '--verbosity' is set to 'DEBUG', no data is kept when
                        all metadata could be successfully determined, and only
                        results are kept otherwise (default: default)
  --verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}
                        logging verbosity level (default: INFO)
  -h, --help            show this help message and exit
  --version             show version information and exit
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
