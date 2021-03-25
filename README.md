# HTSinfer

[![license][badge-license]][badge-url-license]
[![ci][badge-ci]][badge-url-ci]
[![docs][badge-docs]][badge-url-docs]
[![coverage][badge-coverage]][badge-url-coverage]

HTSinfer infers metadata from High Throughput Sequencing (HTS) data.

## Usage

```sh
htsinfer [--output-directory PATH] [--temporary-directory PATH]
         [--cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}]
         [--records INT] [--adapter-designation-min-match-percentage FLOAT]
         [--adapter-designation-frequency-ratio FLOAT] [--adapters FILE]
         [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}] [-h]
         [--version]
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
  --cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}
                        determine which data to keep after each run; in default
                        mode, both temporary data and results are kept when
                        '--verbosity' is set to 'DEBUG', no data is kept when
                        all metadata could be successfully determined, and only
                        results are kept otherwise (default: default)
  --records INT         number of records to process; if set to ``0`` or if the
                        specified value equals or exceeds the number of
                        available records, all records will be processed
                        (default: 0)
  --adapter-designation-min-match-percentage FLOAT
                        minimum percentage of reads that contain a given adapter
                        in order for that adapter sequence to be considered as the
                        resulting sequence
  --adapter-designation-frequency-ratio FLOAT
                        the minimum frequency ratio between the first and second most
                        frequent adapter in order for an adapter sequence to be returned
                        as the resulting sequence
  --adapters FILE       adapter file containing the list of all adapter sequences
                        that neeeds to be searched in the FASTQ files
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

## API documentation

Auto-built API documentation is hosted on [ReadTheDocs][badge-url-docs].

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
[badge-docs]: <https://readthedocs.org/projects/htsinfer/badge/?version=latest>
[badge-license]: <https://img.shields.io/badge/license-Apache%202.0-blue.svg>
[badge-url-ci]: <https://travis-ci.com/zavolanlab/htsinfer>
[badge-url-coverage]: <https://coveralls.io/github/zavolanlab/htsinfer>
[badge-url-docs]: <https://htsinfer.readthedocs.io/en/latest/?badge=latest>
[badge-url-license]: <http://www.apache.org/licenses/LICENSE-2.0>
[conda]: <https://docs.conda.io/en/latest/miniconda.html>
[contact]: <https://zavolan.biozentrum.unibas.ch/>
[issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
