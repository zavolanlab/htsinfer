# HTSinfer

[![license][badge-license]][badge-url-license]
[![docs][badge-docs]][badge-url-docs]
[![release_gh][badge-release-gh]][badge-url-release-gh]
[![release_docker][badge-release-docker]][badge-url-release-docker]
[![ci][badge-ci]][badge-url-ci]
[![coverage][badge-coverage]][badge-url-coverage]

HTSinfer infers metadata from Illumina high-throughput sequencing (HTS) data.

## Examples

**Single-ended library***

```sh
htsinfer tests/files/adapter_single.fastq
```

**Paired-ended library***

```sh
htsinfer tests/files/adapter_1.fastq tests/files/adapter_2.fastq
```

Output is written to `STDOUT` in JSON format. The log is written to `STDERR`.

### Example output

This is the output (`STDOUT`) of the above-mentioned call on a paired-ended
example library:

```json
{
   "library_source": {
      "file_1": {
         "short_name": "hsapiens",
         "taxon_id": "9606"
      },
      "file_2": {
         "short_name": "hsapiens",
         "taxon_id": "9606"
      }
   },
   "library_stats": {
      "file_1": {
         "read_length": {
            "max": 75,
            "min": 75
         }
      },
      "file_2": {
         "read_length": {
            "max": 75,
            "min": 75
         }
      }
   },
   "library_type": {
      "file_1": "first_mate",
      "file_2": "second_mate",
      "relationship": "split_mates"
   },
   "read_layout": {
      "file_1": {
         "adapt_3": "AATGATACGGCGACC",
         "polyA_frac": 10.0
      },
      "file_2": {
         "adapt_3": "AATGATACGGCGACC",
         "polyA_frac": 10.0
      }
   },
   "read_orientation": {
      "file_1": "SF",
      "file_2": "SR",
      "relationship": "ISF"
   }
}
```

To better understand the output, please refer to the [`Results`
model][docs-api-results] in the [API documentation][badge-url-docs]. Note that
`Results` model has several nested child models, such as enumerators of
possible outcomes. Simply follow the references in each parent model for
detailed descriptions of each child model's attributes.

## General usage

```sh
htsinfer [--output-directory PATH]
         [--temporary-directory PATH]
         [--cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}]
         [--records INT]
         [--threads INT]
         [--transcripts FASTA]
         [--read-layout-adapters PATH]
         [--read-layout-min-match-percentage FLOAT]
         [--read-layout-min-frequency-ratio FLOAT]
         [--library-source-min-match-percentage FLOAT]
         [--library-source-min-frequency-ratio FLOAT]
         [--library-type-max-distance INT]
         [--library-type-mates-cutoff FLOAT]
         [--read-orientation-min-mapped-reads INT]
         [--read-orientation-min-fraction FLOAT]
         [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}]
         [-h] [--version]
         PATH [PATH]
```

## Installation

In order to use the HTSinfer, clone the repository and install the
dependencies via [Conda][conda]:

```sh
git clone https://github.com/zavolanlab/htsinfer
cd htsinfer
conda env create --file environment.yml
conda env update --file environment-dev.yml  # optional: install development/testing dependencies
```

Note that creating the environment takes non-trivial time and it is strongly
recommended that you install [Mamba][mamba] and replace `conda` with `mamba`
in the previous commands.

Then, activate the `htsinfer` Conda environment with:

```sh
conda activate htsinfer
```

If you have installed the development/testing dependencies, you may first want
to verify that HTSinfer was installed correctly by executing the tests shipped
with the package:

```sh
python -m pytest
```

Otherwise just go ahead and try one of the [examples](#Examples).

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
[badge-coverage]: <https://codecov.io/gh/zavolanlab/htsinfer/branch/dev/graph/badge.svg?token=KYGJ9MUPHT>
[badge-docs]: <https://readthedocs.org/projects/htsinfer/badge/?version=latest>
[badge-license]: <https://img.shields.io/badge/license-Apache%202.0-blue.svg>
[badge-release-docker]: <https://img.shields.io/docker/image-size/zavolab/htsinfer?color=C39BD3&label=docker>
[badge-release-gh]: <https://img.shields.io/github/v/tag/zavolanlab/htsinfer?color=C39BD3>
[badge-url-ci]: <https://travis-ci.com/zavolanlab/htsinfer>
[badge-url-coverage]: <https://codecov.io/gh/zavolanlab/htsinfer>
[badge-url-docs]: <https://htsinfer.readthedocs.io/en/latest/?badge=latest>
[badge-url-license]: <http://www.apache.org/licenses/LICENSE-2.0>
[badge-url-release-docker]: <https://hub.docker.com/repository/docker/zavolab/htsinfer>
[badge-url-release-gh]: <https://github.com/zavolanlab/htsinfer/releases>
[conda]: <https://docs.conda.io/en/latest/miniconda.html>
[contact]: <https://zavolan.biozentrum.unibas.ch/>
[docs-api-results]: <https://htsinfer.readthedocs.io/en/latest/modules/htsinfer.html#htsinfer.models.Results>
[issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
[mamba]: <https://mamba.readthedocs.io/en/latest/installation.html>
