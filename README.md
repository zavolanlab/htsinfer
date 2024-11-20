# HTSinfer

[![license][badge-license]][badge-url-license]
[![docs][badge-docs]][badge-url-docs]
[![release_gh][badge-release-gh]][badge-url-release-gh]
[![ci][badge-ci]][badge-url-ci]
[![coverage][badge-coverage]][badge-url-coverage]
[![release_biocontainer][badge-release-biocontainer]][badge-url-release-biocontainer]
[![DOI:zenodo][badge-doi-zenodo]][badge-url-doi-zenodo]

HTSinfer infers RNA-Seq metadata from Illumina high-throughput sequencing (HTS) data.

## Quick start

For a more in-depth guide please refer to the [HTSinfer documentation][docs-documentation].

### Installation

HTSinfer is available on [Bioconda][bioconda-release]. To install it in your currently active [Conda][conda] environment, run:

```sh
conda install bioconda::htsinfer
```

### General usage

```sh
htsinfer [-h] [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}] [--version]       
         PATH [PATH]
```

### Examples

**Single-ended library**

```sh
htsinfer tests/files/adapter_single.fastq
```

**Paired-ended library**

```sh
htsinfer tests/files/adapter_1.fastq tests/files/adapter_2.fastq
```

Output is written to `STDOUT` in JSON format. The log is written to `STDERR`.

### Example output

This is the output (`STDOUT`) of the above-mentioned call on a paired-ended
example library:

```json
{
   "library_stats": {
      "file_1": {
         "read_length": {
            "min": 75,
            "max": 75,
            "mean": 75.0,
            "median": 75,
            "mode": 75
         }
      },
      "file_2": {
         "read_length": {
            "min": 75,
            "max": 75,
            "mean": 75.0,
            "median": 75,
            "mode": 75
         }
      }
   },
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
   "library_type": {
      "file_1": "first_mate",
      "file_2": "second_mate",
      "relationship": "split_mates"
   },
   "read_orientation": {
      "file_1": "SF",
      "file_2": "SR",
      "relationship": "ISF"
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
   }
}
```

To better understand the output, please refer to the [`Results`][docs-results]
section in the [documentation][badge-url-docs].

## Versioning

The project follows the [Semantic Versioning][semver] guidelines for version management. 
Currently, the service is in its beta phase, meaning API breaking changes or updates may occur without prior notice.

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

[badge-ci]: <https://github.com/zavolanlab/htsinfer/workflows/ci/badge.svg?branch=dev>
[badge-coverage]: <https://codecov.io/gh/zavolanlab/htsinfer/branch/dev/graph/badge.svg?token=KYGJ9MUPHT>
[badge-docs]: <https://readthedocs.org/projects/htsinfer/badge/?version=latest>
[badge-license]: <https://img.shields.io/badge/license-Apache%202.0-blue.svg>
[badge-release-biocontainer]: <https://img.shields.io/badge/BioContainer-%20htsinfer-blue?style=flat.svg>
[badge-release-gh]: <https://img.shields.io/github/v/tag/zavolanlab/htsinfer?color=C39BD3>
[badge-doi-zenodo]: <https://zenodo.org/badge/265279928.svg>
[badge-url-ci]: <https://github.com/zavolanlab/htsinfer/actions?query=workflow%3Aci>
[badge-url-coverage]: <https://codecov.io/gh/zavolanlab/htsinfer>
[badge-url-docs]: <https://htsinfer.readthedocs.io/en/latest/?badge=latest>
[badge-url-license]: <http://www.apache.org/licenses/LICENSE-2.0>
[badge-url-release-biocontainer]: <https://quay.io/repository/biocontainers/htsinfer>
[badge-url-release-gh]: <https://github.com/zavolanlab/htsinfer/releases>
[badge-url-doi-zenodo]: <https://doi.org/10.5281/zenodo.13985958>
[conda]: <https://docs.conda.io/en/latest/miniconda.html>
[bioconda-release]: <https://anaconda.org/bioconda/htsinfer>
[semver]: <https://semver.org/>
[contact]: <https://zavolan.biozentrum.unibas.ch/>
[docs-documentation]: <https://htsinfer.readthedocs.io/>
[docs-results]: <https://htsinfer.readthedocs.io/en/latest/guides/examples.html#results>
[issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
[mamba]: <https://mamba.readthedocs.io/en/latest/installation.html>
