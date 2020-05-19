# HTSinfer

[![license][badge-license]][badge-url-license]
[![CI][badge-ci]][badge-url-ci]
[![coverage][badge-coverage]][badge-url-coverage]

HTSinfer infers metadata from High Throughput Sequencing (HTS) data

## Usage

```sh
htsinfer [-hv]
```

## Parameters

```console
optional arguments:
  -h, --help     show this help message and exit
  --version      show version information and exit
  -v, --verbose  print logging messages to STDERR
  --debug        also print debugging messages to STDERR
```

## Extended usage

### Run locally

In order to use the package, clone the repository and install the dependencies:

```sh
git clone https://github.com/zavolanlab/htsinfer
cd htsinfer
pip install -r requirements.txt
python setup.py install
```

> **NOTE:** You may want to install dependencies inside a virtual environment,
> e.g., using [`virtualenv`][res-virtualenv].

To verify the correct installation, you can run the tests that are shipped with
the package:

```sh
python -m pytest
```

You can then run the package via the CLI script `htsinfer` as described in the
[Usage](#Usage) section.

### Run inside container

If you have [Docker][res-docker] installed, you can also pull the Docker
image:

```sh
docker pull zavolab/htsinfer:latest
```

The script can be found in directory `/home/user/htsinfer/src` inside the
Docker container.

The Docker image contains an entry point for the CLI script `htsinfer` so that
you can run it, e.g., with:

```sh
docker run --rm -it zavolab/htsinfer:latest --help
```

> **NOTE:** To run the tool on your own data in that manner, you will probably
> need to [mount a volume][res-docker-volume] to allow the container read input
> files and write persistent output from/to the host file system.

## Contributing

This project lives off your contributions, be it in the form of bug reports,
feature requests, discussions, or fixes and other code changes. Please refer
to the [contributing guidelines](CONTRIBUTING.md) if you are interested to
contribute. Please mind the [code of conduct](CODE_OF_CONDUCT.md) for all
interactions with the community.

## Contact

For questions or suggestions regarding the code, please use the
[issue tracker][res-issue-tracker]. For any other inquiries, please contact us
by email: <zavolab-biozentrum@unibas.ch>

(c) 2020 [Zavolan lab, Biozentrum, University of Basel][res-zavolab]

[badge-ci]: <https://travis-ci.com/zavolanlab/htsinfer.svg?branch=master>
[badge-coverage]: <https://img.shields.io/coveralls/github/zavolanlab/htsinfer/master>
[badge-license]: <https://img.shields.io/badge/license-Apache%202.0-orange.svg?style=flat&color=important>
[badge-url-ci]: <https://travis-ci.com/zavolanlab/htsinfer>
[badge-url-coverage]: <https://coveralls.io/github/zavolanlab/htsinfer>
[badge-url-license]: <http://www.apache.org/licenses/LICENSE-2.0>
[res-docker]: <https://www.docker.com/>
[res-docker-vol]: <https://docs.docker.com/storage/volumes/>
[res-issue-tracker]: <https://github.com/zavolanlab/htsinfer/issues>
[res-virtualenv]: <https://virtualenv.pypa.io/en/latest/>
[res-zavolab]: <https://zavolan.biozentrum.unibas.ch/>
