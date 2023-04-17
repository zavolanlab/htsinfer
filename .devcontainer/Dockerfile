# Include mamba from the base image
FROM condaforge/mambaforge:latest

# Install sudo & C compiler in dev container
RUN apt-get update && apt-get install sudo gcc --yes && apt-get clean

# These pkgs need to be installed in base env.
# https://github.com/conda/conda-build/issues/3813
RUN conda install -c conda-forge \
  "conda-build >= 3.23.3" "boa >= 0.14.0" "conda-verify >= 3.1.1"
