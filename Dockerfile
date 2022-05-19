###### BASE IMAGE ######
FROM continuumio/miniconda3:4.11.0

####### METADATA #######
LABEL base_image="continuumio/miniconda3:4.11.0"
LABEL version="1.0"
LABEL software="HTSinfer"
LABEL software.version="v0.9.0"
LABEL about.summary="HTSinfer infers metadata from High Throughput Sequencing (HTS) data"
LABEL about.home="https://github.com/zavolanlab/htsinfer"
LABEL about.documentation="https://htsinfer.readthedocs.io/"
LABEL about.license_file="https://spdx.org/licenses/Apache-2.0"
LABEL about.license="Apache License 2.0"
LABEL about.tags="bioinformatics, ngs"

###### MAINTAINER ######
LABEL maintainer="Alexander Kanitz <alexander.kanitz@alumni.ethz.ch>"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Spitalstrasse 41, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"

##### INSTALLATION #####

# COPY THE YAML & INSTALL SOFTWARE WITH CONDA
WORKDIR /usr/src/app
COPY ./ ./
RUN conda env create --file environment.yml \
    && conda clean --all

# VARIABLES
ARG WORKDIR="/home/USER"
ARG USER="USER"
ARG GROUP="GROUP"
ENV PATH="${WORKDIR}:${PATH}"

# CREATE USER
RUN groupadd -r ${GROUP} && useradd --no-log-init -r -g ${GROUP} ${USER}

# SET ENVIRONMENT
WORKDIR ${WORKDIR}
RUN chown -R ${USER}:${GROUP} ${WORKDIR} && chmod 700 ${WORKDIR}
USER ${USER}
RUN echo "source activate htsinfer" > ~/.bashrc
ENV PATH /opt/conda/envs/htsinfer/bin:$PATH
