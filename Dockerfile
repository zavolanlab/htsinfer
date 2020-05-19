## Set base image
FROM python:3.7.4

## Set metadata
LABEL base.image="python:3.7.4"
LABEL software="HTSinfer"
LABEL software.version="0.1.0"
LABEL software.description="HTSinfer infers metadata from High Throughput Sequencing (HTS) data"
LABEL software.website="https://github.com/zavolanlab/htsinfer"
LABEL software.documentation="https://github.com/zavolanlab/htsinfer/blob/master/README.md"
LABEL software.license="https://github.com/zavolanlab/htsinfer/blob/master/LICENSE"
LABEL software.tags="bioinformatics, ngs, high-throughput sequencing, inference"
LABEL maintainer="Rohan Kandhari"
LABEL maintainer.email="rohan.kandhari.bme16@iitbhu.ac.in"
LABEL maintainer.organisation="Zavolan lab, Biozentrum, University of Basel"

## Set variables
ARG USER="user"
ARG GROUP="usergroup"
ARG WORKDIR="/home/${USER}/htsinfer/"
ENV PATH="${WORKDIR}:${PATH}"

## Create and set working directory
RUN mkdir -p ${WORKDIR}
WORKDIR ${WORKDIR}

## Copy and install required packages
COPY ./requirements.txt ${WORKDIR}
RUN pip install -r requirements.txt

## Copy remaining files
COPY src/ ${WORKDIR}/src
COPY tests/ ${WORKDIR}/tests
COPY LICENSE README.md ${WORKDIR}

## Set up environment
RUN groupadd -r ${GROUP} && \
    useradd --no-log-init -r -g ${GROUP} ${USER} && \
    chown -R ${USER}:${GROUP} ${WORKDIR} && \
    chmod 700 ${WORKDIR}
USER ${USER}
ENTRYPOINT ["src/htsinfer.py"]
