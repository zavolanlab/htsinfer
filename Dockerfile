###### BASE IMAGE ######
FROM python:3.9-slim

####### METADATA #######
LABEL base_image="python:3.9-slim"
LABEL version="1.0"
LABEL software="HTSinfer"
LABEL software.version="v0.9.0"
LABEL about.summary="HTSinfer infers metadata from Illumina high-throughput sequencing (HTS) data"
LABEL about.home="https://github.com/zavolanlab/htsinfer"
LABEL about.documentation="https://htsinfer.readthedocs.io/"
LABEL about.license_file="https://spdx.org/licenses/Apache-2.0"
LABEL about.license="Apache License 2.0"
LABEL about.tags="Transcriptomics"

###### MAINTAINER ######
LABEL maintainer="Alexander Kanitz <alexander.kanitz@alumni.ethz.ch>"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Spitalstrasse 41, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"

##### INSTALLATION #####

WORKDIR /htsinfer

COPY ./ ./

RUN apt-get update && apt-get install -y kallisto rna-star && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip && pip install biopython pandas pyahocorasick pydantic pysam cutadapt && \
    pip install -e . && pip cache purge

ENTRYPOINT ["htsinfer"]
CMD ["-h"]
