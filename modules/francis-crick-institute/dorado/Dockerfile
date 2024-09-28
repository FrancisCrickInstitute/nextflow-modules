FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

RUN apt-get update
RUN apt-get install -y --no-install-recommends curl
RUN apt-get install -y --no-install-recommends git
RUN apt-get install -y --no-install-recommends ca-certificates
RUN apt-get install -y --no-install-recommends nvidia-cuda-toolkit
RUN apt-get install -y --no-install-recommends libhdf5-dev
RUN apt-get install -y --no-install-recommends libssl-dev
RUN apt-get install -y --no-install-recommends libzstd-dev
RUN apt-get install -y --no-install-recommends python3.9
RUN apt-get install -y --no-install-recommends procps
RUN apt-get clean && rm -rf /var/lib/apt/lists/*


WORKDIR /home

ARG DORADO_VERSION=0.8.0
RUN curl -O https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-linux-x64.tar.gz
RUN tar -xzvf dorado-${DORADO_VERSION}-linux-x64.tar.gz

ENV PATH=/home/dorado-${DORADO_VERSION}-linux-x64/bin:$PATH
#RUN dorado download --model all
