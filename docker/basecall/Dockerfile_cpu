FROM nvidia/cuda:10.1-base-ubuntu18.04 

MAINTAINER Leszek Pryszcz <l.p.pryszcz+docker@gmail.com>

ARG GUPPY_VERSION=3.1.5
ARG ALBACORE_VERSION=2.1.7
ARG PORECHOP_VERSION=0.2.4

# Install Python3 & pip, curl
RUN apt update && apt install -y python3-setuptools curl
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py && rm get-pip.py

RUN mkdir -p /project
WORKDIR /project

# Install guppy & its dependencies
RUN apt update && apt install -y p7zip-full libidn11 libgssapi-krb5-2
COPY programs/ont-guppy_${GUPPY_VERSION}_linux64.7z ./
RUN 7z x ont-guppy_${GUPPY_VERSION}_linux64.7z
RUN rm ont-guppy_${GUPPY_VERSION}_linux64.7z

# Install albacore
COPY programs/ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl ./
RUN pip3 install ./ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl
RUN rm ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl

# Installing Porechop
RUN apt update && apt install -y g++ make
RUN bash -c 'curl -k -L https://github.com/rrwick/Porechop/archive/v${PORECHOP_VERSION}.tar.gz > porechop.tar.gz'
RUN tar -zvxf porechop.tar.gz; cd Porechop-${PORECHOP_VERSION}; python setup.py install; cd ../
RUN rm -fr Porechop-0.2.4 porechop.tar.gz

# Cleanup
RUN apt-get clean && apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* 
RUN locale-gen en_US.UTF-8

ENV PATH "/project/ont-guppy_3.1.5_linux64/ont-guppy/bin/:${PATH}"
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/cuda-10.1/compat/
WORKDIR /
