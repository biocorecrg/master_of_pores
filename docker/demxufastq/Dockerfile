#This will do biocorecrg/demxufastq:0.1
FROM biocorecrg/centos-perlbrew-pyenv23-java:centos7

# File Author / Maintainer
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com>

#upgrade pip
RUN pip install --upgrade pip
RUN pip install pandas
RUN pip install biopython
RUN pip install datetime

COPY *.py /usr/local/bin/
