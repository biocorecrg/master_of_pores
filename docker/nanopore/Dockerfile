FROM biocorecrg/centos-perlbrew-pyenv23-java:centos7

# File Author / Maintainer
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

ARG MULTIQC_VERSION=1.7
ARG SAMTOOLS_VERSION=1.4.1
#ARG ALBACORE_VERSION=2.1.7
ARG JVARKIT_VERSION=f39087b
ARG NANOFILT_VERSION=2.3.0
ARG GRAPHMAP_VERSION=eb8c75d
ARG LASTZ_VERSION=1.04.00
ARG PICARD_VERSION=2.20.0
ARG R_VERSION=3.6.0
ARG POREQUALITY_VERSION=4fb2f6e
ARG MINIONQC_VERSION=1.4.1
ARG EPINANO_VERSION=72b67b9 
ARG NANOPOLISH_VERSION=v0.11.1
ARG MINIMAP_VERSION=2.17
#ARG GUPPY_VERSION=3.1.5
ARG NANOPLOT_VERSION=1.26.3
ARG FASTQC_VERSION=0.11.8
ARG HTSEQ_VERSION=0.11.1
ARG NANOCOUNT_VERSION=0.1.a4
ARG TOMBO_VERSION=1.5 
ARG ONTFAST5API_VERSION=1.4.4
ARG TAILFINDR_VERSION=6865f84

# Environment variables
ENV POREQUAL_PATH=$PWD/POREquality/
ENV MINIONQC=/usr/local/bin/MinIONQC.R
ENV EPINANOMAIN=$PWD/EpiNano/scripts/main
ENV EPINANOMISC=$PWD/EpiNano/scripts/misc

#upgrade pip
RUN pip install --upgrade pip

#INSTALLING FASTQC
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip > fastqc.zip'
RUN unzip fastqc.zip; chmod 775 FastQC/fastqc; ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

#Installing HTSEQ
RUN pip install numpy matplotlib  
RUN pip install HTSEQ==${HTSEQ_VERSION}

# Installing nanocount
RUN pip install NanoCount==${NANOCOUNT_VERSION}

# Installing Tombo
RUN pip install ont-tombo==${TOMBO_VERSION}

# Installing samtools
RUN yum install -y xz-devel.x86_64
RUN bash -c 'curl -k -L https://downloads.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 > samtools.tar.bz2'
RUN tar -jvxf samtools.tar.bz2
RUN cd samtools-${SAMTOOLS_VERSION}; ./configure; make; make install; cd ../ 
RUN rm samtools.tar.bz2

# Installing Albacore
#RUN yum install -y epel-release-7-11.noarch
#RUN yum install -y p7zip.x86_64 p7zip-plugins.x86_64
#COPY programs/ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl ./
#RUN pip install ./ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl
#RUN rm ont_albacore-${ALBACORE_VERSION}-cp36-cp36m-manylinux1_x86_64.whl

# Installing Guppy
#RUN yum install -y epel-release-7-11.noarch 
#RUN yum install -y p7zip.x86_64 p7zip-plugins.x86_64
#RUN yum update -y
#COPY programs/ont-guppy_${GUPPY_VERSION}_linux64.7z ./
#RUN 7z x ont-guppy_${GUPPY_VERSION}_linux64.7z
#RUN rm ont-guppy_${GUPPY_VERSION}_linux64.7z
# don't install nvidia-drivers & cuda
# as it should use system drivers https://stackoverflow.com/a/44578838/632242 http://www.morethantechnical.com/2018/01/27/an-automatic-tensorflow-cuda-docker-jupyter-machine-on-google-cloud-platform/
#RUN yum install -y kernel-devel kernel-headers
#RUN curl http://developer.download.nvidia.com/compute/cuda/repos/rhel7/x86_64/7fa2af80.pub > 7fa2af80.pub; rpm --import 7fa2af80.pub
#RUN yum-config-manager --add-repo http://developer.download.nvidia.com/compute/cuda/repos/rhel7/x86_64/cuda-rhel7.repo
#RUN yum clean all
#RUN yum -y install nvidia-driver-latest-dkms cuda
#COPY programs/cuda-repo-rhel7-10.1.168-1.x86_64.rpm ./
#RUN rpm -i cuda-repo-rhel7-10.1.168-1.x86_64.rpm; yum clean all; yum -y install cuda
#RUN rm cuda-repo-rhel7-10.1.168-1.x86_64.rpm
# alternatively you can just run
#RUN yum install -y kernel-devel kernel-headers
#RUN yum-config-manager --add-repo http://developer.download.nvidia.com/compute/cuda/repos/rhel7/x86_64/cuda-rhel7.repo
#RUN yum clean all
#RUN yum -y install nvidia-driver-latest-dkms cuda

# install sam2tsv // version is github tag
RUN git clone "https://github.com/lindenb/jvarkit.git"; cd jvarkit; git checkout ${JVARKIT_VERSION}; ./gradlew sam2tsv

# Installing nanofilt
RUN pip install nanofilt==${NANOFILT_VERSION}

# Installing graphmap -- ONT reads aligner // version is github tag
RUN yum install -y libstdc++-static
RUN git clone "https://github.com/isovic/graphmap.git"; cd graphmap; git checkout ${GRAPHMAP_VERSION}; make modules; make; make install

# Installing minimap2
RUN bash -c 'curl -k -L https://github.com/lh3/minimap2/releases/download/v${MINIMAP_VERSION}/minimap2-${MINIMAP_VERSION}_x64-linux.tar.bz2 > minimap.tar.bz2'
RUN tar -jvxf minimap.tar.bz2; cd minimap2-${MINIMAP_VERSION}_x64-linux; ln -s $PWD/minimap2 /usr/local/bin/minimap2; cd ../

# Installing lastz
RUN bash -c 'curl -k -L http://www.bx.psu.edu/~rsharris/lastz/lastz-${LASTZ_VERSION}.tar.gz > lastz.tar.gz'
RUN tar -zvxf lastz.tar.gz; cd lastz-distrib-${LASTZ_VERSION}; make && make install; cd ../

# Installing Picard
RUN bash -c 'curl -k -L https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar > /usr/local/bin/picard.jar'

#Installing R
RUN yum install -y epel-release libxml2-devel libcurl-devel libcurl libcurl-devel
RUN yum install R-${R_VERSION} -y
RUN mkdir -p /usr/share/doc/R-${R_VERSION}/html
# problem with curl 4.1 // must use 4.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/curl/curl_4.0.tar.gz',repo=NULL,type='source')"

RUN Rscript -e "install.packages(c('data.table','flexdashboard','dplyr','plyr','ggExtra','ggplot2','hexbin','knitr','optparse','RColorBrewer','reshape2'), repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('devtools', 'futile.logger','readr','scales','viridis','yaml', 'ggplot2'), dependencies = TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages(c('BiocManager', 'snow'), repos='http://cran.us.r-project.org')"; Rscript -e "BiocManager::install()"


# Installing tailfindr 
RUN bash -c 'curl -k -L https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/src/hdf5-1.8.21.tar.gz > hdf5-1.8.21.tar.gz'
RUN tar -zvxf hdf5-1.8.21.tar.gz; cd hdf5-1.8.21; \
./configure --prefix=/usr/local/hdf5 --enable-fortran --enable-cxx; make; make install
ENV LD_LIBRARY_PATH "/usr/local/hdf5/lib" 
ENV PATH "/usr/local/hdf5/bin":${PATH} 

RUN Rscript -e "BiocManager::install(c('BiocGenerics', 'Biostrings', 'IRanges', 'zlibbioc', 'XVector', 'bit', 'pryr', 'hdf5r', 'foreach', 'doParallel', 'doSNOW', 'rbokeh', 'psych', 'gistr', 'mnormt'))"
RUN Rscript -e "BiocManager::install(c('Rsamtools'))"
RUN Rscript -e "BiocManager::install(c('tidyr'))"
RUN Rscript -e "devtools::install_github(\"adnaniazi/tailfindr\", dependencies=FALSE)"

# Installing POREquality // version is github tag
RUN yum install -y pandoc
RUN git clone https://github.com/carsweshau/POREquality; cd POREquality/; git checkout ${POREQUALITY_VERSION}; cd ../ 

# Installing MinIONQC
RUN bash -c 'curl -k -L https://github.com/roblanf/minion_qc/releases/download/v${MINIONQC_VERSION}/MinIONQC.R > /usr/local/bin/MinIONQC.R'
RUN chmod +x /usr/local/bin/MinIONQC.R

# Installing epinano // version is github tag
#RUN matplotlib
#numpy	1.15.4
#pandas	0.23.4
#sklearn	0.20.2

RUN git clone https://github.com/biocorecrg/EpiNano.git 
RUN rm -fr EpiNano/.git

# Installing Nanopolish
RUN yum install -y wget 
RUN yum group install "Development Tools" -y
RUN git clone --recursive https://github.com/jts/nanopolish.git
RUN cd nanopolish; git checkout ${NANOPOLISH_VERSION}; make all

# Installing ont-fast5-api
RUN pip install ont-fast5-api==${ONTFAST5API_VERSION}

#install Nanoplot
RUN pip install pandas 
RUN pip install NanoPlot==${NANOPLOT_VERSION}

# Installing MULTIQC_VERSION 
RUN pip install multiqc==${MULTIQC_VERSION}

#Installing dependencies for EpiNano python2 (temporary)
RUN pip2 install --upgrade pip; pip2 install matplotlib numpy pandas h5py

#Adding a module for Tombo
RUN pip install pyfaidx

# Clean cache
RUN yum clean all 

#cleaning
RUN rm -fr *.tar.gz; rm -fr *.bz2; rm -fr ./tmp
RUN rm -rf /var/cache/yum
RUN rm -fr fastqc.zip  samtools-* graphmap lastz-distrib-* nanopolish/*.tar.*  

ENV PATH "/project/ont-guppy_3.1.5_linux64/ont-guppy/bin/:${PATH}" 
ENV PATH "/project/nanopolish/:${PATH}" 
ENV PATH "/project/EpiNano/scripts/misc/:/project/EpiNano/scripts/main/:${PATH}"
ENV SAM2TSV "java -jar /project/jvarkit/dist/sam2tsv.jar"
ENV PICARD "java -jar /usr/local/bin/picard.jar"
ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8
