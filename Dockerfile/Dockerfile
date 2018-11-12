FROM debian:jessie

MAINTAINER Bruno Andrade <bgabriel.andrade.visconti@gmail.com>



RUN echo "deb http://http.debian.net/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git \
    wget \
    curl \
    build-essential \
    python-dev \
    python-pip \
    python-numpy \
    python-pandas \
    python-matplotlib \
    python3 \
    python3-pandas \
    python3-numpy \  
    bowtie2 \
    git \
    ca-certificates \
    bzip2 \
    unzip \
    sudo \
    locales \
    libatlas3-base \
    libfreetype6-dev \
    mercurial \
    && apt-get clean

#Install OpenJDK8
RUN apt-get install -t jessie-backports openjdk-8-jre-headless --yes

#Install biopython
RUN pip install biopython

#Install R libraries
#RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com', dependencies = TRUE);" && \
 #       R -e "install.packages('ggfortify', repos='http://cran.rstudio.com', dependencies = TRUE);" && \
  #      R -e "install.packages('stringr', repos='http://cran.rstudio.com', dependencies = TRUE);"


#Install Metaphlan2
RUN wget https://bitbucket.org/biobakery/metaphlan2/get/default.zip && \
	unzip default.zip -d /opt/
RUN python /opt/biobakery-metaphlan2-*/metaphlan2.py --input_type fastq

#Install FastQC and BBmap
RUN wget -q -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    rm fastqc.zip

RUN wget -q -O BBMap_37.09.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_37.09.tar.gz && \
	tar -xzf BBMap_37.09.tar.gz && \
	mv bbmap/ /opt/ && \
	rm BBMap_37.09.tar.gz

RUN curl -fsSL get.nextflow.io | bash && \
	mv nextflow /opt/
RUN apt install python3-tk --yes
RUN apt install  python3-pip --yes

RUN pip3 install pandas --upgrade
RUN apt install pkg-config
#RUN pip3 install matplotlib

#Clonning Scripts

RUN git clone https://github.com/brunoand/NMP.git && \
        chmod +x NMP/Scripts/Immunoadept_barchart.py && \
	echo e && \
        cp NMP/Scripts/Immunoadept_barchart.py /opt/

RUN wget http://ftp.cs.stanford.edu/mirrors/ubuntu/pool/universe/m/metaphlan2/metaphlan2_2.6.0+ds.orig.tar.xz && \
	tar -xf metaphlan2_2.6.0+ds.orig.tar.xz && \
	mv biobakery-metaphlan2-c43e40a443ed /opt

ENV PATH $PATH:/opt/::/opt/FastQC/:/opt/bbmap:/opt/biobakery-metaphlan2-7df8753d7e5e/
