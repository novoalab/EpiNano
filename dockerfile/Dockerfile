## Parent image (first layer)
FROM biocorecrg/centos-perlbrew-pyenv-java:latest
#FROM centos/python-36-centos7
ARG PYTHON_VERSION=3.6.3
USER root 
RUN bash


## Install wget
RUN pwd 
RUN yum -y upgrade
RUN yum install -y which
RUN yum -y install wget libcurl-devel

## Install python3.7
WORKDIR /usr/local/bin
RUN yum -y install gcc openssl-devel bzip2-devel libffi-devel zlib-devel xz-devel
RUN wget https://www.python.org/ftp/python/3.6.3/Python-3.6.3.tgz
RUN tar xzf Python-3.6.3.tgz

WORKDIR /usr/local/bin/Python-3.6.3
RUN ./configure --enable-optimizations
RUN yum install make -y
RUN make altinstall
RUN echo 'current location'

#RUN rm python3
WORKDIR /usr/bin
RUN ln -fs /usr/local/bin/Python-3.6.3/python python3
#RUN ln -fs /usr/local/bin/Python-3.6.3/python python
RUN chmod +x python; chmod +x python3
#RUN ln -fs /usr/local/bin/pip pip3

# Intall pip
RUN yum -y install epel-release
RUN curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py
RUN python3 get-pip.py
RUN python -m pip install --upgrade pip
RUN pip3 install requests paho-mqtt

## Install samtools
WORKDIR /usr/local/bin 
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
RUN tar xjf samtools-1.3.tar.bz2
RUN ls -lrt 
RUN pwd
RUN chmod +x samtools-1.3/configure
RUN cd samtools-1.3; ./configure ; make; make install
RUN cd samtools-1.3/htslib-1.3; ./configure ; make
RUN cd /usr/local/bin ;  ln -fs samtools-1.3/samtools .
RUN cd /usr/bin/; ln -fs /usr/local/bin/samtools-1.3/samtools samtools 

## Install sam2tsv
WORKDIR /usr/local/bin
RUN wget https://github.com/enovoa/EpiNano/raw/master/misc/sam2tsv.jar

## isntall dependencies
RUN pip3 install --upgrade pip
RUN pip3 install atomicwrites==1.4.0 attrs==21.2.0 biopython==1.76 cloudpickle==1.3.0
RUN pip3 install dask==2.5.2 fsspec==2021.6.1 future==0.17.1 h5py==2.10.0 importlib-metadata==4.6.1
RUN pip3 install locket==0.2.1  more-itertools==8.8.0 numpy==1.17.2 pandas==0.24.2
RUN pip3 install partd==1.2.0 pluggy==0.13.1 py==1.10.0 pysam==0.15.3 pytest==4.4.1 python-dateutil==2.8.1
RUN pip3 install pytz==2021.1 scikit-learn==0.20.2 scipy==1.5.4 six==1.16.0 toolz==0.11.1 typing-extensions==3.10.0.0 zipp==3.5.0




### indtall Epinnao
WORKDIR /usr/local/bin
RUN git clone --depth 2 https://github.com/enovoa/EpiNano.git

# Remove all archives
RUN rm -f *.zip *.bz2 *.gz

# create a data folder to process 
RUN mkdir -p /project/data

#env path
ENV PATH /usr/local/bin:$PATH

