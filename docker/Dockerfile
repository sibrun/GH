FROM ubuntu

ARG DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /root

# Base system
RUN apt update && apt -y install git wget

# sage
RUN apt -y install sagemath

# Download packages
RUN wget https://github.com/xianyi/OpenBLAS/releases/download/v0.3.20/OpenBLAS-0.3.20.tar.gz && \
    tar -xvzf OpenBLAS-0.3.20.tar.gz && \
    wget https://github.com/linbox-team/givaro/releases/download/v4.2.0/givaro-4.2.0.tar.gz && \
    tar -xvzf givaro-4.2.0.tar.gz  && \
    wget https://github.com/linbox-team/fflas-ffpack/releases/download/v2.5.0/fflas-ffpack-2.5.0.tar.gz && \
    tar -xzvf fflas-ffpack-2.5.0.tar.gz &&\
    wget https://github.com/linbox-team/linbox/releases/download/v1.7.0/linbox-1.7.0.tar.gz &&\
    tar -xzvf linbox-1.7.0.tar.gz && \ 
    wget https://pallini.di.uniroma1.it/nauty27r3.tar.gz && \
    tar -xzvf nauty27r3.tar.gz 

ENV LD_LIBRARY_PATH=/opt/OpenBLAS/lib


WORKDIR /root/OpenBLAS-0.3.20

RUN make && make install

WORKDIR /root/givaro-4.2.0

RUN ./configure && make && make install


WORKDIR /root/fflas-ffpack-2.5.0

RUN ./configure --with-blas-libs="-L/opt/OpenBLAS/lib -lopenblas" --with-blas-cflags="-I/opt/OpenBLAS/include" && \
    make && make install

WORKDIR /root/linbox-1.7.0
RUN ./configure && make && make install

WORKDIR /root/nauty27r3
RUN ./configure && make && ln -s /root/nauty27r3/genbgL /usr/local/bin

RUN pip install tqdm pandas

# clone into GH git 
WORKDIR /root
RUN git clone https://github.com/sibrun/GH.git

# compile special rank program
WORKDIR /root/linbox-1.7.0/examples
RUN mv rank.C rank.C.old && \
    cp /root/GH/linbox_rheinfall_modified_rank/rank.C . && \
    make rank && \
    cp rank /root/GH/source/ && \
    cp -r .libs /root/GH/source/

CMD [ /bin/bash ]



### to build and run container (assuming working directory is where this Dockerfile is)
# docker build -t ubuntu-gh .
# docker run -ti ubuntu-gh /bin/bash

### then, to test:
# cd
# cd GH
# sage --python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_b -build_op -rank -linbox rational
 