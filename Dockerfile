FROM python:2.7.18
WORKDIR /app

# update apt to use archive repositories as this is a very outdated image
RUN sed -i 's|deb.debian.org/debian|archive.debian.org/debian|g' /etc/apt/sources.list \
    && sed -i 's|security.debian.org/debian-security|archive.debian.org/debian-security|g' /etc/apt/sources.list
RUN apt update
RUN apt install -y libboost-all-dev gfortran libopenblas-dev liblapack-dev cython libtool swig
RUN pip install --upgrade pip
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

# install xraylib
RUN wget https://github.com/tschoonj/xraylib/releases/download/xraylib-3.3.0/xraylib-3.3.0.tar.gz
RUN tar -xf xraylib-3.3.0.tar.gz
RUN (cd xraylib-3.3.0 && ./configure && make && make install)
# no idea why this is needed, but it seems to fix the xraylib install
RUN apt install -y vim
RUN printf "%s\n" 'wq' > cmds.vim && vim -s cmds.vim -es /usr/local/lib/python2.7/site-packages/xraylib.py

COPY . .
RUN git submodule update --init --recursive
RUN pip install -e .
