FROM biosustain/cameo:latest

RUN export PS1="$USER@marsi"

COPY . /marsi

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 0C49F3730359A14518585931BC711F9BA15703C6
RUN apt-get update
RUN apt-get install -y libopenbabel-dev libeigen3-dev cython mongodb python-rdkit
RUN pip install --upgrade pip
RUN pip install Cython numpy

RUN pip install gevent

RUN python --version

RUN cd /marsi && python setup.py install

RUN pytest tests/

RUN marsi init download && marsi init build-database

ENTRYPOINT ['marsi']