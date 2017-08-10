FROM ubuntu:latest

RUN export PS1="$USER@marsi"

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y postgresql wget build-essential

RUN wget http://repo.continuum.io/miniconda/Miniconda-3.5.5-Linux-x86_64.sh -O miniconda.sh;
RUN bash miniconda.sh -b -p $HOME/conda
RUN export PATH="$HOME/conda/bin:$PATH"
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN conda create -q -n marsi27 python=2.7 pip cmake
RUN conda create -q -n marsi35 python=3.5 pip cmake

COPY . /marsi

RUN source activate marsi27
RUN conda install -q -c rdkit boost rdkit=2016.09.4
RUN conda install -q -c openbabel openbabel
RUN pip install pip --upgrade

RUN pip install flake8 cython numpy scipy pyzmq pandas pytest pytest-cov pytest-benchmark swiglpk optlang
RUN cd /marsi && python setup.py install

RUN psql -c 'create database "marsi-db";' -U postgres
RUN psql -U postgres marsi-db -f tests/fixtures/marsi-db-schema.sql
RUN python bin/restore_db.py
RUN psql -d marsi-db -c 'SELECT COUNT(*) FROM metabolites;' -U postgres

RUN marsi
RUN cd /marsi & pytest tests

RUN source activate marsi35
RUN conda install -q -c rdkit boost rdkit=2016.09.4
RUN conda install -q -c openbabel openbabel
RUN pip install pip --upgrade

RUN pip install flake8 cython numpy scipy pyzmq pandas pytest pytest-cov pytest-benchmark swiglpk optlang
RUN cd /marsi && python setup.py install

RUN psql -c 'create database "marsi-db";' -U postgres
RUN psql -U postgres marsi-db -f tests/fixtures/marsi-db-schema.sql
RUN python bin/restore_db.py
RUN psql -d marsi-db -c 'SELECT COUNT(*) FROM metabolites;' -U postgres

RUN marsi
RUN cd /marsi & pytest tests
