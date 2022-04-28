FROM ubuntu:focal

ENV TZ=Etc/UTC
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update
RUN apt upgrade -y
RUN apt install -y postgresql wget build-essential cmake python3-pip swig 

# Build and install RDKit
RUN apt install -y libboost-dev libboost-python-dev libboost-numpy-dev libboost-system-dev libboost-math-dev libboost-iostreams-dev libboost-serialization-dev libboost-program-options-dev libfreetype-dev libfreetype6 libcairo2 libcairo2-dev libeigen3-dev python3-numpy
RUN wget https://github.com/rdkit/rdkit/archive/Release_2022_03_1.tar.gz
RUN tar -xzf Release_2022_03_1.tar.gz
ENV RDBASE=$HOME/RDKit
RUN mv rdkit-Release_2022_03_1 $RDBASE
ENV LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=$RDBASE:$PYTHONPATH
RUN mkdir $RDBASE/build
WORKDIR $RDBASE/build
 
RUN cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_AVALON_SUPPORT=ON -DRDK_BUILD_CAIRO_SUPPORT=ON -DRDK_BUILD_PYTHON_WRAPPERS=ON ..
RUN make
RUN make install

# Build and install openbabel
WORKDIR /root
RUN wget https://github.com/openbabel/openbabel/releases/download/openbabel-3-1-1/openbabel-3.1.1-source.tar.bz2
RUN tar -xjf openbabel-3.1.1-source.tar.bz2
RUN mkdir openbabel-3.1.1/build
WORKDIR /root/openbabel-3.1.1/build
COPY ./openbabel-python.i /root/openbabel-3.1.1/scripts/
RUN pip install eigen
RUN cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON ..
RUN make
RUN make install

RUN mkdir /marsi
WORKDIR /marsi
COPY ./requirements.txt /marsi/
RUN pip install -r requirements.txt
RUN pip install --no-deps pytest==7.1.2

COPY . /marsi
RUN python3 setup.py install
RUN cp -R build/lib.linux-x86_64-3.*/* ./
RUN chmod +x ./entrypoint.sh

# RUN service postgresql start
# RUN su postgres

# RUN psql -c 'create database "marsi-db";' -U postgres
# RUN psql -U postgres marsi-db -f tests/fixtures/marsi-db-schema.sql
# RUN python3 bin/restore_db.py
# RUN psql -d marsi-db -c 'SELECT COUNT(*) FROM metabolites;' -U postgres

# RUN su root

# RUN pytest --disable-warnings tests/test_io_db.py || true

# RUN pytest --disable-warnings tests/test_bigg_api.py || true
# RUN pytest --disable-warnings tests/test_chemistry.py || true
# RUN pytest --disable-warnings tests/test_cobra.py || true
# RUN pytest --disable-warnings tests/test_data_retrieval.py || true
# RUN pytest --disable-warnings tests/test_optmet_design.py || true
# RUN pytest --disable-warnings tests/test_post_processing.py || true
# RUN pytest --disable-warnings tests/test_targets.py || true
# RUN pytest --disable-warnings tests/test_utils.py || true

ENTRYPOINT ./entrypoint.sh
