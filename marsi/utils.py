# Copyright 2016 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import gzip
import logging
import os
import pickle
import re
import shutil
import time

import numpy as np
from IProgress import ProgressBar, Percentage
from cameo import fba
from cameo.flux_analysis.analysis import n_carbon
from cobra.core.reaction import Reaction

from marsi import config

__all__ = ['data_dir', 'log_dir', 'pickle_large', 'unpickle_large', 'frange', 'src_dir']

data_dir = os.path.join(config.prj_dir, "data")
models_dir = os.path.join(config.prj_dir, "models")
log_dir = os.path.join(config.prj_dir, "log")
src_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")

INCHI_KEY_TYPE = np.dtype("a27")

BIOMASS_RE = re.compile("biomass", re.IGNORECASE)

MAX_BYTES = 2 ** 31 - 1

logger = logging.getLogger(__name__)


def pickle_large(obj, file_path, progress=False):
    with open(file_path, 'wb') as model_handler:
        bytes_out = pickle.dumps(obj)
        output_size = len(bytes_out)
        if progress:
            pbar = ProgressBar(maxval=output_size, widgets=["Writing ", Percentage()])
            for idx in pbar(range(0, output_size, MAX_BYTES)):
                model_handler.write(bytes_out[idx:idx + MAX_BYTES])
        else:
            for idx in range(0, output_size, MAX_BYTES):
                model_handler.write(bytes_out[idx:idx + MAX_BYTES])


def unpickle_large(file_path, progress=False):
    input_size = os.path.getsize(file_path)
    logger.debug("Input size: %f bytes" % input_size)
    with open(file_path, 'rb') as file_handler:
        bytes_in = bytearray(0)

        if progress:
            pbar = ProgressBar(maxval=input_size, widgets=["Loading ", Percentage()])
            for _ in pbar(range(0, input_size, MAX_BYTES)):
                bytes_in += file_handler.read(MAX_BYTES)
        else:
            for _ in range(0, input_size, MAX_BYTES):
                bytes_in += file_handler.read(MAX_BYTES)

    return pickle.loads(bytes_in)


def frange(start, stop=None, steps=10):
    """
    Float range generator.

    Generates *steps* equally separated between *start* and *stop*.
    If *stop* is None, the values are between 0 and *start*

    Parameters
    ----------
    start : float
        The initial value.
    stop : float
        The final value.
    steps : int
        Number of values to generate.

    Returns
    -------
    generator
        A generator that yields float.
    """
    if stop is None:
        stop = start
        start = 0

    step_size = (stop - start) / steps
    logger.debug("Step size %f" % step_size)
    for i in range(steps):
        logger.debug("Iteration %i: %f" % (i + 1, i * step_size))
        yield start + i * step_size


def unique(l):
    """
    Removes repeated values from a list.

    Parameters
    ----------
    l: list

    Returns
    -------
    list
        The same list with only unique values.
    """
    s = set()
    n = 0
    for x in l:
        if x not in s:
            s.add(x)
            l[n] = x
            n += 1
    del l[n:]


def timing(debug=False):  # pragma: no cover
    def function_wrapper(func):
        if debug:
            def debug_wrap_func(*args, **kwargs):
                start = time.time()
                ret = func(*args, **kwargs)
                stop = time.time()
                if config.log.level >= config.Level.DEBUG:
                    print('%s function took %0.3f ms' % (func.__name__, (stop - start) * 1000.0))
                return ret
            return debug_wrap_func
        else:
            def wrap_func(*args, **kwargs):
                start = time.time()
                ret = func(*args, **kwargs)
                stop = time.time()
                print('%s function took %0.3f ms' % (func.__name__, (stop - start) * 1000.0))
                return ret
            return wrap_func
    return function_wrapper


def default_carbon_sources(model):
    solution = fba(model)
    carbon_sources = []

    for ex in model.exchanges:
        assert isinstance(ex, Reaction)
        if ex.lower_bound < 0 and solution[ex.id] < 0 < n_carbon(ex):
            logger.debug("Found carbon source: %s")
            carbon_sources.append(ex)

    return carbon_sources


def gunzip(file):
    assert file[-3:] == ".gz"
    in_name = file
    out_name = file[0:-3]
    with gzip.open(in_name, 'rb') as f_in, open(out_name, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def search_metabolites(model, species_id, ignore_external=True):
    if ignore_external:
        return model.metabolites.query(lambda mid: mid[:-2] == species_id and mid[-2:] != "_e", attribute='id')
    else:
        return model.metabolites.query(lambda mid: mid[:-2] == species_id, attribute='id')
