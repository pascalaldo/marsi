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
import logging

import os
import six

from openbabel import obErrorLog, obError, obWarning, obInfo, obDebug

__all__ = ['Level', 'log', 'prj_dir']


logger = logging.getLogger(__name__)


class Level:
    ERROR = obError
    WARNING = obWarning
    INFO = obInfo
    DEBUG = obDebug


class LogConf(object):
    __level_to_logger = {
        Level.ERROR: logging.ERROR,
        Level.WARNING: logging.WARNING,
        Level.DEBUG: logging.DEBUG,
        Level.INFO: logging.INFO
    }

    def __init__(self):
        self._app_logger = logging.getLogger('marsi.*')
        self._level = Level.WARNING

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, level):
        obErrorLog.SetOutputLevel(level)
        self._app_logger.setLevel(self.__level_to_logger[level])
        self._level = level


log = LogConf()

log.level = Level.INFO

config = six.moves.configparser.ConfigParser()
default = {'marsi': {'prj_dir': "%s/.marsi" % os.getenv('HOME'), 'db_name': 'marsi-db'}}

try:
    logger.debug("Looking for setup.cfg")
    with open('setup.cfg') as file:
        config.readfp(file)
except IOError as e:
    logger.debug("Not available %s" % str(e))
    config = default

prj_dir = os.path.abspath(config['marsi']['prj_dir'])
logger.info("Working dir %s" % prj_dir)

db_name = config['marsi']['db_name']

if not os.path.exists(prj_dir):
    os.mkdir(prj_dir)


log.level = Level.ERROR



