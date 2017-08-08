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
import getpass
import logging
import os

import six
from openbabel import obErrorLog, obError, obWarning, obInfo, obDebug
from sqlalchemy.event import listens_for
from sqlalchemy.exc import DisconnectionError

__all__ = ['Level', 'log', 'prj_dir', 'db_url']

TRAVIS = os.environ.get("TRAVIS", False)


logger = logging.getLogger(__name__)


def _add_process_guards(engine):
    """Add multiprocessing guards.

    Forces a connection to be reconnected if it is detected
    as having been shared to a sub-process.

    """

    @listens_for(engine, "connect")
    def connect(connection, connection_record):
        connection_record.info['pid'] = os.getpid()

    @listens_for(engine, "checkout")
    def checkout(connection, connection_record, connection_proxy):
        pid = os.getpid()
        if connection_record.info['pid'] != pid:
            connection_record.connection = connection_proxy.connection = None
            raise DisconnectionError(
                "Connection record belongs to pid %s, "
                "attempting to check out in pid %s" %
                (connection_record.info['pid'], pid))


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

# TODO: specify database connection configuration

try:
    logger.debug("Looking for setup.cfg")
    with open('setup.cfg') as file:
        config.readfp(file)
except IOError as e:
    logger.debug("Not available %s" % str(e))
    config = default

if isinstance(config, dict) or six.PY3:
    prj_dir = os.path.abspath(config['marsi']['prj_dir'])
    db_name = config['marsi']['db_name']

else:
    prj_dir = os.path.abspath(config.get('marsi', 'prj_dir'))
    db_name = config.get('marsi', 'db_name')

logger.info("Working dir %s" % prj_dir)

if not os.path.exists(prj_dir):
    os.mkdir(prj_dir)


log.level = Level.ERROR

db_config = {}

try:
    db_engine = db_config.get('db_engine', "postgresql")

    if TRAVIS:
        username = db_config.get('db_user', 'postgres')
        password = None
        host = None
        port = None
    else:  # TODO: needs documentation
        username = db_config.get('db_user', getpass.getuser())
        password = db_config.get('db_pass', None)
        host = db_config.get("db_host", "localhost")
        port = int(db_config.get("db_port", 5432))

    if password is None:
        user_access = username
    else:
        user_access = "%s:%s" % (username, password)

    if port is None:
        host_port = host
    else:
        host_port = "%s:%i" % (host, port)

    db_url = "%s://%s@%s/%s" % (db_engine, user_access, host_port, db_name)

except Exception as e:
    logger.error(e)
    default_session = None
    engine = None
    db_url = None
    logger.warning("You are running MARSI without a database connection. \n"
                   "You will not be able to use many functionalities including: \n"
                   "1. Query the metabolite database\n"
                   "2. Search for metabolite analogs")
else:
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker

    engine = create_engine(db_url, client_encoding='utf8', pool_size=10)

    _add_process_guards(engine)

    Session = sessionmaker(engine)
    default_session = Session()
