###########################################################################################
# aeropt/__init__.py
#
# (c) copyright 2022- ecmwf.
#
# this software is licensed under the terms of the apache licence version 2.0
# which can be obtained at http://www.apache.org/licenses/license-2.0.
#
# in applying this licence, ecmwf does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
#
# author:
#    ramiro checa-garcia. ecmwf
#
# modifications:
#    30-Oct-2022   ramiro checa-garcia    1st. version
#
##########################################################################################

import os
import os.path
import sys
try:
    import tomllib as toml
except ModuleNotFoundError:
    import toml

import numpy as np
from datetime import datetime

from aeropt import cmd   as cmd
import aeropt.store_nc   as store
import aeropt.aer        as aer
import aeropt.ifs        as ifs
import aeropt.modes      as modes
import aeropt.process    as run
import aeropt.show       as show
import aeropt.engine     as engine
import aeropt.parsetoml  as parse
import aeropt.plt        as plt
