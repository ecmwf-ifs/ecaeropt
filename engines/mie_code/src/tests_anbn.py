
#  +----------------------------------------------------------------------------------------+
#  | engines/mie_code/src/tests_anby.py                                                     |
#  |                                                                                        |
#  | (C) Copyright 2022- ECMWF.                                                             |
#  |                                                                                        |
#  | This software is licensed under the terms of the Apache Licence Version 2.0            |
#  | which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                   |
#  |                                                                                        |
#  | In applying this licence, ECMWF does not waive the privileges and immunities           |
#  | granted to it by virtue of its status as an intergovernmental organisation             |
#  | nor does it submit to any jurisdiction.                                                |
#  |                                                                                        |
#  |                                                                                        |
#  | Author:                                                                                |
#  |    Ramiro Checa-Garcia. ECMWF                                                          |
#  |                                                                                        |
#  | Modifications:                                                                         |
#  |    25-Dec-2022   Ramiro Checa-Garcia    1st. version                                   |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+

import numpy as np
from mie_code import mie as mie

# Testing using the new structured Mie code

x  = 300.0*2.0*np.pi/375.0
#wl = 0.6328
m  = 1.77 - 0.63j

a  = np.zeros(3000, dtype='D')
b  = np.zeros(3000, dtype='D')

Nmax=30

x, a1, a2 = mie.mie_core_anbn(x,m, Nmax, a, b)



