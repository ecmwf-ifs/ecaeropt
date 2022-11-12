

import numpy as np
from mie_BB import mie_boucher_bozzo as mieBB


# test 1

x  = 300.0*2.0*np.pi/375.0
#wl = 0.6328
m  = 1.77 - 0.63j

a  = np.zeros(3000, dtype='D')
b  = np.zeros(3000, dtype='D')

Nmax=30

x, a1, a2 = mieBB.mie_core_anbn(x,m, Nmax, a, b)
