# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:51:34 2019

@author: Alberto Vecchiato
"""

import numpy as np
from astropy import constants as const
from gsrconst import ppn_gamma, GM_sun

print(ppn_gamma)
print(GM_sun)

rPB = np.arange(7.0)+1.0
nelem = rPB.shape[0]
#ppn_gamma = np.ones(nelem)
ppn_gamma = 1.0
GM_sun = (const.G * const.M_sun).value

h00 = (ppn_gamma + 1) * GM_sun / (const.c.value ** 2 * rPB)

print(h00)
"""
hij = 2 * ppn_gamma * GM_sun / (const.c.value ** 2 * rPB) * np.identity(3)
h01 = 0.0
h02 = 0.0
h03 = 0.0
"""