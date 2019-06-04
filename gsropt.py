#!/usr/bin/env python3
# ----------------------------------

##############################################
from gsrconst import rad2arcsec

unix = 1
projv = 'b'
debug = 0
num_parts = 1

relat = 1

# out and aux
# if (unix == 0):
#     outdir = '/att/nobackup/sberton2/MLA/out/'
#     auxdir = '/att/nobackup/sberton2/MLA/aux/'
# else:
# #    outdir = '/home/sberton2/Works/NASA/Mercury_tides/out/'
#     outdir = '/home/sberton2/Works/NASA/LOLA/out/'
#     auxdir = '/home/sberton2/Works/NASA/LOLA/aux/'
sigma_pert = {'ra':1e-2/rad2arcsec,'dec':1e-2/rad2arcsec} #,'par':1e-2/rad2arcsec,'mu_a':1e-4/rad2arcsec,'mu_d':1e-4/rad2arcsec}
meas_err_sigma = 1e-3/rad2arcsec