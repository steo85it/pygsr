#!/usr/bin/env python3
# ----------------------------------

##############################################
from gsrconst import rad2arcsec

unix = 1
projv = 'b'
debug = 0
num_parts = 1

#relat = 1

# out and aux
# if (unix == 0):
#     outdir = '/att/nobackup/sberton2/MLA/out/'
#     auxdir = '/att/nobackup/sberton2/MLA/aux/'
# else:
# #    outdir = '/home/sberton2/Works/NASA/Mercury_tides/out/'
#     outdir = '/home/sberton2/Works/NASA/LOLA/out/'
#     auxdir = '/home/sberton2/Works/NASA/LOLA/aux/'


class opt:

    def __init__(self):

        self.relat = None
        self.sigma_pert = None
        self.meas_err_sigma = None
        # self.debug = None

    def set_relat(self,k):
        self.relat = k

    # def set_debug(self,k):
    #     self.debug = k

    def set_cat_err(self,k):
        self.sigma_pert = k

    def set_meas_err_sigma(self,k):
        self.meas_err_sigma = k

    # def get_debug(self):
    #     return self.debug