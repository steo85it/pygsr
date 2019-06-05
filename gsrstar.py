import os
import pickle

from gsr_obs_eq import obs_eq
import numpy as np
import random
import pandas as pd

from gsrconst import rad2arcsec
from gsropt import debug


class star:

    def __init__(self, id=0, cat=0,opt=0):

        self.id = id
        self.cat = cat
        self.cat_orig = cat
        self.pert = None
        self.obs_df = None
        self.eph_df = None
        self.att_df = None
        self.obs_eq = None
        self.num_part = None
        self.opt = opt

        # print("Processing star #", self.id)

    def perturb(self,sigma_pert):

        if debug:
            print("Perturbing star #",self.id)
            print(self.cat)
            print(sigma_pert)

        random.seed(self.id)
        self.pert = dict(zip(sigma_pert.keys(),[random.gauss(0, parstd) for parstd in sigma_pert.values()]))
        self.cat = self.cat_orig.reset_index(drop=True).add(pd.DataFrame(self.pert,index=[0]),fill_value=0)

        if debug:
            print(pd.DataFrame(self.pert,index=[0]))
            print(self.cat)
            print("star #",self.id," perturbed")


    def numeric_partials(self):

        pert = {'ra':0.01/rad2arcsec,'dec':0.01/rad2arcsec,'par':1.e-4/rad2arcsec,'mu_a':0.01/rad2arcsec,'mu_d':0.01/rad2arcsec}

        num_part = []
        for par,pval in pert.items():

            # perturb catalog (+)
            # self.pert = dict(zip(sigma_pert.keys(),p))
            if debug:
                print("Processing", par,':',pval)

            # restore initial values
            self.cat = self.cat_orig

            self.cat = self.cat.reset_index(drop=True).add(pd.DataFrame(dict(zip([par],[pval])),index=[0]),fill_value=0)
            # print("compute +",pval,self.cat)
            # compute obs with + perturbation
            self.set_obs_eq()
            b_plus = self.obs_eq.b
            # restore initial values
            self.cat = self.cat_orig

            # perturb catalog (-)
            self.cat = self.cat.reset_index(drop=True).subtract(pd.DataFrame(dict(zip([par],[pval])),index=[0]),fill_value=0)
            # print("compute -",pval,self.cat)
            # compute obs with - perturbation
            self.set_obs_eq()
            b_minus = self.obs_eq.b
            # restore initial values
            self.cat = self.cat_orig

            if debug:
                print("numerical partials test")
                print(b_plus-b_minus,2*abs(pval))
                print((b_plus - b_minus)/(2*abs(pval)))

            num_part.append((b_plus - b_minus)/(2*abs(pval)))

        if debug:
            print("Resulting Num Parts:")
            print(np.column_stack(num_part))

        self.num_part = pd.DataFrame(np.column_stack(num_part),columns=['ra','dec','par','mu_a','mu_d']).fillna(value=0)

    def set_obs_eq(self,simobs=False):

        self.obs_eq = obs_eq(self,simobs,self.opt)
        self.obs_eq.setup(self)

        if simobs:
            if debug:
                print("phi_obs gsrstar", self.obs_eq.auxdf.phi_obs)
                self.obs_df.eta[:2] = self.obs_eq.auxdf.phi_obs.values
            else:
                if len(self.obs_eq.auxdf) != len(self.obs_df) and debug:
                    print("Inconsistent number of obs in star#", self.id,len(self.obs_eq.auxdf), len(self.obs_df))
                # print(self.obs_eq.auxdf.phi_obs)
                self.obs_df = self.obs_df[:len(self.obs_eq.auxdf.phi_obs)]

                meas_err_sigma = self.opt.meas_err_sigma # in radians
                if meas_err_sigma != 0:
                    measurm_err = [random.gauss(0, meas_err_sigma) for i in self.obs_eq.auxdf.phi_obs.values]
                    self.obs_df['eta'] = self.obs_eq.auxdf.phi_obs.values + np.rad2deg(measurm_err)
                else:
                    self.obs_df['eta'] = self.obs_eq.auxdf.phi_obs.values

                if debug:
                    print(self.obs_eq.auxdf.phi_obs.values, np.rad2deg(measurm_err))
                    # exit()

    # save object to pkl
    def save(self, filnam):
        pklfile = open(filnam, "wb")
        pickle.dump(self, pklfile)
        pklfile.close()

    # load from file
    def load(self, filnam):

        if os.path.isfile(filnam):
            pklfile = open(filnam, 'rb')
            self = pickle.load(pklfile)
            pklfile.close()
        else:
            print("No " + filnam + " found")
            self = None
        # print('Object loaded from '+filnam)
        # print(self.ladata_df)
        # print(self.MGRx.tck)

        return self
