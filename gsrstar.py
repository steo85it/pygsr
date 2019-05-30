from gsr_obs_eq import obs_eq
import numpy as np
import random
import pandas as pd

from gsropt import debug


class star:

    def __init__(self, id, cat):

        self.id = id
        self.cat = cat
        self.pert = None
        self.obs_df = None
        self.eph_df = None
        self.att_df = None
        self.obs_eq = None

    def perturb(self,sigma_pert):

        print(self.cat)
        print(sigma_pert)

        random.seed(self.id)
        self.pert = dict(zip(sigma_pert.keys(),[random.gauss(0, parstd) for parstd in sigma_pert.values()]))
        print(pd.DataFrame(self.pert,index=[0]))

        self.cat = self.cat.add(pd.DataFrame(self.pert,index=[0]),fill_value=0)
        print(self.cat)

    def set_obs_eq(self,simobs=False):

        self.obs_eq = obs_eq(self,simobs)
        self.obs_eq.setup(self)

        if simobs:
            print("phi_obs gsrstar", self.obs_eq.auxdf.phi_obs)
            if debug:
                self.obs_df.eta[:2] = self.obs_eq.auxdf.phi_obs.values
            else:
                self.obs_df.eta[:] = self.obs_eq.auxdf.phi_obs.values
