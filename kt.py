#!/usr/bin/env python3
# ----------------------------------
# pygsr.py
#
# Description:
#
# ----------------------------------------------------
# Author: Stefano Bertone
# Created: 18-Feb-2019
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
import pandas as pd
from calcephpy import *

# import time
# from scipy import interpolate
# import pickle
# import re
# import matplotlib.pyplot as plt
# import subprocess

class star:

    def __init__(self, id, data_in):

        self.id = id
        self.data_df = data_in
        self.obs_eq = None

    def set_obs_eq(self):

        self.obs_eq = obs_eq(self.data_df)
        self.obs_eq.setup()

class obs_eq:

    def __init__(self, data_in):

        self.df = data_in
        self.A = None
        self.b = None
        self.x = None
        self.weights = None

    def setup(self):

        self.plapos(2.451545E+06,10)
        self.set_kt()

        part = ['']
        for p in part:
            self.set_partials(p)


    def set_kt(self):

        df = self.df.copy()

        kGM_TTF = df.factor.values*(df.term1.values-df.term2.values)
        # print(df)
        k_hat = np.vstack(df['gaia2star'].values + kGM_TTF)

        self.df.loc[:,'kt'] = self.proj(k_hat)

        print(self.df)
        # kGM_TTF = np.vstack(df['kGM_TTF'].values)

        # print(self.properties)
        #
        # kGM_TTF = self.properties['factor']*(np.array(self.properties['term1'])-np.array(self.properties['term2']))
        #
        # self.k_hat = self.properties['gaia2star'] + kGM_TTF

    def plapos(self,jd0, target,center=12):

        # perform the computation
        PV = peph.compute_unit(jd0, 0, target, center,Constants.UNIT_KM + Constants.UNIT_SEC)
        print('PV',PV)

    def proj(self, k_hat):

        return np.linalg.norm(k_hat,axis=1)

    def set_partials(self, p):
        pass


def read_parse(infil):

    df = pd.read_csv(infil, sep='\t', header=0)
    # df['orbID'] = infil.split('.')[0][-10:]
    # self.name = df['orbID'].unique().squeeze()

    # strip and lower case all column names
    df.columns = ['gtime','rstarUpd','gaiapos','factor','term1','term2','gaia2star']

    df['rstarUpd'] = df.rstarUpd.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    df['gaiapos'] = df.gaiapos.apply(lambda x: np.array([float(y) for y in (x.strip('()').split(','))]))
    df['term1'] = df.term1.apply(lambda x: np.array([float(y) for y in (x.strip('()').split(','))]))
    df['term2'] = df.term2.apply(lambda x: np.array([float(y) for y in (x.strip('()').split(','))]))
    df['gaia2star'] = df.gaia2star.apply(lambda x: np.array([float(y) for y in (x.strip('()').split(','))]))
    df['star_id'] = '123456789'
    # df[['rstarUpd_x','rstarUpd_y','rstarUpd_z']] = pd.DataFrame(df.rstarUpd.values.tolist(), index= df.index)
    # df[['gaiapos_x','gaiapos_y','gaiapos_z']] = pd.DataFrame(df.gaiapos.values.tolist(), index= df.index)
    # df[['rstarUpd_x','rstarUpd_y','rstarUpd_z']] = pd.DataFrame(df.rstarUpd.values.tolist(), index= df.index)
    # df[['rstarUpd_x','rstarUpd_y','rstarUpd_z']] = pd.DataFrame(df.rstarUpd.values.tolist(), index= df.index)


    # df['kGM_TTF'] = df.factor.values*(df.term1.values-df.term2.values)
    # print(df)
    #
    # k_hat = np.vstack(df['gaia2star'].values + df['kGM_TTF'].values)
    # kGM_TTF = np.vstack(df['kGM_TTF'].values)
    #
    # print(k_hat)
    # print(kGM_TTF)

    # exit()
    # df.columns = df.columns.str.lower()

    # only select the required data (column)
    return df

if __name__ == '__main__':

    infil = 'auxdir/input.in'

    # open the ephemeris file
    peph = CalcephBin.open("auxdir/inpop17a_TDB_m100_p100_tt.dat")

    df = read_parse(infil)
    stars = [star(x,df.loc[df['star_id']==x]) for x in df.star_id.unique()]

    [s.set_obs_eq() for s in stars]


    # df = pd.read_csv(infil, sep='\t', header=0)
    # df.columns = ['gtime', 'rstarUpd', 'gaiapos', 'factor', 'term1', 'term2', 'gaia2star']
    #
    # df['rstarUpd'] = df.rstarUpd.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    # df['gaiapos'] = df.gaiapos.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    # df['term1'] = df.term1.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    # df['term2'] = df.term2.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    # df['gaia2star'] = df.gaia2star.apply(lambda x: [float(y) for y in (x.strip('()').split(','))])
    #
    # print(df.to_dict('records'))
    # # strip and lower case all column names
    # obs = [obs(i) for i in df.to_dict('records')]
    #
    # print('obs_epo= ',obs[0].epoch)
    # obs[0].set_k()
    # print('k_hat = ',obs[0].k_hat)

    exit()