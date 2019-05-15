#!/usr/bin/env python3
# ----------------------------------
# pygsr.py
#
# Description:
#
# ----------------------------------------------------
# Author: Stefano Bertone
# Created: 18-Feb-2019
import glob
import warnings
import weakref

from astropy import units as u
from astropy.coordinates import SkyCoord
# from astropy.constants import au, GM_sun
# from astropy.constants import c as clight
from astropy import constants as const

warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
import pandas as pd
from calcephpy import *
from sklearn import preprocessing as prp
# import time
# from scipy import interpolate
# import pickle
# import re
# import matplotlib.pyplot as plt
# import subprocess

class star:

    def __init__(self, id, cat):

        self.id = id
        self.cat = cat
        self.obs_df = None
        self.eph_df = None
        self.att_df = None
        self.obs_eq = None

    def set_obs_eq(self):

        self.obs_eq = obs_eq(self)
        self.obs_eq.setup(self)


def norm(param, df):

    cols = [param+'_x',param+'_y',param+'_z']
    tmp = np.vstack(df[cols].values)
    return np.linalg.norm(tmp)

def normalize(param, df):
    cols = [param+'_x',param+'_y',param+'_z']
    tmp = np.vstack(df[cols].values)
    norm = np.linalg.norm(tmp)
    return tmp/norm


class obs_eq:

    def __init__(self, source):

        # self.source = weakref.ref(source)
        self.A = None
        self.b = None
        self.x = None
        self.weights = None

    def setup(self,source):

        # self.source = weakref.ref(source)
        print("Processing star #",source.id)
        self.set_kt(source)

        part = ['']
        for p in part:
            self.set_partials(p)


    def set_kt(self,source):

        df = pd.merge(source.obs_df, source.eph_df, on='frameID')
        df = pd.merge(df, source.att_df, on='frameID')

        print(source.cat)
        cstar = SkyCoord(ra=source.cat.ra.values * u.rad, dec=source.cat.dec.values * u.rad,
                     distance = 1/source.cat.par.values * u.au,frame='icrs')
                     # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')
        cstar.representation = 'cartesian'

        ppn_gamma = 1

        print(df.columns)
        print(cstar)

        df['RAB_x'],df['RAB_y'],df['RAB_z'] = np.transpose(np.subtract(np.hstack([cstar.x,cstar.y,cstar.z]),
                                                                    df[['Sat_x','Sat_y','Sat_z']].values))
        # rAB = np.linalg.norm(RAB,axis=1)
        df['RPB_x'], df['RPB_y'], df['RPB_z'] = df.Sun_x - df.Sat_x, df.Sun_y - df.Sat_y,df.Sun_z - df.Sat_z

        df['RPA_x'], df['RPA_y'], df['RPA_z'] = np.transpose(np.subtract(df[['Sun_x','Sun_y','Sun_z']].values,
                                                            np.hstack([cstar.x,cstar.y,cstar.z])))
        NAB = normalize('RAB',df)
        rAB = norm('RAB',df)

        NPB = normalize('RPB',df)
        rPB = norm('RPB',df)

        NPA = normalize('RPA',df)
        rPA = norm('RPA',df)

        GM_sun = (const.G*const.M_sun).value

        k_hat = NAB - np.dot( (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) * (1+np.einsum('ik,jk->j', NPA, NPB))**(-1) , (rAB/rPA * NPB - (1+rPB/rPA)*NAB))
        print(k_hat)

        # self.df.loc[:,'kt'] = self.proj(k_hat)

        # print(self.df)
        # kGM_TTF = np.vstack(df['kGM_TTF'].values)

        # print(self.properties)
        #
        # kGM_TTF = self.properties['factor']*(np.array(self.properties['term1'])-np.array(self.properties['term2']))
        #
        # self.k_hat = self.properties['gaia2star'] + kGM_TTF

    def plapos(self,jd0, target,center=12):

        # perform the computation
        PV = peph.compute_unit(np.floor(jd0), jd0 - np.floor(jd0), target, center,Constants.UNIT_KM + Constants.UNIT_SEC)
        # print('PV',PV)
        return PV

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

def read_parse_b(infil,cols=[]):
    df = pd.read_csv(infil, delim_whitespace=True, header=None)
    if len(df.columns) == len(cols):
        df.columns = cols
    return df

if __name__ == '__main__':

    if True:

        infils = glob.glob('auxdir/plan_b/*.txt')
        print(infils)
        cols = {'Ephem':['frameID','epo','Sun_x','Sun_y','Sun_z','Sat_x','Sat_y','Sat_z','Sat_vx','Sat_vy','Sat_vz'],
                'Catalog':['sourceID','ra','dec','par','mu_a','mu_d'],
                'Observ':['sourceID','frameID','fovID','eta','zeta'],
                'Scan':['frameID','epo','SRS_X_x','SRS_X_y','SRS_X_z','SRS_Y_x','SRS_Y_y','SRS_Y_z','SRS_Z_x','SRS_Z_y','SRS_Z_z',
                        'FOV_m_x','FOV_m_y','FOV_m_z','FOV_p_x','FOV_p_y','FOV_p_z'] }

        dfnam = ['cat','obs','scan','eph']
        dfs = []
        for f in infils:
            dfs.append(read_parse_b(f, cols=cols[f.split('/')[-1].split('.')[0]]))

        dfs = dict(zip(dfnam, dfs))

        stars = [star(x,
                      cat= dfs['cat'].loc[dfs['cat'].sourceID==x])
                 for x in dfs['cat'].sourceID.unique()]

        for s in stars:
            setattr(s,'obs_df',
                    dfs['obs'].loc[dfs['obs'].sourceID==s.id] )
            setattr(s,'eph_df',
                    dfs['eph'].loc[dfs['eph'].frameID.isin(s.obs_df.frameID)] )
            setattr(s,'att_df',
                    dfs['scan'].loc[dfs['scan'].frameID.isin(s.obs_df.frameID)] )

        # TODO change to len(s.obs_df > 0) when using full dataset
        [s.set_obs_eq() for s in stars if len(s.eph_df > 0)]

    else:

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