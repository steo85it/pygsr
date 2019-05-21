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

from gsrconst import ppn_gamma
from gsropt import unix, projv

warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
import pandas as pd
#from calcephpy import *
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

# Calculates Rotation Matrix given euler angles.
def eulerAnglesToRotationMatrix(theta) :

    # Tested for one single row of theta=[[1,0,0]], gives the same result
    # thetab = theta[0,:]
    # R_x = np.array([[1,         0,                  0                   ],
    #                 [0,         np.cos(thetab[0]), -np.sin(thetab[0]) ],
    #                 [0,         np.sin(thetab[0]), np.cos(thetab[0])  ]
    #                 ])
    # R_y = np.array([[np.cos(thetab[1]),    0,      np.sin(thetab[1])  ],
    #                 [0,                     1,      0                   ],
    #                 [-np.sin(thetab[1]),   0,      np.cos(thetab[1])  ]
    #                 ])
    # R_z = np.array([[np.cos(thetab[2]),    -np.sin(thetab[2]),    0],
    #                 [np.sin(thetab[2]),    np.cos(thetab[2]),     0],
    #                 [0,                     0,                      1]
    #                 ])
    # R = np.dot(R_z, np.dot( R_y, R_x ))
    # print(R)
    
    nrows = theta.shape[0]

    alpha = theta[:,0]
    M1 = np.reshape(np.hstack([np.column_stack((np.cos(alpha), -np.sin(alpha), np.zeros(nrows))),
                               np.column_stack((np.sin(alpha), np.cos(alpha), np.zeros(nrows))),
                               np.column_stack((np.zeros(nrows), np.zeros(nrows), np.ones(nrows)))]), (-1, 3, 3))

    alpha = theta[:,1]
    M2 = np.reshape(np.hstack([np.column_stack((np.ones(nrows), np.zeros(nrows), np.zeros(nrows))),
                               np.column_stack((np.zeros(nrows), np.cos(alpha), -np.sin(alpha))),
                               np.column_stack((np.zeros(nrows), np.sin(alpha), np.cos(alpha)))]), (-1, 3, 3))

    alpha = theta[:,2]
    M3 = np.reshape(np.hstack([np.column_stack((np.cos(alpha), -np.sin(alpha), np.zeros(nrows))),
                               np.column_stack((np.sin(alpha), np.cos(alpha), np.zeros(nrows))),
                               np.column_stack((np.zeros(nrows), np.zeros(nrows), np.ones(nrows)))]), (-1, 3, 3))
    #tsipm = 0  # mtimesx(M1,mtimesx(M2,M3));

    # Combine rotations
    tmp = np.einsum('ijk,ikl->ijl', M2, M3)

    return np.einsum('ijk,ikl->ijl', M1, tmp)

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
        self.set_cosPhi(source)

        part = ['']
        for p in part:
            self.set_partials(p)


    def set_partials(self,parameter):

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB = self.get_auxvar(df)


        dk = dNAB - (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) / (1+np.einsum('ik,jk->j', NPA, NPB)) *(
            - np.einsum('ik,jk->j', dNPA, NPB)/(1+np.einsum('ik,jk->j', NPA, NPB))*(
            NAB*RAB/RPA - NAB*(1+RPB/RPA)   ) +
            NPB/RPA**2 * (RPA*dRAB - RAB*dRPA) - dNAB*(1+rPB/rPA) + NAB*dRPA*rPB/rPA**2
        )

        h00, h01, h02, h03 = self.set_metric()
        E_tetrad = self.get_com_tetrad(h00, h01, h02, h03)

        Etet_dk = np.einsum('lij,lj->li', E_tetrad[:,:,1:], dk)
        Etet_k = np.einsum('lij,lj->li', E_tetrad[:,:,1:], k)
        dcosPsi = - ( Etet_dk[:,1:]*(E_tetrad[:,0,0]+Etet_k[:,0]) - (E_tetrad[:,1:,0]+Etet_k[:,1:])*Etet_dk[:,0] ) / (
                        E_tetrad[:,0,0]+Etet_k[:,0]   )**2

        cosPsi = self.set_cosPsi()

        dcosPhi = dcosPsi[:,1]/np.sqrt(1-cosPsi[:,2]**2) + cosPsi[:,0]*cosPsi[:,2]*dcosPsi[:,2]/(1-cosPsi[:,3]**2)**3./2.


        # khat = NAB - np.dot( (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) * (1+np.einsum('ik,jk->j', NPA, NPB))**(-1) , (rAB/rPA * NPB - (1+rPB/rPA)*NAB))

    def set_khat(self, source):

        df = pd.merge(source.obs_df, source.eph_df, on='frameID')
        df = pd.merge(df, source.att_df, on='frameID')

        print(source.cat)
        cstar = SkyCoord(ra=source.cat.ra.values * u.rad, dec=source.cat.dec.values * u.rad,
                     distance = 1/source.cat.par.values * u.au,frame='icrs')
                     # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')
        cstar.representation = 'cartesian'

        print(df.columns)
        print(cstar)

        df['RAB_x'],df['RAB_y'],df['RAB_z'] = np.transpose(np.subtract(np.hstack([cstar.x,cstar.y,cstar.z]),
                                                                    df[['Sat_x','Sat_y','Sat_z']].values))
        # rAB = np.linalg.norm(RAB,axis=1)
        df['RPB_x'], df['RPB_y'], df['RPB_z'] = df.Sun_x - df.Sat_x, df.Sun_y - df.Sat_y,df.Sun_z - df.Sat_z

        df['RPA_x'], df['RPA_y'], df['RPA_z'] = np.transpose(np.subtract(df[['Sun_x','Sun_y','Sun_z']].values,
                                                            np.hstack([cstar.x,cstar.y,cstar.z])))

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB = self.get_auxvar(df)

        khat = NAB - np.dot( (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) * (1+np.einsum('ik,jk->j', NPA, NPB))**(-1) , (rAB/rPA * NPB - (1+rPB/rPA)*NAB))
        print(khat)


        # exit()
        
        # self.df.loc[:,'kt'] = self.proj(khat)

        # print(self.df)
        # kGM_TTF = np.vstack(df['kGM_TTF'].values)

        # print(self.properties)
        #
        # kGM_TTF = self.properties['factor']*(np.array(self.properties['term1'])-np.array(self.properties['term2']))
        #
        # self.khat = self.properties['gaia2star'] + kGM_TTF

    def get_auxvar(self, df):
        NAB = normalize('RAB', df)
        rAB = norm('RAB', df)
        NPB = normalize('RPB', df)
        rPB = norm('RPB', df)
        NPA = normalize('RPA', df)
        rPA = norm('RPA', df)
        GM_sun = (const.G * const.M_sun).value
        beta_sq = (df['Sat_vx'] ** 2 + df['Sat_vy'] ** 2 + df['Sat_vz'] ** 2) / (((const.c).value) ** 2)
        beta_x = df['Sat_vx'] / ((const.c).value)
        beta_y = df['Sat_vy'] / ((const.c).value)
        beta_z = df['Sat_vz'] / ((const.c).value)
        return GM_sun, NAB, NPA, NPB, rAB, rPA, rPB

    def set_metric(self):
        # Compute the metric components
        h00 = (ppn_gamma + 1) * GM_sun / (const.c.value ** 2 * rPB)
        # TODO
        hij = 2 * ppn_gamma * GM_sun / (const.c.value ** 2 * rPB) * np.identity(3)
        h01 = 0.0
        h02 = 0.0
        h03 = 0.0
        return h00, h01, h02, h03

    def set_cosPsi(self):

        h00, h01, h02, h03 = self.set_metric()

        E = self.get_com_tetrad(h00, h01, h02, h03)
        # Compute the direction cosines
        # Ho preso la formula da Bertone et al. A&A608A, 2017, equazione (6);
        # al momento è solo un riferimento, sorry Il problema è che devo prendere solo le parti spaziali
        # delle tetradi per fare il prodotto con hij e khat e sto remando
#        denom = E00 + Ej0 * hij * ki
#        cosPsi = -(E[0]i+Eji * hjl * kl) / denom

        return cosPsi

    def get_com_tetrad(self, h00, h01, h02, h03):
        l_bst = self.get_local_frame(h00, h01, h02, h03)
        # Compute the SRS attitude matrix
        # TODO: define the Euler angles in the Data Model and substitute a,b,c with their values
        rot_mat = eulerAnglesToRotationMatrix([a, b, c])
        # Compute the AL observable (phi_calc)
        E_tetrad = np.einsum('ijk,ikl->ijl', rot_mat, l_bst)
        return E_tetrad

    def get_local_frame(self, h00, h01, h02, h03):
        # Compute the local BCRS and the bootsed tetrads (use hij)
        l_bcrs = np.array([[h01, 1.0 - h00 / 2, 0.0, 0.0],
                           [h02, 0.0, 1.0 - h00 / 2, 0.0],
                           [h03, 0.0, 0.0, 1.0 - h00 / 2]])
        fact = (1.0 + 3 * h00 / 2 + beta_sq / 2)
        l_bst = l_bcrs + np.array([[beta_x * fact, (beta_x ** 2) / 2, beta_x * beta_y / 2, beta_x * beta_z / 2],
                                   [beta_y * fact, beta_x * beta_y / 2, (beta_y ** 2) / 2, beta_y * beta_z / 2],
                                   [beta_z * fact, beta_x * beta_z / 2, beta_y * beta_z / 2, (beta_z ** 2) / 2]])
        return l_bst

    def set_cosPhi(self, source):

        self.set_khat(source)
        self.set_cosPsi()

        # Compute phi_calc, z_calc and kt AL and AC
        # ovviamente anche qui dovrei mettere tutto vettoriale, palle!
        phi = np.arccos(n[0] / numpy.sqrt(1.0 - n[3] ** 2))

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

    if projv == 'b':

        tmp = eulerAnglesToRotationMatrix(np.array([[1,0,0]]))
        print(tmp)

        infils = glob.glob('auxdir/plan_b/*.txt')
        print(infils)
        cols = {'Ephem':['frameID','epo','Sun_x','Sun_y','Sun_z','Sat_x','Sat_y','Sat_z','Sat_vx','Sat_vy','Sat_vz'],
                'Catalog':['sourceID','ra','dec','par','mu_a','mu_d'],
                'Observ':['sourceID','frameID','fovID','eta','zeta'],
                'Scan':['frameID','epo','SRS_X_x','SRS_X_y','SRS_X_z','SRS_Y_x','SRS_Y_y','SRS_Y_z','SRS_Z_x','SRS_Z_y','SRS_Z_z',
                        'FOV_m_x','FOV_m_y','FOV_m_z','FOV_p_x','FOV_p_y','FOV_p_z'] }

        dfnam = ['cat','eph','obs','scan']
        dfs = []
        for f in infils:
            if unix:
                dfs.append(read_parse_b(f, cols=cols[f.split('/')[-1].split('.')[0]]))
            else:
                dfs.append(read_parse_b(f, cols=cols[f.split('\\')[-1].split('.')[0]]))

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
        from calcephpy import CalcephBin

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