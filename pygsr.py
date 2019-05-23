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
from gsropt import unix, projv, debug

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
    return np.linalg.norm(tmp,axis=1)

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
#        for p in part:
#            self.set_partials(source,p)


    def set_partials(self,source,parameter):

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar(source)

        dNAB = 0

        dk = dNAB - (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) / (1+np.einsum('ik,jk->j', NPA, NPB)) *(
            - np.einsum('ik,jk->j', dNPA, NPB)/(1+np.einsum('ik,jk->j', NPA, NPB))*(
            NAB*rAB/rPA - NAB*(1+rPB/rPA)   ) +
            NPB/rPA**2 * (rPA*drAB - rAB*drPA) - dNAB*(1+rPB/rPA) + NAB*drPA*rPB/rPA**2
        )

        h00, h01, h02, h03 = self.set_metric(source)
        print(h00)

        E_tetrad = self.get_com_tetrad(h00, h01, h02, h03)

        Etet_dk = np.einsum('lij,lj->li', E_tetrad[:,:,1:], dk)
        Etet_k = np.einsum('lij,lj->li', E_tetrad[:,:,1:], k)
        dcosPsi = - ( Etet_dk[:,1:]*(E_tetrad[:,0,0]+Etet_k[:,0]) - (E_tetrad[:,1:,0]+Etet_k[:,1:])*Etet_dk[:,0] ) / (
                        E_tetrad[:,0,0]+Etet_k[:,0]   )**2

        cosPsi = self.set_cosPsi()

        dcosPhi = dcosPsi[:,1]/np.sqrt(1-cosPsi[:,2]**2) + cosPsi[:,0]*cosPsi[:,2]*dcosPsi[:,2]/(1-cosPsi[:,3]**2)**3./2.


        # khat = NAB - np.dot( (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) * (1+np.einsum('ik,jk->j', NPA, NPB))**(-1) , (rAB/rPA * NPB - (1+rPB/rPA)*NAB))

    def set_khat(self):

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar()

        khat = NAB - np.dot( (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) * (1+np.einsum('ik,jk->j', NPA, NPB))**(-1) ,
                             (np.reshape(rAB/rPA,(-1,1)) * NPB - (1+np.reshape(rPB/rPA,(-1,1)))*NAB))
        # khat = NAB

        return khat


        # exit()
        
        # self.df.loc[:,'kt'] = self.proj(khat)

        # print(self.df)
        # kGM_TTF = np.vstack(df['kGM_TTF'].values)

        # print(self.properties)
        #
        # kGM_TTF = self.properties['factor']*(np.array(self.properties['term1'])-np.array(self.properties['term2']))
        #
        # self.khat = self.properties['gaia2star'] + kGM_TTF

    def set_auxdf(self, source, df):

        print(df.columns)
        cstar = SkyCoord(ra=source.cat.ra.values * u.rad, dec=source.cat.dec.values * u.rad,
                     distance = 1/source.cat.par.values * u.au,frame='icrs')
                     # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')
        cstar.representation = 'cartesian'
        # print(cstar)

        df['RAB_x'],df['RAB_y'],df['RAB_z'] = np.transpose(np.subtract(np.hstack([cstar.x,cstar.y,cstar.z]),
                                                                    df[['Sat_x','Sat_y','Sat_z']].values))
        # rAB = np.linalg.norm(RAB,axis=1)
        df['RPB_x'], df['RPB_y'], df['RPB_z'] = df.Sun_x - df.Sat_x, df.Sun_y - df.Sat_y,df.Sun_z - df.Sat_z

        df['RPA_x'], df['RPA_y'], df['RPA_z'] = np.transpose(np.subtract(df[['Sun_x','Sun_y','Sun_z']].values,
                                                            np.hstack([cstar.x,cstar.y,cstar.z])))

    def get_auxvar(self):

        df = self.auxdf.copy()

        GM_sun = (const.G * const.M_sun).value
        NAB = normalize('RAB', df)
        rAB = norm('RAB', df)
        NPB = normalize('RPB', df)
        rPB = norm('RPB', df)
        NPA = normalize('RPA', df)
        rPA = norm('RPA', df)
        beta_sat = df.filter(regex='Sat_v.*').values / ((const.c).value)
#        beta_sq = np.linalg.norm(beta_sat,axis=1)**2

        return GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat

    def set_metric(self):
        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar()
        # Compute the metric components
        h00 = (ppn_gamma + 1) * GM_sun / (const.c.value ** 2 * rPB)

        nelem = len(rPB)
        h0i = np.zeros(3*nelem).reshape((nelem,3))

        id3d = np.tile(np.identity(3), (nelem, 1)).reshape(nelem,3,3)
        hij = np.reshape(-h00,(-1,1,1)) * id3d

        # some help vectors
        hlp = np.concatenate([np.reshape(h00,(-1,1)),h0i],axis=1)
        hlp2 = np.concatenate([h0i.reshape((-1,3,1)),hij],axis=2)

        return np.concatenate([hlp.reshape((-1,1,4)),hlp2],axis=1)

    def set_cosPsi(self,khat):

        # khat = self.auxdf.khat
        print(khat)

        met_tensor = self.set_metric()

        if debug:
            print("metric:")
            print(self.set_metric())

        h00 = met_tensor[:,0,0]

        #TODO just ok because h0i = 0
        h01 = h02 = h03 = met_tensor[:,0,1]

        E_tetrad = self.get_com_tetrad(h00, h01, h02, h03)

        if debug:
            print("E_tetrad:")
            print(E_tetrad)

        Etet_k = np.einsum('lij,lj->li', E_tetrad[:,:,1:], khat)

        denom = E_tetrad[:,0,0]+Etet_k[:,0]
        cosPsi = -(E_tetrad[:,1:,0]+Etet_k[:,1:]) / denom.reshape((-1,1))

        if debug:
            print("cosPsi=")
            print(cosPsi)

        return cosPsi

    def get_com_tetrad(self, h00, h01, h02, h03):
        l_bst = self.get_local_frame(h00, h01, h02, h03)
        # Compute the SRS attitude matrix
        eulerAngles = self.auxdf[['angle_psi','angle_theta','angle_phi']].values
        rot_mat = eulerAnglesToRotationMatrix(eulerAngles)

        # Compute the spatial part of the tetrad
        E_tetrad = np.einsum('ijk,ikl->ijl', rot_mat, l_bst[:,1:])
        # add temporal components E0i
        E_tetrad = np.concatenate([l_bst[:,:1],E_tetrad],axis=1)

        return E_tetrad

    def get_local_frame(self, h00, h01, h02, h03):

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar()

        # beta_sat = beta_sat*0
        # h00 = h00 * 0

        # beta_x = beta_sat[:,0]
        # beta_y = beta_sat[:,1]
        # beta_z = beta_sat[:,2]
        beta_sq = np.linalg.norm(beta_sat,axis=1)

        # Vectorial form of l_bcrs
        # l_bcrs = np.array([[h01, 1.0 - h00 / 2, 0.0, 0.0],
        #                    [h02, 0.0, 1.0 - h00 / 2, 0.0],
        #                    [h03, 0.0, 0.0, 1.0 - h00 / 2]])
        nelem = len(h00)
        id3d = np.tile(np.identity(3), (nelem, 1)).reshape(nelem,3,3)
        hlp = id3d * (1 - h00/2).reshape(-1,1,1)
        hlp2 = np.transpose([h01,h02,h03]).reshape(nelem,-1,1)

        l_bcrs = np.concatenate([hlp2,hlp],axis=2)

        # Vectorial form of
        # l_bst = l_bcrs + np.array([[beta_x * fact, (beta_x ** 2) / 2, beta_x * beta_y / 2, beta_x * beta_z / 2],
        #                            [beta_y * fact, beta_x * beta_y / 2, (beta_y ** 2) / 2, beta_y * beta_z / 2],
        #                            [beta_z * fact, beta_x * beta_z / 2, beta_y * beta_z / 2, (beta_z ** 2) / 2]])
        fact = (1.0 + 3 * h00 / 2 + beta_sq / 2)
        hlp = np.hstack([fact.reshape(-1,1),beta_sat]).reshape((nelem,1,-1))
        l_bst = l_bcrs + np.einsum('ijk,ikl->ijl', 0.5*beta_sat.reshape((nelem,-1,1)), hlp)

        u_s = np.hstack([np.reshape(1+h00+beta_sq/2,(nelem,1)),beta_sat])
        u_s = u_s.reshape(nelem,1,-1)

        l_bst = np.concatenate([u_s,l_bst],axis=1)

        if debug:
            print("l_bcrs:")
            print(l_bcrs)
            print("l_bst:")
            print(l_bst)

        return l_bst

    def set_cosPhi(self, source):

        print("processing source", source.id," ...")

        df = pd.merge(source.obs_df, source.eph_df, on='frameID')
        df = pd.merge(df, source.att_df, on='frameID')
        self.set_auxdf(source, df)

        if debug:
            self.auxdf = df.loc[:1,:].copy()
        else:
            self.auxdf = df.copy()

        bas_angle = np.rad2deg(self.auxdf.angle_phip - self.auxdf.angle_phi)*2

        khat = self.set_khat()
        cospsi = self.set_cosPsi(khat)

        # Compute phi_calc, z_calc and kt AL and AC
        # ovviamente anche qui dovrei mettere tutto vettoriale, palle!
        cosphi = cospsi[:,0] / np.sqrt(1.0 - cospsi[:,2]**2)
        print("phi_calc = ")
        print(cosphi)
        print(np.rad2deg(np.arccos(cosphi)))
        exit()

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

        infils = glob.glob('auxdir/plan_b/*.txt')
        print(infils)
        cols = {'Ephem':['frameID','epo','Sun_x','Sun_y','Sun_z','Sat_x','Sat_y','Sat_z','Sat_vx','Sat_vy','Sat_vz'],
                'Catalog':['sourceID','ra','dec','par','mu_a','mu_d'],
                'Observ':['sourceID','frameID','fovID','eta','zeta'],
                'Scan':['frameID','epo','angle_psi','angle_theta','angle_phi','angle_phip','angle_phif'] }

        dfnam = ['cat','eph','obs','scan']
        dfs = []
        for f in infils:
            if unix:
                dfs.append(read_parse_b(f, cols=cols[f.split('/')[-1].split('.')[0]]))
            else:
                dfs.append(read_parse_b(f, cols=cols[f.split('\\')[-1].split('.')[0]]))

        dfs = dict(zip(dfnam, dfs))

        # update ephemeris units to m, m/s
        dfs['eph'][dfs['eph'].filter(regex="_[x,y,z]").columns.values] = dfs['eph'].filter(regex="_[x,y,z]").apply(lambda x: (x.values * u.au).to(u.m).value)

        stars = [star(x,
                      cat= dfs['cat'].loc[dfs['cat'].sourceID==x])
                 for x in dfs['cat'].sourceID.unique()][:1]

        if debug:
            stars = stars[:1]

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