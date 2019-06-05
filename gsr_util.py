import os
import pickle

import numpy as np
import pandas as pd
from astropy import units as u
from scipy.sparse.linalg import lsqr

from gsropt import debug, num_parts, unix

def norm(param, df):

    cols = [param+'_x',param+'_y',param+'_z']
    tmp = np.vstack(df[cols].values)
    return np.linalg.norm(tmp,axis=1)


def normalize(param, df):
    cols = [param+'_x',param+'_y',param+'_z']
    tmp = np.vstack(df[cols].values)
    norm = np.reshape(np.linalg.norm(tmp,axis=1),(-1,1))
    return tmp/norm


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
    tmp = np.einsum('ijk,ikl->ijl', M2, M1)

    if debug:
        print("M1")
        print(M1)
        print("M2")
        print(M2)
        print("M3")
        print(M3)
        print("M3xM2xM1")
        print(np.einsum('ijk,ikl->ijl', M3, tmp))

    return np.einsum('ijk,ikl->ijl', M3, tmp)


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


def process(stars,options):
    # TODO change to len(s.obs_df > 0) when using full dataset
    [s.set_obs_eq(simobs=True) for s in stars if len(s.obs_df > 0)]
    if debug:
        print("Simulated observations (phi_obs)")
        print([s.obs_eq.auxdf.phi_obs for s in stars])
    # exit()
    # check numerical ders
    if num_parts:
        [s.numeric_partials() for s in stars]
    # define and apply perturbation to catalog : sigma pars in as, as/y
    # sigma_pert = {'ra': 1.e-2 / rad2arcsec, 'dec': 1.e-2 / rad2arcsec, 'par': 1.e-2 / rad2arcsec, 'mu_a': 1.e-4/ rad2arcsec,
    #               'mu_d': 1.e-4 / rad2arcsec}
    [s.perturb(sigma_pert=options.sigma_pert) for s in stars if len(s.obs_df > 0)]
    # TODO change to len(s.obs_df > 0) when using full dataset
    [s.set_obs_eq() for s in stars if len(s.obs_df > 0)]


def solve_star(s):
    print("Solution for star #", s.id)
    print("Perts to retrieve :", s.pert)
    # print(s.obs_eq.b)
    if debug:
        print("first analyt parts :", s.obs_eq.A.loc[:1, :])
        if num_parts:
            print("first numeric parts :", s.num_part.loc[:1, :])
    # solution based on analytical partials
    subset_anal_part = s.obs_eq.A.loc[:, list(s.pert.keys())].values
    x = np.linalg.lstsq(subset_anal_part * 100., s.obs_eq.b, rcond=None)
    print("analyt sol (rad, parts*100) : ", x)
    s.obs_eq.x = x

    if num_parts:
        # solution based on numerical partials
        subset_num_part = s.num_part.loc[:, list(s.pert.keys())].values
        x = np.linalg.lstsq(subset_num_part, s.obs_eq.b, rcond=None)
        print("num sol (rad):", x)
        s.obs_eq.x = x

        if debug:
            x = lsqr(subset_num_part, s.obs_eq.b, show=True)
            print("num sol scipy (rad):", x)


def load_data(infils):
    # global dfs
    cols = {
        'Ephem': ['frameID', 'epo', 'Sun_x', 'Sun_y', 'Sun_z', 'Sat_x', 'Sat_y', 'Sat_z', 'Sat_vx', 'Sat_vy', 'Sat_vz'],
        'Catalog': ['sourceID', 'ra', 'dec', 'par', 'mu_a', 'mu_d'],
        'Observ': ['sourceID', 'frameID', 'fovID', 'eta', 'zeta'],
        'Scan': ['frameID', 'epo', 'angle_psi', 'angle_theta', 'angle_phi', 'angle_phip', 'angle_phif']}
    dfnam = ['cat', 'eph', 'obs', 'scan']
    dfs = []
    for f in infils:
        if unix:
            dfs.append(read_parse_b(f, cols=cols[f.split('/')[-1].split('.')[0]]))
        else:
            dfs.append(read_parse_b(f, cols=cols[f.split('\\')[-1].split('.')[0]]))
    dfs = dict(zip(dfnam, dfs))
    # update ephemeris units to m, m/s
    dfs['eph'][dfs['eph'].filter(regex="_[x,y,z]").columns.values] = dfs['eph'].filter(regex="_[x,y,z]").apply(
        lambda x: (x.values * u.au).to(u.m).value)
    dfs['eph'][dfs['eph'].filter(regex="_v[x,y,z]").columns.values] = dfs['eph'].filter(regex="_v[x,y,z]").apply(
        lambda x: (x.values * u.cm / u.s).to(u.m / u.s).value)

    return dfs