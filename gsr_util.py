import numpy as np
import pandas as pd

from gsropt import debug


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