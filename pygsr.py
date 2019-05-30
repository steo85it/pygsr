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
import time
import warnings
import numpy as np

from astropy import units as u
from gsr_util import read_parse, read_parse_b
from gsrconst import rad2arcsec
from gsropt import unix, projv, debug
import gsrstar

warnings.filterwarnings("ignore", category=RuntimeWarning)

if __name__ == '__main__':

    ##############################################
    # launch program and clock
    # -----------------------------
    start = time.time()

    if projv == 'b':

        infils = np.sort(glob.glob('auxdir/plan_b/*.txt'))
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
        dfs['eph'][dfs['eph'].filter(regex="_v[x,y,z]").columns.values] = dfs['eph'].filter(regex="_v[x,y,z]").apply(lambda x: (x.values * u.cm/u.s).to(u.m/u.s).value)

        stars = [gsrstar.star(x,
                      cat= dfs['cat'].loc[dfs['cat'].sourceID==x])
                 for x in dfs['cat'].sourceID.unique()][:1]

        if debug:
            stars = stars[:]

        for s in stars:
            setattr(s,'obs_df',
                    dfs['obs'].loc[dfs['obs'].sourceID==s.id] )
            setattr(s,'eph_df',
                    dfs['eph'].loc[dfs['eph'].frameID.isin(s.obs_df.frameID)] )
            setattr(s,'att_df',
                    dfs['scan'].loc[dfs['scan'].frameID.isin(s.obs_df.frameID)] )

        # TODO change to len(s.obs_df > 0) when using full dataset
        [s.set_obs_eq(simobs=True) for s in stars if len(s.eph_df > 0)]
        print(s.obs_eq.auxdf.phi_obs.values)

        # sigma pars in as, as/y
        sigma_pert = {'ra': 1.e-2 / rad2arcsec, 'dec': 1.e-2 / rad2arcsec, 'par': 1.e-2 / rad2arcsec, 'mu_a': 1.e-4/ rad2arcsec,
                      'mu_d': 1.e-4 / rad2arcsec}
        # sigma_pert = {'ra':10/rad2arcsec,'dec':10/rad2arcsec,'par':1/rad2arcsec,'mu_a':0.1/rad2arcsec,'mu_d':0.1/rad2arcsec}
        [s.perturb(sigma_pert=sigma_pert) for s in stars if len(s.eph_df > 0)]

        # TODO change to len(s.obs_df > 0) when using full dataset
        [s.set_obs_eq() for s in stars if len(s.eph_df > 0)]

        exit()

    # stop clock and print runtime
    # -----------------------------
    end = time.time()
    print('----- Runtime = ' + str(end - start) + ' sec -----' + str((end - start) / 60.) + ' min -----')
