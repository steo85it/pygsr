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

from astropy import units as u
from gsr_util import read_parse, read_parse_b
from gsropt import unix, projv, debug
from gsrstar import star

warnings.filterwarnings("ignore", category=RuntimeWarning)

if __name__ == '__main__':

    if projv == 'b':

        infils = glob.glob('auxdir/plan_b_full/*.txt')
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
        # print(dfs['eph'].filter(regex="_[x,y,z]"))
        dfs['eph'][dfs['eph'].filter(regex="_[x,y,z]").columns.values] = dfs['eph'].filter(regex="_[x,y,z]").apply(lambda x: (x.values * u.au).to(u.m).value)
        dfs['eph'][dfs['eph'].filter(regex="_[vx,vy,vz]").columns.values] = dfs['eph'].filter(regex="_[vx,vy,vz]").apply(lambda x: (x.values * u.cm/u.s).to(u.m/u.s).value)

        stars = [star(x,
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
        [s.set_obs_eq() for s in stars if len(s.eph_df > 0)]

    else:
        from calcephpy import CalcephBin

        infil = 'auxdir/input.in'

        # open the ephemeris file
        peph = CalcephBin.open("auxdir/inpop17a_TDB_m100_p100_tt.dat")

        df = read_parse(infil)
        stars = [star(x, df.loc[df['star_id'] == x]) for x in df.star_id.unique()]

        [s.set_obs_eq() for s in stars]

    exit()