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
import pandas as pd
import time
import warnings
import numpy as np

from gsr_util import solve_star, load_data, process
from gsrconst import rad2arcsec
from gsropt import projv, opt, debug
import gsrstar

warnings.filterwarnings("ignore", category=RuntimeWarning)

if __name__ == '__main__':

    ##############################################
    # launch program and clock
    # -----------------------------
    start = time.time()

    if projv == 'b':

        options = opt()
        options.set_cat_err({'ra': 1e-200 / rad2arcsec,
                      'dec': 1e-200 / rad2arcsec})  # ,'par':1e-2/rad2arcsec,'mu_a':1e-4/rad2arcsec,'mu_d':1e-4/rad2arcsec}
        options.set_meas_err_sigma(0) #1e-3 / rad2arcsec)
        # options.set_debug(0)
        options.set_relat(1)

        infils = np.sort(glob.glob('auxdir/plan_b/*.txt'))
        print(infils)
        dfs = load_data(infils)

        #         print([(s_id, len(dfs['obs'].loc[dfs['obs'].sourceID == s_id]),
        #           dfs['cat'].loc[dfs['cat'].sourceID == s_id][['ra', 'dec']]) for s_id in dfs['cat'].sourceID.unique()])
        nobs_per_star = pd.DataFrame([(s_id,len(dfs['obs'].loc[dfs['obs'].sourceID==s_id])) for s_id in dfs['cat'].sourceID.unique()], columns=["id", "num_obs"])
        print(nobs_per_star.max())
        print(nobs_per_star.sort_values(by=["num_obs"]))
            # print((s_id,len(dfs['obs'].loc[dfs['obs'].sourceID==s_id])))

        print("Processing stars #", dfs['cat'].sourceID.unique()[811:812])
        # exit()
        chosen_stars = dfs['cat'].sourceID.unique()[811:812]

        stars = [gsrstar.star(x,
                      cat= dfs['cat'].loc[dfs['cat'].sourceID==x],opt=options)
                 for x in chosen_stars]

        if debug:
            stars = stars[:]

        for s in stars:
            setattr(s,'obs_df',
                    dfs['obs'].loc[dfs['obs'].sourceID==s.id] )
            setattr(s,'eph_df',
                    dfs['eph'].loc[dfs['eph'].frameID.isin(s.obs_df.frameID)] )
            setattr(s,'att_df',
                    dfs['scan'].loc[dfs['scan'].frameID.isin(s.obs_df.frameID)] )
            s.obs_df.reset_index(inplace=True)
            s.eph_df.reset_index(inplace=True)
            s.att_df.reset_index(inplace=True)

        # s.save('test.pkl')

        # s1 = star()

        process(stars,options)

        for s in stars:

            _ = solve_star(s)

    # stop clock and print runtime
    # -----------------------------
    end = time.time()
    print('----- Runtime = ' + str(end - start) + ' sec -----' + str((end - start) / 60.) + ' min -----')
