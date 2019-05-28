import numpy as np
import pandas as pd
from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord

from gsrconst import ppn_gamma
from gsropt import debug
from gsr_util import norm, normalize, eulerAnglesToRotationMatrix


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
        cosphi = self.set_cosPhi(source)
        print("cosPhi = ")
        print(cosphi)
        # exit()
        #
        # print("Processing partials")
        # cstar = SkyCoord(ra=source.cat.ra.values * u.rad, dec=source.cat.dec.values * u.rad,
        #              distance = 1/source.cat.par.values * u.au,frame='icrs')
        #              # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')
        # print(cstar)
        # print(cstar.to(u.m))
        # exit()
        # cstar.representation = 'cartesian'
        # cstar = np.hstack([cstar.x, cstar.y, cstar.z])
        # rA = np.linalg.norm(cstar)
        # src = source.cat
        # part = {'': 1,
        #         'pi': -u.AU/src.par * cstar/rA,
        #         'ra': rA*(-np.sin(src.ra)*np.cos(src.dec),
        #                   np.cos(src.ra) * np.cos(src.dec),
        #                  0),
        #         'dec': rA*(-np.cos(src.ra)*np.sin(src.dec),
        #                   np.sin(src.ra) * np.sin(src.dec),
        #                   np.cos(src.dec))}
        # dt = self.auxdf.epo_x
        # part['mua'] = dt * part['ra']
        # part['mud'] = dt * part['dec']
        #
        # for p in part:
        #     self.set_partials(source,p)


    def set_partials(self,partial):

        print("Now processing partial w.r.t ", partial.key)
        dxA = partial.value

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar()

        dNAB = - (NAB) / rAB #- (NAB x dxA x NAB) / rAB
        drAB = - np.einsum('ik,jk->j', NAB, dxA)
        drPA = np.einsum('ik,jk->j', NPA, dxA)
        dNPA = (dxA-NPA*drPA)/rPA

        dk = dNAB - (ppn_gamma+1)*GM_sun/(const.c.value**2 * rPB) / (1+np.einsum('ik,jk->j', NPA, NPB)) *(
            - np.einsum('ik,jk->j', dNPA, NPB)/(1+np.einsum('ik,jk->j', NPA, NPB))*(
            NAB*rAB/rPA - NAB*(1+rPB/rPA)   ) +
            NPB/rPA**2 * (rPA*drAB - rAB*drPA) - dNAB*(1+rPB/rPA) + NAB*drPA*rPB/rPA**2
        )

        met_tensor = self.set_metric()
        h00 = met_tensor[:,0,0]
        #TODO just ok because h0i = 0
        h01 = h02 = h03 = met_tensor[:,0,1]

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
        if debug:
            khat = NAB

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

        print("parallax")
        print(source.cat.par.values)

        cstar = SkyCoord(ra=source.cat.ra.values * u.rad, dec=source.cat.dec.values * u.rad,
                     distance = 1/(206265*source.cat.par.values) * u.pc,frame='icrs')
                     # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')

        if debug:
            print(df.columns)
            print("cstar radec =", cstar)

        cstar.representation = 'cartesian'
        rstar = np.hstack([cstar.x.to(u.m),cstar.y.to(u.m),cstar.z.to(u.m)])

        if debug:
            print("cstar =", cstar)
            print("cstar =", rstar, np.linalg.norm(rstar))
            print("cstar normalized =", rstar/np.linalg.norm(rstar))

            print(rstar)
            print(df[['Sat_x','Sat_y','Sat_z']].values)
            print(np.subtract(rstar,
                                                                    df[['Sat_x','Sat_y','Sat_z']].values))
            print([2.266305e+19, -3.200638e+15,  2.111004e+19]/np.linalg.norm([2.266305e+19, -3.200638e+15,  2.111004e+19]))

        df['RAB_x'],df['RAB_y'],df['RAB_z'] = np.transpose(np.subtract(rstar,
                                                                    df[['Sat_x','Sat_y','Sat_z']].values))
        # rAB = np.linalg.norm(RAB,axis=1)
        df['RPB_x'], df['RPB_y'], df['RPB_z'] = df.Sun_x - df.Sat_x, df.Sun_y - df.Sat_y,df.Sun_z - df.Sat_z

        df['RPA_x'], df['RPA_y'], df['RPA_z'] = np.transpose(np.subtract(df[['Sun_x','Sun_y','Sun_z']].values,
                                                            rstar))

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

        if debug:
            print("NAB")
            print(NAB)

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

        met_tensor = self.set_metric()

        if debug:
            print("khat")
            print(khat)
            print("metric:")
            print(self.set_metric())

        h00 = met_tensor[:,0,0]

        #TODO just ok because h0i = 0
        h01 = h02 = h03 = met_tensor[:,0,1]

        E_tetrad = self.get_com_tetrad(h00, h01, h02, h03)

        if debug:
            print("E_tetrad:")
            print(E_tetrad)

        Etet_k = np.einsum('lkj,lj->lk', E_tetrad[:,:,1:], khat[:,:])

        denom = E_tetrad[:,0,0]+Etet_k[:,0]
        cosPsi = (E_tetrad[:,0,1:]+Etet_k[:,1:]) / denom.reshape((-1,1))

        if debug:
            print(E_tetrad[0,:,1:].shape,khat[0,:].shape)
            print(E_tetrad[0,:,1:])
            print(khat[0,:])
            print(Etet_k[:])
            print(denom)

        if debug:
            print("cosPsi=")
            print(cosPsi)
        # exit()

        return cosPsi

    def get_com_tetrad(self, h00, h01, h02, h03):

        if debug:
            print(self.auxdf)

        l_bst = self.get_local_frame(h00, h01, h02, h03)
        # Compute the SRS attitude matrix
        eulerAngles = self.auxdf[['angle_psi','angle_theta','angle_phi']].values
        # TODO ACHTUNG!!!
        # eulerAngles1 = eulerAngles + np.vstack([[0,0,np.deg2rad(53.25)],[0,0,np.deg2rad(53.25)]])
        # # exit()
        # rot_mat = eulerAnglesToRotationMatrix(eulerAngles1)
        # print("rot_mat +53")
        # print(rot_mat[:,0,:])
        # eulerAngles1 = eulerAngles - np.vstack([[0,0,np.deg2rad(53.25)],[0,0,np.deg2rad(53.25)]])
        # # exit()
        # rot_mat = eulerAnglesToRotationMatrix(eulerAngles1)
        # print("rot_mat -53")
        # print(rot_mat[:,0,:])
        #
        # # print(cartesian_to_spherical(rot_mat[:,0,0], rot_mat[:,1,0], rot_mat[:,2,0]))
        # #
        # # cstar = SkyCoord(x= rot_mat[0,0,0] * u.m, y=rot_mat[0,1,0] * u.m, z=rot_mat[0,2,0] * u.m,frame='icrs')
        # # print(cstar)
        # #              # pm_ra=source.cat.mu_a.values * u.mas / u.yr, pm_dec=source.cat.mu_d.values*u.mas/u.yr, frame='icrs')
        # # cstar.representation = 'spherical'
        # # print(cstar)
        rot_mat = eulerAnglesToRotationMatrix(eulerAngles)
        # print(rot_mat)

        # Compute the spatial part of the tetrad
        E_tetrad = np.einsum('ijk,ikl->ijl', rot_mat, l_bst[:,1:])

        if debug:
            print("rotmat x l_bst")
            print(rot_mat)
            print(l_bst[:, 1:])
            print(E_tetrad)
        # exit()
        # add temporal components E0i
        E_tetrad = np.concatenate([l_bst[:,:1],E_tetrad],axis=1)

        return E_tetrad

    def get_local_frame(self, h00, h01, h02, h03):

        GM_sun, NAB, NPA, NPB, rAB, rPA, rPB, beta_sat = self.get_auxvar()

        if debug:
            beta_sat = beta_sat*0
            h00 = h00 * 0

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
            print(self.auxdf)
        else:
            self.auxdf = df.copy()

        if debug:
            bas_angle = np.rad2deg(self.auxdf.angle_phip - self.auxdf.angle_phi)*2
            print("bas_angle = ",bas_angle)

        khat = self.set_khat()
        cospsi = self.set_cosPsi(khat)

        # Compute phi_calc, z_calc and kt AL and AC
        # ovviamente anche qui dovrei mettere tutto vettoriale, palle!
        cosphi = cospsi[:,0] / np.sqrt(1.0 - cospsi[:,2]**2)
        print("cos phi_calc = ")
        print(cosphi)
        print("phi_calc deg = ")
        print(np.rad2deg(np.arccos(cosphi)))
        print(np.rad2deg(np.arccos(cosphi))-53.25)
        print("eta = ")
        print(np.rad2deg(self.auxdf.eta.values))

        return cosphi

    # def plapos(self,jd0, target,center=12):
    #
    #     # perform the computation
    #     PV = peph.compute_unit(np.floor(jd0), jd0 - np.floor(jd0), target, center,Constants.UNIT_KM + Constants.UNIT_SEC)
    #     # print('PV',PV)
    #     return PV

    def proj(self, k_hat):

        return np.linalg.norm(k_hat,axis=1)

    def set_partials(self, p):
        pass