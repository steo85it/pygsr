from gsr_obs_eq import obs_eq


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