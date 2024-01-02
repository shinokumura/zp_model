# from .params import ACN, ZCN
try:
    from .mass import mass_excess
except:
    from mass import mass_excess


class FFPair:
    def __init__(self, ZCN, ACN, zh, ah, mcn, tke_ah_dict, fractional_df):
        ## Compund
        self.ZCN = ZCN
        self.ACN = ACN
        ## Fission fragment pair
        self.zh = int(zh)
        self.ah = int(ah)
        self.mcn = mcn
        self.tke_ah_dict = tke_ah_dict
        self.fractional_df = fractional_df

        if zh > 0:
            self.zl = self.get_pair_z()
        if ah > 0:
            self.al = self.get_pair_a()
            
        self.k = self.get_pair_id()
        self.ffy = self.get_ffyield()


    def __repr__(self) -> str:
        return "Nuclide: " + self.zh + "-" + self.ah


    def get_pair_id(self):
        return self.fractional_df.loc[
            (self.fractional_df.charge_h == float(self.zh))
            & (self.fractional_df.mass_h == float(self.ah))
        ]["k"].to_string(index=False)


    def get_pair_z(self):
        return int(self.ZCN - self.zh)


    def get_pair_a(self):
        return int(self.ACN - self.ah)


    def get_ave_tke(self):
        return self.tke_ah_dict[self.ah]


    def get_ffyield(self):
        return self.fractional_df.loc[
            (self.fractional_df.charge_h == float(self.zh))
            & (self.fractional_df.mass_h == float(self.ah))
        ]["pair_yield"].to_string(index=False)


    def update_pair_yield(self,y):
        self.fractional_df.loc[(self.fractional_df.charge_h == float(self.zh)) & (self.fractional_df.mass_h == float(self.ah)), "pair_yield"] = y
        

    def pair_mass_excess(self):
        return mass_excess(self.zl, self.al) + mass_excess(self.zh, self.ah)


    def q_value(self):
        # print(self.mcn, self.pair_mass_excess())
        return self.mcn - self.pair_mass_excess()

