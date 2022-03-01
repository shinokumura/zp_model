from params import ACN, ZCN
from func import mass_excess

# fractional_df = pd.DataFrame()
class FFPair:
    def __init__(self, zh, ah, mcn, tke_ah_dict, fractional_df):
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
        return int(ZCN - self.zh)


    def get_pair_a(self):
        return int(ACN - self.ah)


    def get_ave_tke(self):
        # df = pd.DataFrame(columns = ['charge_h', 'mass_h', 'pair_yield'])
        return self.tke_ah_dict[self.ah]


    def get_ffyield(self):
        # df = pd.DataFrame(columns = ['charge_h', 'mass_h', 'pair_yield'])
        return self.fractional_df.loc[
            (self.fractional_df.charge_h == float(self.zh))
            & (self.fractional_df.mass_h == float(self.ah))
        ]["pair_yield"].to_string(index=False)


    def update_pair_yield(self,y):
        self.fractional_df.loc[(self.fractional_df.charge_h == float(self.zh)) & (self.fractional_df.mass_h == float(self.ah)), "pair_yield"] = y
        

    def pair_mass_excess(self):
        return mass_excess(self.zl, self.al) + mass_excess(self.zh, self.ah)


    def q_value(self):
        return self.mcn - self.pair_mass_excess()

    # def get_paircharge_product(self):
    #     zh_charge_list = self.fractional_df.loc[self.fractional_df.mass_h == float(self.ah)][
    #         "charge_h"
    #     ].to_list()

    #     zl_charge_list = [ZCN - self.ah for self.ah in zh_charge_list]

    #     print(zh_charge_list)
    #     print(zl_charge_list)

    #     charge_prod = [zh_charge_list[i] * zl_charge_list[i] for i in range(len(zh_charge_list))]

    #     print(charge_prod)

    #     s = math.prod(charge_prod)
    #     n = len(zh_charge_list)
    #     ss = self.tke_ah_dict[self.ah] * n / s

    #     return ss
