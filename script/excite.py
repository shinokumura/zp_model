import numpy as np
import pandas as pd
import math

try:
    from .pair import FFPair
    from .func import n_ind_sep_en, normalization
    from .mass import mass_excess
    from .ya import _gen_h_mass_range
    from .level import ld_shell_correction, kckAsymptoticLevelDensity, ld_data
except:
    from pair import FFPair
    from func import n_ind_sep_en, normalization
    from mass import mass_excess
    from ya import _gen_h_mass_range
    from level import ld_shell_correction, kckAsymptoticLevelDensity, ld_data


"""
function to devide excitation energy to two fragments
"""


def _hf3d_tke_systematics(ffmassrange):
    """
    currently only U-236 (compound)
    """
    a = 335.284
    b = 1.17404
    c = 0.187596
    d = 69.0834
    return list(((a - b * x) * (1 - c * math.exp(-((x - 118) ** 2) / d)) for x in ffmassrange))



def fragment_excitation_energy(params_dict, tke_ah_dict = None, width_ah_dict = None, fractional_df = pd.DataFrame()):
    """
    convert tke to tex
    txe(zl,al,zh,ah) = Einc + Bn + [M(zcn,acn) - Mn(zl,al) - Mn(zh,ah)] x c^2 - TKE(zl,al,zh,ah)
    but generate only for heavy fragment data
    """

    if not params_dict:
        from .params import ACN, ZCN, EIN, RT, TKE_WIDTH_CONST, TKE_WIDTH, TKE_WIDTH_GIVEN

    else:
        ZCN = params_dict["ZCN"]
        ACN = params_dict["ACN"]
        EIN = params_dict["EIN"]
        RT = params_dict["RT"]
        TKE_WIDTH_CONST = params_dict["TKE_WIDTH_CONST"]
        TKE_WIDTH = params_dict["TKE_WIDTH"]
        TKE_WIDTH_GIVEN = params_dict["TKE_WIDTH_GIVEN"] # must be True

    fragment_dict = {}

    mass_ex_tn = mass_excess(ZCN, ACN - 1)  # only the case of the first chance fission
    mass_ex_cn = mass_excess(ZCN, ACN)

    ## separation energy of fissioning system
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)
    mcn = mass_ex_cn + sep_cn + EIN


    # should loop over ah
    for ah in tke_ah_dict.keys():
        s = 0.0
        n = 0.0
        ss = 0.0
        
        ## here similler to ffpyield.cpp: FFPTXEDistribution
        ## filter with mass number a and make a list of charge number
        zh_charge_list = fractional_df.loc[fractional_df.mass_h == float(ah)]["charge_h"].to_list()
        zl_charge_list = [ZCN - zh for zh in zh_charge_list]

        ## charge product for the Coulomb repulsion
        pair_charge_prod = [
            zh_charge_list[i] * zl_charge_list[i] for i in range(len(zh_charge_list))
        ]

        ## for the pair (but loop over heavy fragments)
        s = sum(pair_charge_prod)
        n = len(pair_charge_prod)
        if s < 1E-30:
            continue
        else:
            ss = tke_ah_dict[ah] * n / s

        ## assume TKE for each separation is proportional to the Coulomb repulsion
        for zh in zh_charge_list:

            ex = 0.0
            pair = FFPair(ZCN, ACN, zh, ah, mcn, tke_ah_dict, fractional_df)
            
            q = pair.q_value()
            ek = zh * pair.get_pair_z() * ss
            ex = q - ek


            if ex < 0.5:
                # need to modify dataframe
                pair.update_pair_yield(0.0)
                continue
            
            ## Excitation energy divide into two fragments
            ratio = _energy_devide(pair, ex, RT)
            # print(pair.get_pair_z(), pair.get_pair_a(), zh, ah, ek, ex, pair.q_value(), s, ss, n)

            if ratio > 0.0:
                eH_mean = ex * (1.0 - ratio)
                eL_mean = ex * ratio


            ## Width of TXE convert from TKE width           
            if TKE_WIDTH_CONST:
                eH_width = TKE_WIDTH
                eL_width = TKE_WIDTH
            

            elif TKE_WIDTH_GIVEN and width_ah_dict is not None:
                c = gen_tke_width(width_ah_dict[ah], eH_mean, eL_mean)
                eH_width = c * eH_mean
                eL_width = c * eL_mean


            else:
                c = gen_tke_width(TKE_WIDTH, eH_mean, eL_mean)
                eH_width = c * eH_mean
                eL_width = c * eL_mean


            if eH_mean < 0 or eL_mean < 0:
                continue

            else:
                fragment_dict[pair.k] ={
                    "al": pair.get_pair_a(),
                    "zl": pair.get_pair_z(),
                    "ffy":float(pair.ffy),
                    "ah": ah,
                    "zh": zh,
                    "ke": ek,
                    "ex": ex,
                    "ex_H": eH_mean, 
                    "ex_w_H": eH_width,
                    "ex_L": eL_mean,
                    "ex_w_L":  eL_width}

    ## normalize to 1.0 again, just in case
    ffy_l = [ ffy["ffy"] for ffy in fragment_dict.values() ]
    if sum(ffy_l) != 1.0:
        ffy_l2 = normalization(ffy_l,1.0)
        i = 0
        for key, vals in fragment_dict.items():
            fragment_dict[key]["ffy"] = ffy_l2[i]
            i += 1 
        
    df = pd.DataFrame.from_dict(fragment_dict, orient='index')

    return fragment_dict




def _energy_devide(pair, ex, RT):
    '''
    so far, RT model only
    See details in JNST 2018, 55, 9, 1009–1023 https://doi.org/10.1080/00223131.2018.1467288
    '''
    # pair

    ldpa_h = kckAsymptoticLevelDensity(pair.ah)
    ldpa_l = kckAsymptoticLevelDensity(pair.get_pair_a())

    def _rt_model():
        eps = 1.0e-6
        max_iter = 100
        x = 1.0

        ld_h = ld_data(pair.zh, pair.ah)
        ld_l = ld_data(pair.get_pair_z(), pair.get_pair_a())
        

        if any(len(l)==0 for l in (ld_h, ld_l)):
            ## if there is no data in kcksyst.dat then ratio should be 0.5 (devide by half)
            ratio = 0.5
            return ratio

        ## memo; paring: ld[2],  eshell: ld[3], lda: ld[4], temp: ld[5], e0: ld[6]

        eH = ex * 0.5
        eL = ex * 0.5

        aH = ld_shell_correction(eH - ld_h[2], ldpa_h, ld_h[3], pair.ah)
        aL = ld_shell_correction(eL - ld_l[2], ldpa_l, ld_l[3], pair.get_pair_a())

        # print("start:   ",eH, eL, aH, aL)
        p0 = aL/aH * RT * RT

        for i in range(max_iter):
            eH = (ex - (ld_l[2]+ld_h[2]))/(1.0 + p0) + ld_h[2]
            eL = ex - eH

            aH = ld_shell_correction(eH - ld_h[2], ldpa_h, ld_h[3], pair.ah)
            aL = ld_shell_correction(eL - ld_l[2], ldpa_l, ld_l[3], pair.get_pair_a())

            p1 = aL/aH * RT * RT
            x = abs(1.0-p1/p0)

            if x < eps:
                break

            p0 = p1

        ratio = eL / (eL + eH)
        # print(pair.zh, pair.ah, ratio)

        return ratio


    if RT == 1.0:
        ratio = 0.5

    elif RT != 1.0:
        ratio = _rt_model()

    else:
        ratio = 1.0

    return ratio




def gen_tke_a(ACN, MAX_A) -> dict:
    """
    assueming the tke as a function of compound mass is given in the following format
    """
    # fission fragment mass range
    ah_range = _gen_h_mass_range(ACN, MAX_A)
    # will make option from Viola (average over compound nucleus mass) systematics and analytical functions
    tke_ah = _hf3d_tke_systematics(ah_range)

    return dict(zip(ah_range, tke_ah))




def gen_tke_width(tke_width, eH_mean, eL_mean):
    return tke_width / math.sqrt(eH_mean * eH_mean + eL_mean * eL_mean)  




if __name__ == "__main__":
    gen_tke_a()
