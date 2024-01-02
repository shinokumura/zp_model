import math
import pandas as pd


# from .params import CN, ACN, ZCN, EIN
try:
    from .func import n_ind_sep_en
    from .mass import mass_excess
    from .params import VY58, VX58
except:
    from func import n_ind_sep_en
    from mass import mass_excess
    from params import VY58, VX58, VX20, VX21


pd.reset_option("display.max_columns")
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("max_colwidth", None)
pd.set_option("display.width", 1200)


#######################################################
###
###              Wahl's Y(A) model  
###
#######################################################


# Z_f = ZCN 
# A_f = ACN - 1  # since the ZF AF are for U235
# PE  = EIN
# PA  = ACN
# PZ  = ZCN     # precursor (target + projectile) values, PA; PZ, and PE,

def ya_parameter_function3(p, Z_f, PE):
    # print(p,Z_f, PE)
    p1 = p[0] + p[3] * (Z_f - 92)
    p2 = p[1] + p[4] * (Z_f - 92)
    p3 = p[2] + p[5] * (Z_f - 92)

    p_par = p1 + (p2 - p1) * (1.0 - math.exp(p3 * PE))
    return p_par

def ya_parameter_function4(p, Z_f, PE):
    p1 = p[0] + p[3] * (Z_f - 92)
    p2 = p[1] + p[4] * (Z_f - 92)

    p_par = p1 + p2 * PE
    return p_par


def chain_yield(params_dict):
    CN = params_dict["elem"] + params_dict["ACN"]
    Z_f = params_dict["ZCN"]
    A_f = params_dict["ACN"] - 1
    PE = params_dict["EIN"]
    PA =params_dict["ACN"]
    PZ = params_dict["ZCN"]

    ## take equation (3) or (4) based on the table in P.12 of LA-13928
    ## y: curve intensity, sig: sigma, delta:displacement
    sig_15  = ya_parameter_function3([  2.808,   8.685,   -0.0454,   0.372,    -0.620,   -0.0122,   63  ], Z_f, PE)
    sig_24  = ya_parameter_function4([   2.45,   0.000,      None,   0.000,     0.000,      None,   38  ], Z_f, PE)
    sig_3   = ya_parameter_function4([    8.6,   0.000,      None,   0.000,     0.000,      None,   10  ], Z_f, PE)
    sig_67  = ya_parameter_function4([   3.17,   0.000,      None,   0.303,     0.000,      None,    6  , Z_f, PE])
    delta_5 = ya_parameter_function3([  25.34,   18.55,   -0.0402,  -1.220,    -1.732,       0.0,   63  ], Z_f, PE)
    a_4     = ya_parameter_function4([ 136.66,  -0.177,      None,   0.060,    -0.038,      None,   45  ], Z_f, PE)  ## A4, at the maxima of the heavy inner{peak curves
    delta_7 = ya_parameter_function4([  30.31,     0.0,      None,     0.0,       0.0,      None,    5  ], Z_f, PE)
    y_24    = ya_parameter_function4([  43.00,   -1.91,      None,   -3.41,       0.0,      None,   47  ], Z_f, PE)

    if PE > 8.0:
        y_67    = 6.8 - (6.8/12) * (PE -8.0)
        if Z_f == 93:
            y_67   = y_67 / 2.0
        elif Z_f < 93 or PE > 20.0:
            y_67 = 0.0
    else:  
        y_67    = ya_parameter_function4([   6.80,     0.0,      None,     0.0,       0.0,      None,   10  ])


    if PE < 11.96:
        y_3 = 4.060 * math.exp(0.470 * (PE - 11.96))

    else:
        T = -0.0030 + 0.0050 * (A_f - 236.0)
        y_3 = 4.060 + 86.02 * (1.0 - math.exp(T * (PE - 11.96)))

    nt      = ya_parameter_function3([  1.563,   16.66,  -0.00804,  0.0918,       0.0,       0.0,   45  ])

    a_3 = (A_f-nt)/2.0  # PA being the sum of the projectile and target mass number
    delta_4 = a_4 - (A_f - nt)/2.0
    delta_1 = -delta_5
    delta_2 = -delta_4
    delta_6 = -delta_7
    delta_3 = 0.0



    fissile = CN
    GAUSS1 =   [1.0,  sig_15, delta_1]      # principal peak curves
    GAUSS2 =   [y_24, sig_24, delta_2]      # inner peak curves
    GAUSS3 =   [y_67, sig_67, delta_6]      # wing curves y_67 delta_7 for  Z_f>94
    GAUSSSYM = [y_3,  sig_3,  delta_3]      # central peak curves, delta_3, a_3, y_3: curve intensity


    # x_range = gen_mass_range_k(70,165)
    ###### create two hampt Y(A) from 4 Gaussian functions ######

    from ya import gen_massdist
    gen_massdist(fissile, GAUSS1, GAUSS2, GAUSSSYM)




#######################################################
###
###              Wahl's Zp model  
###
#######################################################

def zp_parameter_function(p, ZCN, ACN, EIN):
    dz = ZCN - 92.0
    x  = 0.0

    mass_ex_tn = mass_excess(ZCN, ACN - 1)  # only the case of the first chance fission
    mass_ex_cn = mass_excess(ZCN, ACN)
    
    ## separation energy of fissioning system
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)

    if(EIN < 8.0):
        ## low energty < 8 MeV
        da = ACN - 236.0
        de = sep_cn + EIN - 6.551
        x = p[0] + p[1]*dz + p[2]*da + p[3]*de + p[4]*da*da

    else:
        ## high energy upto 20 MeV
        p1 = p[0] + p[2] * dz
        p2 = p[1] + p[3] * dz
        x = p1 + (p2 - p1)*(1.0 - math.exp(-p[4]*de))

    return x



def mass_region_boundary(ZCN, ACN, delta_zmax, sl50, delta_z140, deltaz_SL):
    bundary_dict = {}

    def _amax(delta_zmax, sl50):
        ## Formula in the end pf P.23
        f1 = (250.0 - ACN) / (250.0 - 236.0)

        if(f1 < 0.0):
            f1 = 0.0
        elif(f1 > 1.0):
            f1 = 1.0

        f2 = 1.0 - f1

        ak1 = 50.0 * (ACN / ZCN) - delta_zmax / sl50
        ak2 = (50.0 - delta_zmax) * (ACN / ZCN)

        return f1 * ak1 + f2 * ak2

    a_max = _amax(delta_zmax, sl50)

    b1 = 70
    b2 = 77 + 0.036 * (ACN - 236)
    b4 = (delta_zmax - delta_z140 + a_max * sl50 + 140 * deltaz_SL) / (sl50 + deltaz_SL)
    b3 = ACN - b4
    b5 = ACN - b2
    b6 = ACN - b1
    # it seems ba <--> bb, not sure why. 
    # ba = a_max
    # bb = ACN - a_max
    bb = a_max
    ba = ACN - a_max


    bundary_dict = {"a_max": a_max, "b1": b1, "b2": b2, "b3": b3, "b4": b4, "b5": b5, "b6": b6, "ba": ba, "bb": bb }


    return bundary_dict


def wahl_parameter_calc(params_dict):
    CN = params_dict["CN"]
    ZCN = params_dict["ZCN"]
    ACN = params_dict["ACN"]
    EIN = params_dict["EIN"]

    wahl_dict_low = {}
    wahl_dict_he = {}


    ## parameters for regions near peaks (B2-3, B4-5)
    wahl_dict_low["sigma_z140"] = zp_parameter_function([  0.566 ,   0.0,   0.0064,  0.0109,      0.0], ZCN, ACN, EIN) # PF(1) = SIGZ(140) = P(2) IN ORA2, CALC2
    wahl_dict_low["delta_z140"] = zp_parameter_function([ -0.487 ,   0.0,   0.0180,     0.0, -0.00203], ZCN, ACN, EIN) # PF(2) = DELTZ(140) = -P(1) IN ORA2, CALC2
    wahl_dict_low["fz_z140"]    = zp_parameter_function([  1.207 ,   0.0,  -0.0420,     0.0,  0.0022 ], ZCN, ACN, EIN) # PF(3) = EOZ(140) = P(3) IN ORA2, CALC2
    wahl_dict_low["fn_z140"]    = zp_parameter_function([  1.076 ,   0.0,      0.0,     0.0,      0.0], ZCN, ACN, EIN) # PF(4) = EON(140) = P(4) ORA2, CALC2
    wahl_dict_low["sigz_SL"]    = zp_parameter_function([ -0.0038,   0.0,      0.0,     0.0,      0.0], ZCN, ACN, EIN) # PF(5) = SIGZSL  = P(16) IN ORA2, CALC2
    wahl_dict_low["deltaz_SL"]  = zp_parameter_function([ -0.0080,   0.0,      0.0,     0.0,      0.0], ZCN, ACN, EIN) # PF(6) = DELTZSL = -P(5) IN ORA2, CALC2
    wahl_dict_low["d_fz"]       = zp_parameter_function([  0.0030,   0.0,      0.0,     0.0,      0.0], ZCN, ACN, EIN) # PF(10) = EOZSL = P(18) IN ORA2, CALC2 - 12/30/99

    ## parameters for  regions near symmetry (B3-4)
    wahl_dict_low["sl50"]       = zp_parameter_function([  0.191,     0.0,  -0.0076,    0.0,     0.0], ZCN, ACN, EIN) # PF(7) = SL50 = P(6) IN ORA2, CALC2
    wahl_dict_low["sigma_z50"]  = zp_parameter_function([  0.356,   0.060,     0.0,     0.0,     0.0], ZCN, ACN, EIN) # PF(8) = SIG50 = P(8) IN ORA2, CALC2
    wahl_dict_low["delta_zmax"] = zp_parameter_function([  0.699,     0.0,     0.0,     0.0,     0.0], ZCN, ACN, EIN) # PF(9) = DELTZmax = P(7) IN ORA2, CALC2


    # if wahl_dict_low["sl50"] < 0.5:
    #     ## temporary solution to avoid the value gets too small
    #     ## in BeoH, it restricts to 0.0
    #     wahl_dict_low["sl50"] = 0.5

    ## parameters for wing regions (B1-2, B 5-6)
    wahl_dict_low["sigma_zwing"] = zp_parameter_function([   -0.045,  0.0094,    0.0,    0.0,     0.0  ], ZCN, ACN, EIN) # PF(11) = SIGWSL = P(21) IN ORA2, CALC2
    wahl_dict_low["delta_zwing"] = zp_parameter_function([    0.0  , -0.0045,    0.0,    0.0,     0.0  ], ZCN, ACN, EIN) # PF(12) = DELWSL = P(24) IN ORA2, CALC2
    wahl_dict_low["fz_zwing"]    = zp_parameter_function([    0.159, -0.028 ,    0.0,    0.0,     0.0  ], ZCN, ACN, EIN) # PF(13) = EOZWSL = P(22) IN ORA2, CALC2
    wahl_dict_low["fn_zwing"]    = zp_parameter_function([    0.039,     0.0,    0.0,    0.0,     0.0  ], ZCN, ACN, EIN) # PF(14) = EONWSL = P(23) IN ORA2, CALC2

    if EIN < VY58:
        wahl_dict = { key: value for key, value in wahl_dict_low.items() }

    elif VY58 < EIN < VX58:
        wahl_dict_he["sigma_z140"] = zp_parameter_function([  0.542,  1.310,   0.033,   0.0,    -0.005  ], ZCN, ACN, EIN) # PF(1) 
        wahl_dict_he["delta_z140"] = zp_parameter_function([ -0.428,    0.0,    0.0,    0.164,  -0.0116 ], ZCN, ACN, EIN) # PF(2) 
        wahl_dict_he["fz_z140"]    = zp_parameter_function([    1.0,    0.0,    0.0,    0.0,     0.0    ], ZCN, ACN, EIN) # PF(3) 
        wahl_dict_he["fn_z140"]    = zp_parameter_function([    1.0,    0.0,    0.0,    0.0,     0.0    ], ZCN, ACN, EIN) # PF(4) 
        wahl_dict_he["sigz_SL"]    = 0.0    # PF(5) 
        wahl_dict_he["deltaz_SL"]  = 0.0    # PF(6) 
        wahl_dict_he["d_fz"]       = 0.0    # PF(10)

        ## parameters for  regions near symmetry (B3-4) for HE
        wahl_dict_he["sl50"]       = wahl_dict_low["sl50"]          # PF(7)
        wahl_dict_he["sigma_z50"]  = wahl_dict_low["sigma_z140"]    # PF(8)
        wahl_dict_he["delta_zmax"] = 0.0                            # PF(9)

        ## parameters for wing regions (B1-2, B 5-6) for HE
        wahl_dict_he["sigma_zwing"] = 0.0   # PF(11)
        wahl_dict_he["delta_zwing"] = 0.0   # PF(12)
        wahl_dict_he["fz_zwing"]    = 0.0   # PF(13)
        wahl_dict_he["fn_zwing"]    = 0.0   # PF(14)

        ## Linear parameter EG functions for transition between low and high energy
        FRH = (EIN - VY58) / (VX58 - VY58)
        FRL = 1.0 - FRH

        wahl_dict = { key: FRL * l + FRH * wahl_dict_he[key] for key, l in wahl_dict_low.items() }

    else:
        ## Above 20 MeV won't work
        wahl_dict = { key: value for key, value in wahl_dict_he.items() }


    bd = mass_region_boundary(ZCN, 
                              ACN, 
                              wahl_dict["delta_zmax"], 
                              wahl_dict["sl50"], 
                              wahl_dict["delta_z140"], 
                              wahl_dict["deltaz_SL"]
                              )
    
    # print("# params", CN, ZCN, ACN, "{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}{:10.3F}".format(*wahl_dict.values()))

    wahl_dict["bd"] = bd

    return wahl_dict


def wahl_systematics(params_dict, wahl_dict, a):

    ZCN = params_dict["ZCN"]
    ACN = params_dict["ACN"]
    EIN = params_dict["EIN"]

    ap = a                  # fragment A, light
    ah = ACN - ap           # A-light since the a is always heavy fragment

    sigma_z140 =  wahl_dict["sigma_z140"]          # PF(1) 
    delta_z140 =  wahl_dict["delta_z140"]          # PF(2) 
    fz_z140  =  wahl_dict["fz_z140"]               # PF(3) 
    fn_z140  =  wahl_dict["fn_z140"]               # PF(4)
    sigz_SL =  wahl_dict["sigz_SL"]                # PF(5)
    deltaz_SL =  wahl_dict["deltaz_SL"]            # PF(6) 
    d_fz     =  wahl_dict["d_fz"]                  # PF(10) 

    sl50 = wahl_dict["sl50"]                       # PF(7) 
    sigma_z50 = wahl_dict["sigma_z50"]             # PF(8) 
    delta_zmax = wahl_dict["delta_zmax"]           # PF(9) 

    sigma_zwing  = wahl_dict["sigma_zwing"]        # PF(11)
    delta_zwing  = wahl_dict["delta_zwing"]        # PF(12) 
    fz_zwing   = wahl_dict["fz_zwing"]             # PF(13) 
    fn_zwing   = wahl_dict["fn_zwing"]             # PF(14) 

    bd = wahl_dict["bd"]


    def _far_wing():
        sigma_z = sigma_z140 + sigz_SL * (bd["b5"] - 140.0)     # SIG5
        delta_z = delta_z140 + deltaz_SL * (bd["b5"] - 140.0)   # DEL5
        # if(ap < ah):
        if (ap < ACN/2):
            delta_z *= -1.0
        fz     = fz_z140
        fn     = fn_z140
        fa = wahl_evenodd_term(ZCN, ACN, fz, fn)

        return (delta_z, sigma_z, fa)


    def _wing():
        deltaB5 = delta_z140 + deltaz_SL * (bd["b5"] - 140.0) 
        sigmaB5 = sigma_z140 + sigz_SL * (bd["b5"] - 140.0) 

        if(ap < ACN/2):
            ## light fragment, but loop is always for heavy fragments only, so nver comes here
            delta_z = - deltaB5 - delta_zwing * (bd["b5"] - ah)  # APC-B5 = ACN-AP-B5
            sigma_z =   sigmaB5 + sigma_zwing * (bd["b2"] - ap)
            fz = fz_z140 + fz_zwing * (bd["b2"] - ap)
            fn = fn_z140 + fn_zwing * (bd["b2"] - ap)
        
        else:
            ## heavy fragment
            delta_z = deltaB5   - delta_zwing * (ap - bd["b5"])
            sigma_z = sigmaB5   + sigma_zwing * (ap - bd["b5"])
            fz = fz_z140 + fz_zwing * (ap - bd["b5"])
            fn = fn_z140 + fn_zwing * (ap - bd["b5"])
            
        fa = wahl_evenodd_term(ZCN, ACN, fz, fn)

        return (delta_z, sigma_z, fa)
    

    def _peak():
        if (ap < ah):
            a =  ah

        else:
            a =  ap

        delta_z = delta_z140 + deltaz_SL * (a - 140.0)

        if(ap < ACN/2):
            delta_z *= -1.0

        sigma_z = sigma_z140 + sigz_SL * (a - 140.0)

        fz = fz_z140 + d_fz * (ap - 140.0)
        fn = fn_z140
        fa = wahl_evenodd_term(ZCN, ACN, fz, fn)

        return (delta_z, sigma_z, fa)


    def _symm():
        if(ap < bd["ba"]):
            delta_z = - delta_zmax + sl50 * (ah - bd["bb"])
            sigma_z = sigma_z50

        elif(ap < bd["bb"]):
            ## this condition means   bd["ba"] < ap < bd["bb"]
            delta_z = - delta_zmax + 2.0 * delta_zmax * (bd["bb"] - ah) / (bd["bb"] - bd["ba"])
            sigma_z = sigma_z140 - sigz_SL * (140.0 - bd["bb"])

        else:
            delta_z = delta_zmax - sl50*(ap - bd["bb"])
            sigma_z = sigma_z50

        ## accoring to CFYP.FOR
        # sigma_z = sigma_z50
        # delta_z = delta_zmax - sl50 * (ap - bd["a_max"])

        # if ap <= (ACN- bd["a_max"]):
        #     delta_z = -delta_zmax + sl50 * (ACN - ap - bd["a_max"])

        # if bd["ba"] < ap < bd["bb"]:
        #     delta_z = -delta_zmax + (bd["a_max"] - ACN - ap) * 2.0 * delta_zmax / (bd["bb"] - bd["ba"])
        
        fz = 1.0
        fn = 1.0

        fa = wahl_evenodd_term(ZCN, ACN, fz, fn)

        return (delta_z, sigma_z, fa)

    if(ap < bd["b1"]):
       return  _far_wing()

    elif(ap < bd["b2"]):
        return _wing()

    elif(ap < bd["b3"]):
        return _peak()

    elif(ap < bd["b4"]):
        return _symm()

    elif(ap < bd["b5"]):
        return _peak()

    elif(ap < bd["b6"]):
        return _wing()
        
    else:
        return _far_wing()


def wahl_evenodd_term(ACN, ZCN, fz, fn):
    f = 1.0

    eoz = ZCN % 2
    eon = (ACN - ZCN)%2

    if  ((eoz == 0) and (eon == 0)):
        f = fz * fn
    elif((eoz == 0) and (eon == 1)):
        f = fz / fn
    elif((eoz == 1) and (eon == 0)):
        f = fn / fz
    else:
        f = 1.0/(fz * fn)
    return f


if __name__ == "__main__":
    # print(wahl_systematics(130))
    chain_yield()