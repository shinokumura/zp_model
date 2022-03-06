#######################################################
### 
###              Wahl's Zp model  
### 
#######################################################

from params import ACN, ZCN, EIN
from func import mass_excess, n_ind_sep_en
import math

def zp_parameter_function(p):
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
        p1 = p[0] + p[2]*dz
        p2 = p[1] + p[3]*dz
        x = p1 + (p2 - p1)*(1.0 - math.exp(-p[4]*de))

    return x

sigma_z140 = zp_parameter_function([  0.566 ,  0.0,  0.0064,  0.0109,  0.0    ])
delta_z140 = zp_parameter_function([ -0.487 ,  0.0,  0.0180,  0.0   , -0.00203])
fz_z140    = zp_parameter_function([  1.207 ,  0.0, -0.0420,  0.0   ,  0.0022 ])
fn_z140    = zp_parameter_function([  1.076 ,  0.0,  0.0   ,  0.0   ,  0.0    ])
d_sigmaz   = zp_parameter_function([ -0.0038,  0.0,  0.0   ,  0.0   ,  0.0    ])
d_deltaz   = zp_parameter_function([ -0.0080,  0.0,  0.0   ,  0.0   ,  0.0    ])
d_fz       = zp_parameter_function([  0.0030,  0.0,  0.0   ,  0.0   ,  0.0    ])

d50        = zp_parameter_function([  0.191,  0.0  , -0.0076,  0.0,  0.0])
sigma_z50  = zp_parameter_function([  0.356,  0.060,  0.0   ,  0.0,  0.0])
delta_zmax = zp_parameter_function([  0.699,  0.0  ,  0.0   ,  0.0,  0.0])

if d50 < 0.0:
    d50 = 0.0

sigma_zwing = zp_parameter_function([   -0.045,  0.0094,  0.0,  0.0,  0.0])
delta_zwing = zp_parameter_function([    0.0  , -0.0045,  0.0,  0.0,  0.0])
fz_zwing    = zp_parameter_function([    0.159, -0.028 ,  0.0,  0.0,  0.0])
fn_zwing    = zp_parameter_function([    0.039,  0.0   ,  0.0,  0.0,  0.0])

def mass_region_boundary():
    bundary_dict = {}

    def _amax(delta_zmax, d50):
        f1 = (250.0 - ACN)/(250.0 - 236.0)

        if(f1 < 0.0):
            f1 = 0.0
        elif(f1 > 1.0):
            f1 = 1.0

        f2 = 1.0 - f1
        r  = ACN / ZCN
        ak1 = 50.0*r - delta_zmax/d50
        ak2 = (50.0 - delta_zmax)*r

        return f1 * ak1 + f2 * ak2

    a_max = _amax(delta_zmax, d50)
    delta_max = 2.0
    b1 = 70
    b2 = 77 + 0.036 * (ACN - 236)
    b4 = (delta_zmax - delta_z140 + a_max * d50 + 140 * d_deltaz)/(d50 + d_deltaz)
    b3 = ACN - b4
    b5 = ACN - b2
    b6 = ACN - b1
    ba = a_max
    bb = ACN - a_max

    bundary_dict = { "b1": b1, "b2": b2, "b3": b3, "b4": b4, "b5": b5, "b6": b6, "ba": ba, "bb": bb }

    return bundary_dict

def wahl_systematics(a):

    bd = mass_region_boundary()
    ap = a
    ah = ACN - ap

    def _far_wing():
        sigma_z = sigma_z140 + d_sigmaz * (bd["b5"] - 140.0)
        delta_z = delta_z140 + d_deltaz * (bd["b5"] - 140.0) 
        if(ap < ah):
            delta_z *= -1.0
        fz     = fz_z140
        fn     = fn_z140
        fa = wahl_evenodd_term(fz, fn)

        return (delta_z, sigma_z, fa)

    def _wing():
        deltaB5 = delta_z140 + d_deltaz * (bd["b5"] - 140.0) 
        sigmaB5 = sigma_z140 + d_sigmaz * (bd["b5"] - 140.0) 

        if(ap < ah):
            delta_z = - deltaB5 - delta_zwing * (bd["b5"] - ah)
            sigma_z =   sigmaB5 + sigma_zwing * (bd["b2"] - ap)
            fz = fz_z140 + fz_zwing * (bd["b2"] - ap)
            fn = fn_z140 + fn_zwing * (bd["b2"] - ap)
        
        else:
            delta_z = deltaB5   - delta_zwing * (ap - bd["b5"])
            sigma_z = sigmaB5   + sigma_zwing * (ap - bd["b5"])
            fz = fz_z140 + fz_zwing * (ap - bd["b5"])
            fn = fn_z140 + fn_zwing * (ap - bd["b5"])
        fa = wahl_evenodd_term(fz, fn)

        return (delta_z, sigma_z, fa)
    

    def _peak():
        if (ap < ah):
            a =  ah 
        else:
            a =  ap

        delta_z = delta_z140 + d_deltaz * (a - 140.0)
        if(ap < ah):
            delta_z *= -1.0
        sigma_z = sigma_z140 + d_sigmaz * (a - 140.0)

        fz = fz_z140 + d_fz * (a - 140.0)
        fn = fn_z140
        fa = wahl_evenodd_term(fz, fn)

        return (delta_z, sigma_z, fa)


    def _symm():
        if(ap < bd["ba"]):
            delta_z = - delta_zmax + d50*(ah - bd["bb"])
            sigma_z = sigma_z50
        
        elif(ap < bd["bb"]):
            delta_z = - delta_zmax +  2.0*delta_zmax * (bd["bb"] - ah) / (bd["bb"] - bd["ba"])
            sigma_z = sigma_z140 - d_sigmaz * (140.0 - bd["bb"])
        
        else:
            delta_z = delta_zmax - d50*(ap - bd["bb"])
            sigma_z = sigma_z50
        
        fz = 1.0
        fn = 1.0

        fa = wahl_evenodd_term(fz, fn)
        
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


def wahl_evenodd_term(fz,fn):
    f = 1.0

    eoz = ZCN % 2
    eon = (ACN - ZCN)%2

    if     ((eoz == 0) and (eon == 0)):
        f = fz * fn
    elif((eoz == 0) and (eon == 1)):
        f = fz / fn
    elif((eoz == 1) and (eon == 0)):
        f = fn / fz
    else:
        f = 1.0/(fz * fn)
    return f


if __name__ == "__main__":
    print(wahl_systematics(130))