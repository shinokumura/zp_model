import pandas as pd
import scipy
import re
import numpy as np
import math

try:
    from .elem import elemtoz
except:
    from elem import elemtoz

def slices(s, *args):
    position = 0
    for length in args:
        yield s[position : position + length]
        position += length



def round_half_up(n, decimals=0):
    multiplier = 10**decimals
    return math.floor(n * multiplier + 0.5) / multiplier


def linear(x_range, a, b):
    return a * x_range + b



def gaussian(x_range, amp, mu, sig):
    if amp == 0:
    # amp = 1/( (2 * math.pi) ** (1/2) * sig)
        amp = 1/(np.sqrt(2*math.pi) * sig)
    return amp * np.exp(-np.power(x_range - mu, 2.0) / (2 * np.power(sig, 2.0)))

    # return amp / (sig*np.sqrt(2)*pi) * np.exp(-np.power(x_range - mu, 2.) / (2 * np.power(sig, 2.)))
    # gauss(x)=amp/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))



def gaussian_math(x_range, amp, mu, sig):
    import math
    if amp == 0:
        amp = 1/( (2 * math.pi) ** (1/2) * sig)
    return [amp * math.exp( (-(x - mu) ** 2) / (2 * sig ** 2) ) for x in x_range]



def erf(x_range, amp, d, sig):
    b = 0.5
    return (
        0.5
        * amp
        * (
            scipy.special.erf((x_range - d + b) / sig * np.sqrt(2))
            - scipy.special.erf((x_range - d - b) / sig * np.sqrt(2))
        )
    )




# ------------------------------------------------------------------------------
# Element
#


def sep_nuclide(nuc):

    if nuc[0].isdigit():
        zcn = nuc[0:3]
        acn = nuc[2:6]
    else:
        zcn = elemtoz(re.sub(r"\d{2,3}", "", nuc))
        acn = re.sub(r"\D{1,2}", "", nuc)

    return int(zcn), int(acn)




def ucd(nuc, a):
    
    ## should check if the nuclide input is compound or fissile
    ## for the safe, call def sep_nuclide(nuc)
    zcn, acn = sep_nuclide(nuc)

    return (int(zcn) / int(acn)) * a



def normalization(dist, normfactor):
    sum = pd.Series(dist).sum(min_count=1)
    if sum != normfactor:
        for i in range(0, len(dist)):
            dist[i] *= normfactor / sum
    return dist



def n_ind_sep_en(mass_ex_1, mass_ex_2):  # neutron induced case of separation energy
    try:
        from .params import ENEUTRON
    except:
        from params import ENEUTRON
    
    mass_ex = mass_ex_1 + float(ENEUTRON) - mass_ex_2
    if mass_ex_1 == 0:
        return 1e10
    else:
        return mass_ex

