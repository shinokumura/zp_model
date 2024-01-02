###############################################################
#
# Katakura systematics based on JAERI-Research 2003-004
#
#                                          April 2022 SO
###############################################################

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .func import gaussian, normalization, mass_excess, n_ind_sep_en, sep_nuclide
from .params import MIN_A, MAX_A

pd.reset_option("display.max_columns")
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("max_colwidth", None)
pd.set_option("display.width", 1200)


'''
Shell energy formula
The S factor is represented as follows by Meyers and Swiatecki: 
Meyers W. D. and Swiatecki W., Nucl. Phys. 81, 1 (1966). eq.(5)
'''
def F(nucleon, N):
    '''
    # N: number of nucleon here
    # M is the nucleon mass: it means M = N for neutrons, M = Z for proton case
    #
    # M_i is the magic number. For the S factor adopting the magic
    # number of 
    # N = 126, 164 and 184 for neutrons and 
    # Z = 82, 100 and 114 for protons. 
    # The numbers N = 164 and Z 100 are introduced to take into account the shell effect of fission fragment.
    '''
    if nucleon == "neutron":
        if N <= 126:
            M_i  = 126
            M_i1 = 82
        elif N <= 164:
            M_i  = 164
            M_i1 = 126
        elif N > 164:
            M_i  = 184 
            M_i1 = 164         
    elif nucleon == "proton":
        if N <= 82:
            M_i  = 82
            M_i1 = 50
        elif N <= 100:
            M_i  = 100
            M_i1 = 82
        elif N > 100:
            M_i  = 114
            M_i1 = 100
    else:
        print("M_i, M_i1 not defined")

    q_i = (3/5) * (M_i ** (5/3) - M_i1 ** (5/3)) / (M_i - M_i1)

    return q_i * (N - M_i1) - 3/5 * (N ** (5/3) - M_i1 ** (5/3) )

# print(F("neutron",126))

def S_factor(N_f, Z_f):
    '''
    ## S is the shell factor defined in Eq. (2.17) but the magic numbers, 
    ## N=164 and Z=I00 introduced by Moriyama and Ohnishi to consider 
    ## the shell effect of fission fragment are removed here.
    '''
    A = Z_f + N_f  # A_f 
    s = ((F("neutron", N_f) + F("proton", Z_f)) / ((1/2) * A) ** (2/3)) - 0.26 * A ** (1/3)
    print(s)
    return 5.8 * s


def gen_massdist_k(fissile=None, N_a = 0, N_s = 0, F = 0, GAUSS1=None, GAUSS2=None, GAUSSSYM=None):
    """
    generate gauss type mass distribution with 7 gaussians reading parameter from congig.py
    """

    if GAUSS1 is None and GAUSS2 is None and GAUSSSYM is None:
        from .params import GAUSS1, GAUSS2, GAUSSSYM

    if fissile is not None:
        _, ACN = sep_nuclide(fissile)
        x_range = _gen_mass_range()

    ###### create two hampt Y(A) from 4 Gaussian functions ######
    ## call gaussian(x_range, amp, mu, sig)
    ## [amp(sig_h1), sig_h1, A_h1 - A_f/2]
    ass_massdist1 = (
          gaussian(x_range, GAUSS1[0], ACN / 2 + GAUSS1[2], GAUSS1[1])
        + gaussian(x_range, GAUSS1[0], ACN / 2 - GAUSS1[2], GAUSS1[1])
    )
    ass_massdist2 = (
          gaussian(x_range, GAUSS2[0], ACN / 2 + GAUSS2[2], GAUSS2[1])
        + gaussian(x_range, GAUSS2[0], ACN / 2 - GAUSS2[2], GAUSS2[1])
    )
    symm_massdist = gaussian(x_range, GAUSSSYM[0], GAUSSSYM[2], GAUSSSYM[1])


    ass_massdist1 = np.array(ass_massdist1) * N_a
    ass_massdist2 = np.array(ass_massdist2) * N_a * F
    symm_massdist = np.array(symm_massdist) * N_s

    massdist = list(ass_massdist1 + ass_massdist2 + symm_massdist)

    ###### normalize total Y(A) to 2.0  ######
    massdist = normalization(massdist, 2)
    print(len(x_range))
    print(len(massdist))
    df = pd.DataFrame(columns=["mass", "yield"])
    df["mass"] = x_range
    df["yield"] = massdist

    print(df)

    ##### Show Y(A) plot ######
    plt.plot(x_range, ass_massdist1, label = "ass1")
    plt.plot(x_range, ass_massdist2, label = "ass2")
    plt.plot(x_range, symm_massdist, label = "symm")
    plt.plot(x_range, massdist, label = "all")
    plt.legend()
    plt.show()

    return df


def _gen_mass_range():
    
    n_sample = MAX_A - MIN_A + 1
    return np.linspace(MIN_A, MAX_A, n_sample)


def amp(sig): 
    return 1 /( (2 * math.pi) ** (1/2) * sig)

# def gen_mass_range_k():
#     from params import MIN_A, MAX_A
#     return [x + 1 for x in range(MIN_A, MAX_A)]



################### 
# B:  parameters from Katakura Systematics
###################
## E* is the excitation energy, that is, the sum of the incident energy (E) and 
## the binding energy (BN) of the incident particle. 


def katakura_systematics(fissile, Ein):
    
    Z_f, A_f = sep_nuclide(fissile)
    N_f = A_f - Z_f

    mass_ex_tn = mass_excess(Z_f, A_f - 1)
    mass_ex_cn = mass_excess(Z_f, A_f)
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)
    print("BN:", sep_cn)

    BN = sep_cn    # 6.545 for U235 as compound
    Ex = BN + Ein
    nu_bar =  1.404 + 0.1067 * (A_f - 236) + (14.986 - 0.1067 * (A_f - 236)) * (1.0 - math.exp(-0.00858 * Ex))
    # nu_bar = 0
    print("nu_bar:", nu_bar)
    R = ( 112.0 + 41.24 * math.sin(3.675 * S_factor(N_f, Z_f)) ) * ( 1.0/( BN ** 0.331 + 0.2067) )/( 1./(Ein ** (0.993) + 0.0951) )
    F = 10.4 - 1.44 * S_factor(N_f, Z_f)
    
    sig_s = 12.6

    ## maybe same as MO systematics
    N_s = 200/(1 + 2 * R)
    N_a = (200 * R)/((1 + F) * (1 + 2 * R))

    sig_h1 = (-25.27 + 0.0345 * A_f + 0.216 * Z_f)* (0.438 + Ein + 0.333 * BN ** 0.333) ** 0.0864 # will retur 2.767 
    sig_h2 = (-30.73 + 0.0394 * A_f + 0.285 * Z_f)* (0.438 + Ein + 0.333 * BN ** 0.333) ** 0.0864 # will return 4.82933
    A_s = (A_f - nu_bar)/2
    A_h1 =  0.5393 * (A_f - nu_bar) + 0.01542 * A_f * (40.2 - (Z_f ** 2)/A_f) ** (1/2)  # will return 133.55
    A_h2 =  0.5612 * (A_f - nu_bar) + 0.01910 * A_f * (40.2 - (Z_f ** 2)/A_f) ** (1/2)  # will return 140.48

    print("S_factor: ", S_factor(N_f, Z_f))
    print("R&F:",  R, F)
    print("N_a, N_s:",  N_a, N_s)
    # print("sig:", sig_h1, sig_h2, sig_s)
    # print("A:  ", A_h1, A_h2, A_s)


    GAUSS1 =   [amp(sig_h1), sig_h1, A_h1 - A_s]
    GAUSS2 =   [amp(sig_h2), sig_h2, A_h2 - A_s]
    GAUSSSYM = [amp(sig_s),  sig_s,  A_s]
    print(GAUSS1)
    print(GAUSS2)
    print(GAUSSSYM)

    # x_range = gen_mass_range_k(70,165)
    ###### create two hampt Y(A) from 4 Gaussian functions ######
    print("Katakura sys")
    ## gen_massdist_k(fissile=None, N_a = 0, N_s = 0, F = 0, GAUSS1=None, GAUSS2=None, GAUSSSYM=None):
    gen_massdist_k(fissile, N_a, N_s, F, GAUSS1, GAUSS2, GAUSSSYM)

    # df = pd.DataFrame(columns=["mass", "yield"])
    # df["mass"] = x_range
    # df["yield"] = massdist
    # plt.plot(x_range, massdist, "ro:", label="data")
    # plt.show()




###################
# A:  parameters from Moriyama-Ohnishi Systematics 
###################
def mo_systematics(fissile, Ein):
    Z_f, A_f = sep_nuclide(fissile)
    N_f = A_f - Z_f

    mass_ex_tn = mass_excess(Z_f, A_f - 1)
    mass_ex_cn = mass_excess(Z_f, A_f)
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)
    print("BN:", sep_cn)

    BN = sep_cn    # 6.545 for U235 as compound
    Ex = BN + Ein

#  alpha and beta represent the pairing effect of the fissioning nuclide (Z_f, A_f))
    ## eq.(2.15) and (2.16)
    if (Z_f % 2 == 0) and (N_f % 2 == 0):
        alpha = +7.55 * A_f ** (-1/2)
    elif (Z_f % 2 != 0) and (N_f % 2 != 0):
        alpha = -7.55 * A_f ** (-1/2)
    elif A_f % 2 != 0:
        alpha = 0
    ## eq.(2.16)
    if A_f % 2 == 0:
        beta = 2.531 - 0.10546 * Ex
    else:
        beta = 2.118 - 0.07990 * Ex

    ## nu_bar systematics Eq.(2.5)
    if Ex <= 7.967:
        nu_bar = 0.107 * A_f - 23.37 + 0.085 * Ex
    elif Ex > 7.967:
        nu_bar = 0.107 * A_f - 23.968 + 0.16 * Ex
    # overwrite nu_bar as 0
    # nu_bar = 0


    ## for light and heavy peaks
    A_h1 =  0.5190 * (A_f - nu_bar) + 0.02840 * A_f * (40.2 - Z_f ** 2 / A_f) ** (1/2)
    A_h2 =  0.5673 * (A_f - nu_bar) + 0.02151 * A_f * (40.2 - Z_f ** 2 / A_f) ** (1/2)
    sig_h1 = 0.07327 * A_f - 14.50 + 0.01783  * Ex ** 0.25 * (S_factor(N_f, Z_f) + 31.908)
    sig_h2 = 0.08089 * A_f - 15.87 + 0.006847 * Ex ** 0.25 * (S_factor(N_f, Z_f) + 31.908)
    A_s = (A_f - nu_bar) /2 

    ## from eq.(2.13) and (2.14)
    ## F the ratio of the asymmet- ric component I to the asymmetric component 2
    R = math.exp( 2.262 * S_factor(N_f, Z_f) + 0.0682 * A_f - 11.542 + alpha + math.exp(beta)) 
    F = sig_h2/sig_h1 * math.exp(-0.08349 * A_f + 19.43)


    ## for the symmetric compornent
    N_s = 200/(1 + 2 * R)
    N_a = 200 * R/((1 + F) * (1 +2 * R))
    sig_s = 0.1759 * A_f - 41.71 + 0.2056 * Ex ** 0.25 * (S_factor(N_f, Z_f) + 31.908)

    print("MO sys")
    print("R&F:",  R, F)
    print("sig:", sig_h1, sig_h2, sig_s)
    print("A:  ", A_h1, A_h2)
    print("N:  ", N_s, N_a)

    GAUSS1 =   [amp(sig_h1), sig_h1, A_h1 - A_s]
    GAUSS2 =   [amp(sig_h2), sig_h2, A_h2 - A_s]
    GAUSSSYM = [amp(sig_s),  sig_s,  A_s]
    print(GAUSS1)
    print(GAUSS2)
    print(GAUSSSYM)

    gen_massdist_k(fissile, N_a, N_s, F, GAUSS1, GAUSS2, GAUSSSYM)
    ###### create two hampt Y(A) from 4 Gaussian functions ######
    ## gaussian(x_range, amp, mu, sig):
    


if __name__ == "__main__":
    from .func import sep_nuclide
    fissile = "U236"
    Ein = 2.53e-8
    katakura_systematics(fissile, Ein)
    mo_systematics(fissile, Ein)

