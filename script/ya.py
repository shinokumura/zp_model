import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def gen_massdist(fissile=None, GAUSS1=None, GAUSS2=None, GAUSSSYM=None):
    """
    generate gauss type mass distribution with 7 gaussians reading parameter from congig.py
    """
    try:
        from .func import gaussian, normalization, sep_nuclide
    except:
        from func import gaussian, normalization, sep_nuclide

    if GAUSS1 is None:
        try:
            from .params import GAUSS1, GAUSS2, GAUSSSYM
        except:
            from params import GAUSS1, GAUSS2, GAUSSSYM

    if fissile is not None:
        _, ACN = sep_nuclide(fissile)
        x_range = _gen_mass_range(int(ACN / 2 - 50), int(ACN / 2 + 50))

    ###### create two hampt Y(A) from 4 Gaussian functions ######
    massdist = (
        gaussian(x_range, GAUSS1[0], ACN / 2 + GAUSS1[2], GAUSS1[1])
        + gaussian(x_range, GAUSS1[0], ACN / 2 - GAUSS1[2], GAUSS1[1])
        + gaussian(x_range, GAUSS2[0], ACN / 2 + GAUSS2[2], GAUSS2[1])
        + gaussian(x_range, GAUSS2[0], ACN / 2 - GAUSS2[2], GAUSS2[1])
        + gaussian(x_range, GAUSSSYM[0], ACN / 2, GAUSSSYM[1])
    )


    ###### normalize total Y(A) to 2.0  ######
    massdist = normalization(massdist, 2)
    df = pd.DataFrame(columns=["mass", "yield"])
    df["mass"] = x_range
    df["yield"] = massdist


    ##### Show Y(A) plot ######
    plt.plot(x_range, massdist)
    plt.show()

    return df


def _gen_mass_range(MIN_A, MAX_A):
    # from .params import MIN_A, MAX_A

    n_sample = MAX_A - MIN_A + 1

    return np.linspace(MIN_A, MAX_A, n_sample)


def _gen_h_mass_range(ACN, MAX_A):
    # from .params import ACN, MAX_A

    n_sample = MAX_A - int(ACN / 2) + 1

    return np.linspace(int(ACN / 2), MAX_A, n_sample)


def read_massdist(filename):
    """
    read mass distribution from file
    """

    col=["mass", "yield", "dyeild"]
    df = pd.read_csv(filename, comment="#", sep= "\s+", header=None, names=col)
    df = df[df['yield'] != 0]

    return df
