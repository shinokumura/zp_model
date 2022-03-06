import numpy as np
import pandas as pd
import scipy
import re


def linear(x_range, a, b):
    return a * x_range + b


def gaussian(x_range, amp, mu, sig):
    return amp * np.exp(-np.power(x_range - mu, 2.0) / (2 * np.power(sig, 2.0)))
    # return amp / (sig*np.sqrt(2)*pi) * np.exp(-np.power(x_range - mu, 2.) / (2 * np.power(sig, 2.)))
    # gauss(x)=amp/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))


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
ELEMS = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
]


def slices(s, *args):
    position = 0
    for length in args:
        yield s[position : position + length]
        position += length


def ztoelem(z):
    if z == 0:
        elem_name = "n"
    else:
        try:
            z = z - 1
            elem_name = ELEMS[z]
        except ValueError:
            elem_name = ""
    return elem_name


def elemtoz(elem):
    try:
        z = ELEMS.index(elem)
        z = z + 1
        # z = str(z).zfill(3)
    except ValueError:
        z = ""
    return z


def sep_nuclide(nuc):
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


def eshellKTUY():
    eshell_df = pd.read_csv(
        "data/KTUY05_m246.dat",
        sep="\s+",
        index_col=None,
        header=None,
        usecols=[0, 1, 3],
        comment="#",
        names=["ZZ", "NN", "Esh"],
    )

    eshell_df["AA"] = eshell_df["ZZ"] + eshell_df["NN"]
    eshell_df = eshell_df[["ZZ", "AA", "NN", "Esh"]]
    # print(eshell_df)
    return eshell_df


def eshellMoeller():
    eshell_df = pd.read_csv(
        "data/pmshell.dat",
        sep="\s+",
        index_col=None,
        header=None,
        usecols=[0, 1, 2, 3],
        comment="#",
        names=["ZZ", "AA", "NN", "Esh"],
    )
    # print(eshell_df)
    return eshell_df


# def frdmmasstable():
#     df = pd.read_csv(
#         "data/mass-frdm95.dat",
#         sep="\s+",
#         index_col=None,
#         header=None,
#         usecols=[0, 1, 3],
#         comment="#",
#         names=["ZZ", "AA", "fl", "Mexp", "Err", "Mth", "Emic", "beta2", "beta3", "beta4", "beta6",],
#     )
#     # print(eshell_df)
#     return df


def read_mass_table():
    from func import slices
    from params import MASSFILE

    z = []
    a = []
    mselect = []
    with open(MASSFILE) as f:
        lines = f.readlines()
        for line in lines[5:]:
            t_z, t_a, t_s, t_fl, t_mexp, t_dmexp, t_mth = slices(line, 4, 4, 3, 2, 10, 10, 10)
            z.append(int(t_z))
            a.append(int(t_a))

            if not t_mexp.isspace():
                mselect.append(float(t_mexp))
            elif t_mexp.isspace() and not t_mth.isspace():
                mselect.append(float(t_mth))
            else:
                mselect.append(0.0)

    df = pd.DataFrame({"Z": z, "A": a, "Mselect": mselect})

    ## check
    if df.empty:
        raise TypeError()

    return df


mass_table_df = read_mass_table()


def mass_excess(z, a):
    try:
        mx = mass_table_df.loc[(mass_table_df.Z == z) & (mass_table_df.A == a)][
            "Mselect"
        ].to_string(index=False)
    except:
        mx = 0.0
    return float(mx)


def n_ind_sep_en(mass_ex_1, mass_ex_2):  # neutron induced case of separation energy
    from params import ENEUTRON

    mass_ex = mass_ex_1 + float(ENEUTRON) - mass_ex_2
    if mass_ex_1 == 0:
        return 1e10
    else:
        return mass_ex

