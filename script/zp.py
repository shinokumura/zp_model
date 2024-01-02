import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


try:
    from .func import gaussian, ucd, normalization, sep_nuclide
    from .mass import eshellKTUY
    from .params import CN, ACN, ZCN, SIGMAZ, Y_CUTOFF
    from .ya import gen_massdist
    from .wahl import wahl_parameter_calc, wahl_systematics
except:
    from func import gaussian, ucd, normalization, sep_nuclide
    from mass import eshellKTUY
    from params import Y_CUTOFF
    from ya import gen_massdist
    from wahl import wahl_parameter_calc, wahl_systematics


def _gen_z_range(MIN_Z, MAX_Z):
    from params import MIN_Z, MAX_Z

    n_sample = MAX_Z - MIN_Z + 1

    return np.linspace(MIN_Z, MAX_Z, n_sample)



def _gen_z_range_in_mass(zucd, deltaZ):

    zucd = zucd + deltaZ

    return np.linspace(int(zucd - 12), int(zucd + 12), 25, dtype = int)




def ucdZp(fissile, mass_df):
    df = pd.DataFrame(columns=["mass", "charge", "yield"])
    for index, row in mass_df.iterrows():
        a = row["mass"]
        y = row["yield"]

        # should be compound nuclide
        zucd = ucd(fissile, a)
        z_range = _gen_z_range_in_mass(zucd)  # charge range in the certain A

        #### def gaussian(x_range, amp, mu, sig):
        chargedist = gaussian(z_range, 1.0, zucd, SIGMAZ)

        chargedist = normalization(chargedist, y)

        df2 = pd.DataFrame(columns=["charge", "mass", "yield"])
        df2["charge"] = z_range
        df2["mass"] = np.linspace(a, a, 11)
        df2["yield"] = chargedist

        df = df.append(df2, ignore_index=True)

    return df



def readZpFromFile():
    file = "/Users/okumuras/Dropbox/calculations/cp_ebata/FP_E1800_SkMs_fit.dat"
    # file = "/Users/okumuras/Dropbox/calculations/cp_ebata/FP_E1793_SkMs_fit.dat"
    
    dict = {}

    with open(file) as f:
        lines = f.readlines()
        # print(lines.split())

    for line in lines:
        if line.startswith("#"):
            continue

        t = line.split()
        dict[t[0]] = t[1]
       
    return dict





def fractionZp(params_dict, mass_df, model = "wahl_zp", eo_model = ""):
    df = pd.DataFrame(columns=["charge_h", "mass_h", "pair_yield"])

    CN = params_dict["CN"]
    ACN = params_dict["ACN"]

    # Modified temporary for the CP from the mean field theory # 2023-02-20
    # deltaZ_A = readZpFromFile()
    
    if model == "wahl_zp":
        wahl_dict = wahl_parameter_calc(params_dict)

    # print("{:10}{:10}{:10}{:10}".format(CN, CN, CN, CN))
    for _, row in mass_df.iterrows():
        a = row["mass"]
        y = row["yield"]

        ## because it is symmetric
        if a >= ACN / 2:

            if model == "wahl_zp":
                ## full Wahl systematics
                if ACN < 260:
                
                    wahl_params = wahl_systematics(params_dict, wahl_dict, a)
                    
                    deltaZ = wahl_params[0] if wahl_params[0] > -1 else -1.0
                    sigZ   = wahl_params[1] if wahl_params[1] < 1 else 1.0
                    fa     = wahl_params[2] if wahl_params[2] < 1.8 else 1.8

                else:
                    deltaZ = 0.0
                    sigZ   = 0.5
                    fa = 1.0

                # print("{:5.0F}{:10.3F}{:10.3F}{:10.3F}".format(a, deltaZ, sigZ, fa))

            elif model == "fixed":
                ## partially Wahl systematics + changeable even-odd model
                from .params import DELTAZ, SIGMAZ

                deltaZ = DELTAZ
                sigZ = SIGMAZ

                if eo_model == "wahl":
                    fa = _wahl_evenodd_term(z, n)

                elif eo_model == "minato":
                    fa = _minato_evenodd_term(z, n)
                    
                elif eo_model == "titech":
                    fa = _titech_evenodd_term(z, n)

            else:
                deltaZ = _wahl_deltaZ_term()
                sigZ = _wahl_sigma_term
                fa = 1.0
            
            ## determin the Z range wothin this A
            zucd = ucd(CN, a)
            z_range_in_a = _gen_z_range_in_mass(zucd, deltaZ)
            zp = ucd(CN, a) + deltaZ
            
            # print(z_range_in_a, a, deltaZ, sigZ, fa, "--> ", zucd, zp)

            fracyields = []
            for z in z_range_in_a:
                n = a - z
                v = (z - zp + 0.5) / (sigZ * math.sqrt(2))
                w = (z - zp - 0.5) / (sigZ * math.sqrt(2))

                ## each mass's fractional yield
                fracyield = 0.5 * fa * (math.erf(v) - math.erf(w))

                ## list of yields in the same mass
                if fracyield >= Y_CUTOFF:
                    fracyields.append(fracyield)

                else:
                    fracyields.append(0.0)

                ## normalize the fractional yield with the mass yield
            fracyields = normalization(fracyields, y)

            ## temporary dataframe
            df2 = pd.DataFrame(columns=["charge_h", "mass_h", "pair_yield"])
            df2["charge_h"] = z_range_in_a
            df2["mass_h"] = np.linspace(a, a, len(z_range_in_a))
            df2["pair_yield"] = fracyields
            df = pd.concat([df, df2])


    ## convert data type to float
    df = df[df['pair_yield'] != 0]
    df["pair_yield"] = df["pair_yield"].astype("float64")

    ## store pair yield into dataframe
    df.insert(0, "k", range(1, len(df) + 1))

    return df




def _minato_evenodd_term(z, n):
    xi = 0.830
    en = 2.53e-08  # in MeV??
    endependterm = 2 / (1 + np.exp(xi * en))

    if (z % 2) == 0 and (n % 2) == 0:
        c = 0.344
        return 1 + (c * endependterm)

    elif (z % 2) == 0 and (n % 2) != 0:
        c = 0.202
        return 1 + (c * endependterm)

    elif (z % 2) != 0 and (n % 2) == 0:
        c = 0.202
        return 1 - (c * endependterm)

    elif (z % 2) != 0 and (n % 2) != 0:
        c = 0.344
        return 1 - (c * endependterm)

    else:
        return 1.0




def _wahl_evenodd_term(z, n):
    fz = 1 + (1.207 - 1.0) * 1.8
    fn = 1 + (1.076 - 1.0) * 0.9
    if (z % 2) == 0 and (n % 2) == 0:
        fa = fz * fn
    elif (z % 2) == 0 and (n % 2) != 0:
        fa = fz / fn
    elif (z % 2) != 0 and (n % 2) == 0:
        fa = fn / fz
    elif (z % 2) != 0 and (n % 2) != 0:
        fa = 1 / (fz * fn)
    else:
        fa = 1.0

    return fa





def _wahl_deltaZ_term():
    deltaZ = -0.487
    return deltaZ




def _wahl_sigma_term():
    sig = 0.566
    return sig



def _titech_evenodd_term(z, n):
    """
    https://doi.org/10.1080/00223131.2020.1813643
    """
    eshell_df = eshellKTUY()
    esh = eshell_df[(eshell_df["ZZ"] == z) & (eshell_df["NN"] == n)]["Esh"]

    if (z % 2) == 0 and (n % 2) == 0:
        epair = -12 / math.sqrt(z + n)
    elif (z % 2) != 0 and (n % 2) != 0:
        epair = 12 / math.sqrt(z + n)
    else:
        epair = 0.0

    ed = 1.0 / 0.37  # Taken from FPY_param Figure 14: beta = 1/Ed

    foe = math.exp(-(esh + epair) / ed)

    return foe



if __name__ == "__main__":
    fractionZp("wahl")


