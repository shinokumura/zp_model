import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from func import gaussian, ucd, normalization, eshellKTUY, sep_nuclide
from params import ACN, ZCN, SIGMAZ, Y_CUTOFF
from ya import gen_massdist


def _gen_z_range():
    from params import MIN_Z, MAX_Z

    n_sample = MAX_Z - MIN_Z + 1

    return np.linspace(MIN_Z, MAX_Z, n_sample)


def _gen_z_range_in_mass(cn, zucd):
    zcn, _ = sep_nuclide(cn)

    if zucd < zcn / 2:
        zucd = zucd + 0.5
    elif zucd > zcn / 2:
        zucd = zucd - 0.5

    return np.linspace(int(zucd - 6), int(zucd + 6), 13)


def ucdZp(fissile, mass_df):
    df = pd.DataFrame(columns=["mass", "charge", "yield"])
    for index, row in mass_df.iterrows():
        a = row["mass"]
        y = row["yield"]

        # print(a, zucd)
        # should be compound nuclide
        zucd = ucd(fissile, a)
        z_range = _gen_z_range_in_mass(zucd)  # charge range in the certain A
        # print(x_range)
        #### def gaussian(x_range, amp, mu, sig):
        chargedist = gaussian(z_range, 1.0, zucd, SIGMAZ)
        # print(chargedist)

        chargedist = normalization(chargedist, y)

        df2 = pd.DataFrame(columns=["charge", "mass", "yield"])
        df2["charge"] = z_range
        df2["mass"] = np.linspace(a, a, 11)
        df2["yield"] = chargedist

        df = df.append(df2, ignore_index=True)

    return df


def fractionZp(cn, mass_df, eo_model="wahl"):
    df = pd.DataFrame(columns=["charge_h", "mass_h", "pair_yield"])

    for _, row in mass_df.iterrows():
        # print (row)
        a = row["mass"]
        y = row["yield"]

        # because it is symmetric
        if a >= ACN / 2:
            # determin the Z range in particular A
            zucd = ucd(cn, a)
            z_range_in_a = _gen_z_range_in_mass(cn, zucd)
            fracyields = []

            for z in z_range_in_a:
                n = a - z

                deltaZ = _wahl_deltaZ_term()
                sigZ = _wahl_sigma_term()

                if a > ACN / 2:  # heavy fragment
                    zp = ucd(cn, a) + deltaZ
                elif a <= ACN / 2:  # light fragment
                    zp = ucd(cn, a) - deltaZ

                if eo_model == "wahl":
                    fa = _wahl_evenodd_term(z, n)
                elif eo_model == "minato":
                    fa = _minato_evenodd_term(z, n)
                elif eo_model == "titech":
                    fa = _titech_evenodd_term(z, n)

                # print (a, z, deltaZ, sigZ, zp, fa)

                v = (z - zp + 0.5) / (sigZ * math.sqrt(2))
                w = (z - zp - 0.5) / (sigZ * math.sqrt(2))

                ## each mass's fractional yield
                fracyield = 0.5 * fa * (math.erf(v) - math.erf(w))
                # print(a, z, fracyield)

                ## list of yields in the same mass
                if fracyield >= Y_CUTOFF:
                    fracyields.append(fracyield)
                else:
                    fracyields.append(0.0)

            ## normalize the fractional yield with the mass yield
            fracyields = normalization(fracyields, y)
            # print(a, fracyields)

            ## temporary dataframe
            df2 = pd.DataFrame(columns=["charge_h", "mass_h", "pair_yield"])
            df2["charge_h"] = z_range_in_a
            df2["mass_h"] = np.linspace(a, a, len(z_range_in_a))
            df2["pair_yield"] = fracyields
            df = pd.concat([df, df2])

    ## just convert data type
    df["pair_yield"] = df["pair_yield"].astype("float64")

    ## pair dataframe version
    df.insert(0, "k", range(1, len(df) + 1))
    ## check if it's ok
    # print(df)
    ## sum of yield should be == 2.0, or 1.0 if only heave
    # print(df.sum())

    return df


def _minato_evenodd_term(z, n):
    xi = 0.830
    en = 2.53e-08  # in MeV??
    endependterm = 2 / (1 + np.exp(xi * en))

    if (z % 2) == 0 and (n % 2) == 0:
        c = 0.344
        return 1 + (c * endependterm)
        # return 1 + c

    elif (z % 2) == 0 and (n % 2) != 0:
        c = 0.202
        return 1 + (c * endependterm)
        # return 1 + c

    elif (z % 2) != 0 and (n % 2) == 0:
        c = 0.202
        return 1 - (c * endependterm)
        # return 1 - c

    elif (z % 2) != 0 and (n % 2) != 0:
        c = 0.344
        return 1 - (c * endependterm)
        # return 1 - c

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

    # print(z, n, fa)
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

    # the dumping energy that defines the extent for washing out of the shell and pairing effects. The parameter Ed(A) was adjusted for each mass number A of fission product to reproduce the odd-even staggering or distortion from the Gaussian distribution.
    ed = 1.0 / 0.37  # Taken from FPY_param Figure 14: beta = 1/Ed

    foe = math.exp(-(esh + epair) / ed)

    return foe


if __name__ == "__main__":
    fractionZp("wahl")


"""
def WahlZp():
    df = pd.DataFrame(columns = ['mass', 'charge', 'yieldWahl'])
    
    for index, row in mass_df.iterrows():
        a = row['mass']
        y = row['yield']
        
        zucd = ucd(a)
        z_range = get_z_range(zucd)

        chargedist = []
        for z in z_range:
            n = a - z

            deltaZ = WahlDeltaZTerm()
            sigZ = WahlSigmaTerm()

            if a > ACN/2:     # heavy fragment
                zp = ucd(a) + deltaZ
            elif a <= ACN/2:  # light fragment
                zp = ucd(a) - deltaZ

            fa = WahlEvenOddTerm(z, n)

            # print (a, z, deltaZ, sigZ, zp, fa)

            v = (z - zp + 0.5) / (sigZ * math.sqrt(2)) 
            w = (z - zp - 0.5) / (sigZ * math.sqrt(2)) 

            fracyield = 0.5 * fa * (math.erf(v) - math.erf(w))

            chargedist.append(fracyield)
        
        chargedist = normalization(chargedist, y)
        # print(chargedist)
        df2 = pd.DataFrame(columns = ['mass', 'charge', 'yieldWahl'])
        df2['mass'] = np.linspace(a, a, 11)
        df2['charge'] = z_range
        df2['yieldWahl'] = chargedist

        df = df.append(df2, ignore_index=True)
    # print(df)

    return df

def MinatoZp():
    df = pd.DataFrame(columns = ['mass', 'charge', 'yieldMinato'])
    for index, row in mass_df.iterrows():
        a = row['mass']
        y = row['yield']
        
        zucd = ucd(a)
        z_range = get_z_range(zucd)

        chargedist = []
        for z in z_range:
            n = a - z

            deltaZ = WahlDeltaZTerm()
            sigZ = WahlSigmaTerm()

            if a > ACN/2:     # heavy fragment
                zp = ucd(a) + deltaZ
            elif a <= ACN/2:  # light fragment
                zp = ucd(a) - deltaZ

            fa = MinatoEvenOddTerm(z, n)

            # print (a, z, deltaZ, sigZ, zp, fa)

            v = (z - zp + 0.5) / (sigZ * math.sqrt(2)) 
            w = (z - zp - 0.5) / (sigZ * math.sqrt(2)) 

            fracyield = 0.5 * fa * (math.erf(v) - math.erf(w))

            chargedist.append(fracyield)
        
        chargedist = normalization(chargedist, y)
        # print(chargedist)
        df2 = pd.DataFrame(columns = ['mass', 'charge', 'yieldMinato'])
        df2['mass'] = np.linspace(a, a, 11)
        df2['charge'] = z_range
        df2['yieldMinato'] = chargedist

        df = df.append(df2, ignore_index=True)
    print (df)

    return df

def TitechZp():
    df = pd.DataFrame(columns = ['mass', 'charge', 'yieldTitech'])
    
    for index, row in mass_df.iterrows():
        a = row['mass']
        y = row['yield']
        
        zucd = ucd(a)
        z_range = np.linspace(int(zucd - 5), int(zucd + 5), 11) 

        chargedist = []
        for z in z_range:
            n = a - z

            deltaZ = WahlDeltaZTerm()
            sigZ = WahlSigmaTerm()

            if a > ACN/2:     # heavy fragment
                zp = ucd(a) + deltaZ
            elif a <= ACN/2:  # light fragment
                zp = ucd(a) - deltaZ

            fa = TitechEvenOddTerm(z, n)

            # print (a, z, deltaZ, sigZ, zp, fa)

            v = (z - zp + 0.5) / (sigZ * math.sqrt(2)) 
            w = (z - zp - 0.5) / (sigZ * math.sqrt(2)) 

            fracyield = 0.5 * fa * (math.erf(v) - math.erf(w))

            chargedist.append(fracyield)
        
        chargedist = normalization(chargedist, y)
        # print(chargedist)
        df2 = pd.DataFrame(columns = ['mass', 'charge', 'yieldTitech'])
        df2['mass'] = np.linspace(a, a, 11)
        df2['charge'] = z_range
        df2['yieldTitech'] = chargedist

        df = df.append(df2, ignore_index=True)
    # print(df)

    return df
"""
