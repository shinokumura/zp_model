import pandas as pd
import os

from decimal import Decimal, ROUND_HALF_UP
from .func import sep_nuclide

def nuA(fissile, energy):
    file = "".join(["data/nuA-", fissile, "-E", energy, ".dat"])

    if os.path.exists(file):
        nuA_df = pd.read_csv(
            file,
            sep="\s+",
            index_col=None,
            header=None,
            comment="#",
            names=["mass", "nubarA"],
            usecols=[0, 4],
        )

        nuA_df = nuA_df.astype({"mass": "int", "nubarA": "float64",})
        return nuA_df
    else:
        print("no nubar(A) file found for ", fissile, energy)
        return pd.DataFrame()


def nuAsubstraction(fissile, energy, df):
    # for fissile in fissiles:
    # get nu_bar(A) systematics (temporary just read the BeoH results)
    nuA_df = nuA(fissile, energy)

    for index, row in nuA_df.iterrows():
        a = row["mass"]
        nu = row["nubarA"]
        # print (df.loc[(df['Apost'] == a) & (df['fissile'] == fissile)])
        # df.loc[df['Apost'] == a, 'A'] = round(float(a) - nu, 1)
        df.loc[(df["Apost"] == a) & (df["fissile"] == fissile), "A"] = Decimal(
            str(a + nu)
        ).quantize(Decimal("0"), rounding=ROUND_HALF_UP)
        # df['Apre'] = df[df['A'] == a]

    acn, zcn = sep_nuclide(fissile)

    # turn heavy to light
    df.loc[df.A > int(acn) / 2, "A"] = int(acn) - df["A"]
    df.loc[df.Z > int(zcn) / 2, "Z"] = int(zcn) - df["Z"]

    df = df.astype({"Apost": "int", "A": "int",})

    return df
