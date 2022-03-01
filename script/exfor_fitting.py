# read exfor data
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import os
import re
import glob

from func import (
    gaussian,
    erf,
    ucd,
    get_z_range,
    eshellKTUY,
    eshellMoeller,
    normalization,
)
from params import SIGMAZ, DELTAZ, ACN, ZCN
from nubar import nuAsubstraction

# https://teratail.com/questions/305680
# https://base64.work/so/python/780806

pd.reset_option("display.max_columns")
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("max_colwidth", None)
pd.set_option("display.width", 1200)


# EXFOR_DATA_PATH = '/Users/sin/Documents/nucleardata/exforpyplot2/libraries/FY/n/'
EXFOR_DATA_PATH = "/Users/sin/Documents/nucleardata/exforpyplot/libraries/FY/n/"


def create_exfy(exfiles, fissile, energy):
    dfs = []
    for e in exfiles:
        # datasetname = os.path.basename(e)
        datasetname = re.split("[-]", os.path.basename(e), 5)
        df = pd.read_csv(
            e,
            sep="\s+",
            index_col=None,
            header=None,
            comment="#",
            names=["Z", "Apost", "Iso", "FPY", "dFPY"],
        )

        df = df[df["Z"] != "0****"]  # to drop error file from C5
        df = df[df["Z"] != 0]
        df = df[df["Apost"] != 0]
        df = df[df["Iso"] == -1]  # data should be without isomeric states
        df = df[df["Z"] < df["Apost"]]
        df["author"] = datasetname[3]
        df["entry"] = datasetname[4]
        df["year"] = re.split("[.]", datasetname[-1])[-1]
        df["A"] = df["Apost"]
        df["fissile"] = fissile
        df["energy"] = energy
        dfs.append(df)

    if dfs:
        df = pd.concat(dfs, ignore_index=True)
        df = df.astype(
            {
                "Z": "int",
                "Apost": "int",
                "Iso": "int",
                "FPY": "float64",
                "dFPY": "float64",
                "A": "int",
            }
        )  # , errors='ignore'

    # substract Nu_bar(A) for each fissile system
    df = nuAsubstraction(fissile, energy, df)

    # sort by mass, otherwise plot would be scattered
    df = df.sort_values(["A", "Z"])
    df = df[["Z", "A", "Apost", "Iso", "FPY", "dFPY", "author", "entry"]]

    print(df)
    return df


# fitting to the normal distribution
def gauss_fit(zucd, z_data, y_data, y_error):
    # p0: the initial guess
    # mean = sum(z_data)/len(z_data)

    # popt: parameter optimized, pcov: covariance
    # should be done with more than 3 data points
    try:
        # to get more presise fitting, the initial guess of zucd should be zucd+/-deltaZ
        popt, pcov = curve_fit(gaussian, z_data, y_data, p0=[1.0, zucd, SIGMAZ], maxfev=1000)
        # print (popt)
    except:
        popt = [0.0, 0.0, SIGMAZ]
        pcov = []

    return popt, pcov


# fitting to the error function
def erf_fit(zucd, z_data, y_data, y_error):
    ######  popt,pcov = erf_fit(zucd, z_data, y_data, y_error)  ######
    # zucd will be given as an initial parameter of Z(A) - Zp(A) where Zp(A) = Zucd - deltaZ

    try:
        # popt, pcov = curve_fit(erf, z_data, y_data, p0=[1.0, zucd, SIGMAZ], maxfev=1000)
        popt, pcov = curve_fit(
            erf,
            z_data,
            y_data,
            p0=[0.01, zucd + DELTAZ, SIGMAZ],
            sigma=1 / y_error ** 2,
            absolute_sigma=True,
            maxfev=1000,
        )

    except:
        popt = [0.0, 0.0, SIGMAZ]
        pcov = []

    return popt, pcov


def isobaric_plot(fissile, a, z_data, y_data, z_range, popt, popt2):
    # original data
    plt.plot(z_data, y_data, "ro:", label="data")

    # get ucdwise Z
    zucd = ucd(fissile, a) + DELTAZ
    ## create 11 points of lineshape of z range
    z_range = get_z_range(zucd)

    # complete gaussian shape
    if len(popt) > 1:
        # fitted data point using z point
        plt.plot(z_data, gaussian(z_data, *popt), "cd")
        plt.plot(z_range, gaussian(z_range, *popt), "b-.", label="fit")

    if len(popt2) > 1:
        plt.plot(z_data, erf(z_data, *popt2), "cd")
        plt.plot(z_range, erf(z_range, *popt2), "g-.", label="fit")

    plt.title("".join([fissile, "(n_th,f) at A=", str(a)]))
    plt.yscale("log")

    return plt.show()


def get_value(func, a, z_data, popt, y_data, y_error):
    # gaussian(x_range, amp, mu, sig)

    # to get fitted yield at z in z_data (to compare with experimental average):
    if func == "gauss":
        value = gaussian(z_data, popt[0], popt[1], popt[2])

    elif func == "erf":
        value = erf(z_data, popt[0], popt[1], popt[2])
        # value2 = normalization(value, 1.0)
    else:
        value = gaussian(z_data, popt[0], popt[1], popt[2])

    # print (value)

    df = pd.DataFrame(
        data={"A": a, "Z": z_data, "value": value, "y_data": y_data, "y_error": y_error,},
        columns=["A", "Z", "value", "y_data", "y_error"],
    )

    # difference between experimental average and fitted data
    df["diff"] = df["y_data"] - df["value"]

    return df


def merge_value_esh(df):
    eshell_df = eshellKTUY()
    # eshell_df =  eshellMoeller()
    # print (df.dtypes)

    ## C = np.where(cond, A, B)
    ## defines C to be equal to A where cond is True, and B where cond is False.
    # df['eshKTUY'] = np.where((df['a'] == eshell_df['AA']) & (df['z'] == eshell_df['ZZ']), eshell_df['Esh'], np.NaN)

    # looping over the A, Z in the experimental dataframe and take eshell_df to merge, comparison with mass of PRE
    for index, row in df.iterrows():
        a = row["A"]
        z = row["Z"]
        try:
            df.loc[(df["A"] == a) & (df["Z"] == z), "Esh"] = eshell_df.loc[
                (eshell_df["AA"] == a) & (eshell_df["ZZ"] == z), "Esh"
            ].item()
        except:
            df.loc[(df["A"] == a) & (df["Z"] == z), "Esh"] = 0.0
    return df


def store_params(fissile, a, popt):
    aa = []
    ucdz = []
    amp = []
    mu = []
    sig = []
    deltaZ = []

    if not popt[0] == 0:
        aa += [a]
        amp += [popt[0]]
        ucdz += [ucd(fissile, a)]
        mu += [popt[1]]
        sig += [popt[2]]
        deltaZ += [popt[1] - ucd(fissile, a)]

    return pd.DataFrame(
        data={"a": aa, "amp": amp, "mu": mu, "ucdz": ucdz, "sig": sig, "deltaZ": deltaZ,},
        columns=["a", "amp", "mu", "ucdz", "sig", "deltaZ"],
    )


def main(func):
    fissiles = ["U235", "Pu239"]
    # fissiles = ['Pu239']
    energy = "2.53E-08"
    exfiles = []

    for fissile in fissiles:
        print("# ", fissile)
        indfpy_files = "".join([EXFOR_DATA_PATH, fissile, "/exfor/454/*", energy, "*"])
        exfiles = glob.glob(indfpy_files)

        # read files
        # read exp. data nad subtract nu_bar(A)
        df = create_exfy(exfiles, fissile, energy)

        # print(df)
        value_df = pd.DataFrame()
        param_df = pd.DataFrame()

        for a in range(df["A"].min(), df["A"].max() + 1):
            if (
                not df[df["A"] == a]["Z"].empty and df[df["A"] == a]["Z"].shape[0] >= 3
            ):  # if data point exist more than 3
                y_data = []
                z_data = []
                y = []
                dy = []

                ## get unique Z values
                z_data = df[df["A"] == a]["Z"].unique()
                ## loop over z to get experimental average
                for z in z_data:
                    y.append(df[(df["A"] == a) & (df["Z"] == z)]["FPY"].mean())
                    dy.append(df[(df["A"] == a) & (df["Z"] == z)]["dFPY"].mean())

                ## convert to numpy array
                y_data = np.array(y)
                y_error = np.array(dy)

                ## get ucdwise Z
                zucd = ucd(fissile, a)  # + DELTAZ
                # generate wide range of z for plot
                z_range = get_z_range(fissile, zucd)

                if func == "gauss":
                    ## fitting to Gauss function on Z distribution
                    popt, pcov = gauss_fit(zucd, z_data, y_data, y_error)

                    ## check by plot
                    # isobaric_plot(fissile, a, z_data, y_data, z_range, popt, '')
                elif func == "erf":
                    ## fitting to Error function on Z distribution
                    popt, pcov = erf_fit(zucd, z_data, y_data, y_error)

                    ## check by plot
                    isobaric_plot(fissile, a, z_data, y_data, z_range, "", popt)
                else:
                    popt, pcov = gauss_fit(zucd, z_data, y_data, y_error)
                    popt2, pcov2 = erf_fit(zucd, z_data, y_data, y_error)
                    ## check by plot
                    # isobaric_plot(fissile, a, z_data, y_data, z_range, popt, popt2)

                # Store parameters for each A
                temp2 = store_params(fissile, a, popt)
                param_df = param_df.append(temp2, ignore_index=True)

                ## get value at fitted gaussian function
                temp = get_value(func, a, z_data, popt, y_data, y_error)
                value_df = value_df.append(temp, ignore_index=True)

        value_df = merge_value_esh(value_df)
        df = merge_value_esh(df)

        outfile = True
        if outfile:

            name = "".join(["output/", fissile, "-", energy, "_exp.dat"])
            with open(name, "w") as d_file:
                d_file.write(df.to_string(index=False))

            name = "".join(["output/", fissile, "-", energy, "_param.dat"])
            with open(name, "w") as p_file:
                p_file.write(param_df.to_string(index=False))

            name = "".join(["output/", fissile, "-", energy, "_average.dat"])
            with open(name, "w") as a_file:
                a_file.write(value_df.to_string(index=False))

        # print(value_df.to_string(index=False))

    return param_df, value_df


def read_yza(fissile, energy):

    file = "".join(["data/Yza-", fissile, "-E", energy, ".dat"])

    if os.path.exists(file):
        df = pd.read_csv(
            file,
            sep="\s+",
            index_col=None,
            header=None,
            usecols=[2, 3, 4, 5, 6],
            comment="#",
            names=["Zl", "Al", "Zh", "Ah", "preyield"],
        )

    df["Nl"] = df["Al"] - df["Zl"]
    df["Nh"] = df["Ah"] - df["Zh"]

    param_df = pd.DataFrame()

    for a in range(df["Al"].min(), df["Al"].max() + 1):
        zucd = ucd(fissile, a)
        z_range = get_z_range(fissile, zucd)

        popt, pcov = erf_fit(
            zucd, df.loc[df["Al"] == a, "Zl"], df.loc[df["Al"] == a, "preyield"], 0.1
        )

        print(popt, pcov)

        temp2 = store_params(fissile, a, popt)
        # print(temp2)
        param_df = param_df.append(temp2, ignore_index=True)

    return df, param_df


if __name__ == "__main__":

    func = "erf"  # or 'erf' 'both'

    param_df, value_df = main(func)
    # df, param_df = read_yza('U235', '2.53E-08')
    # print(param_df)
    # print(param_df.to_string(index=False))
    # print(value_df.to_string(index=False))
