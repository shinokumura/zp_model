import pandas as pd
import numpy as np
import os


def nubar_systematics():
    import matplotlib.pyplot as plt

    fissiles = ["U233", "Np237"]
    # fissiles = ['Cf252']
    energy = "E2.53E-08"
    # energy = 'SF'
    for fissile in fissiles:
        file = "".join(["data/exp/nuA/exp_nuA-", fissile, "-", energy, ".dat"])
        # print(file)

        if os.path.exists(file):
            if energy == "SF":
                df = pd.read_csv(
                    file,
                    sep="\s+",
                    index_col=None,
                    header=None,
                    comment="#",
                    names=["nubarA", "dnubarA", "A"],
                    usecols=[0, 1, 2],
                )
            else:
                df = pd.read_csv(
                    file,
                    sep="\s+",
                    index_col=None,
                    header=None,
                    comment="#",
                    names=["nubarA", "dnubarA", "A"],
                    usecols=[0, 1, 4],
                )
        df = df[df["A"] >= 50]

        df = df.astype({"A": "int", "nubarA": "float64", "dnubarA": "float64",})

        # plt.title(fissile)
        # plt.plot(df['A'],df['nubarA'],'ro:',label='data')
        # plt.show()

        x_data = range(df["A"].min(), df["A"].max() + 1)
        y_data = []
        y = []
        dy = []
        for a in x_data:
            y.append(df[df["A"] == a]["nubarA"].mean())
            dy.append(df[df["A"] == a]["dnubarA"].mean())

        y_data = np.array(y)
        y_error = np.array(dy)

        plt.title(fissile)
        plt.plot(x_data, y_data, "ro:", label="data")
        plt.show()


nubar_systematics()
