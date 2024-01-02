import pandas as pd

from .zp import (
    ucdZp,
    fractionZp,
    WahlZp,
    MinatoZp,
    TitechZp,
    WahlEvenOddTerm,
    MinatoEvenOddTerm,
    TitechEvenOddTerm,
)


def main():

    options = ["wahl", "minato", "titech"]
    for switch in options:
        df = fractionZp(switch)
        df = pd.merge(df)

    # df_wahl = WahlZp()
    # df_gauss = ucdZp()
    # df_minato = MinatoZp()
    # df_titech = TitechZp()

    # df = pd.merge(df_wahl, df_gauss)
    # df = pd.merge(df, df_minato)
    # df = pd.merge(df, df_titech)

    # print(df)

    ###### Show Y(A=135) plot ######
    # df[df['mass']== 135 ]['yieldWahl'].plot()
    # plt.show()

    print(df[(df["mass"] == 99)])
    print(df[(df["mass"] == 100)])

    print(df[(df["mass"] == 135)])
    print(df[(df["mass"] == 136)])

    # nums = ([2,2], [2,3], [3,2], [3,3])
    nums = ([50, 78], [50, 79], [51, 78], [51, 79])
    print("# Wahl               Minato         Titech")
    for n in nums:
        print(
            WahlEvenOddTerm(n[0], n[1]),
            MinatoEvenOddTerm(n[0], n[1]),
            TitechEvenOddTerm(n[0], n[1]),
        )


if __name__ == "__main__":
    main()
