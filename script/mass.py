import pandas as pd

try:
    from .func import slices
except:
    from func import slices



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
    # from func import slices
    try:
        from .params import MASSFILE
    except:
        from params import MASSFILE

    z = []
    a = []
    mselect = []
    with open(MASSFILE) as f:
        lines = f.readlines()
        if any(x in MASSFILE for x in ["mass-frdm95.dat", "mass-hfb14.dat"]):
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

        if "mass_1.mass20.txt" in MASSFILE:

            for line in lines[36:]:
                t_1, t_nz, t_n, t_z, t_a, t_el, t_0, t_mexcess = slices(line, 1, 3, 5, 5, 5, 4, 5, 14)
                z.append(int(t_z))
                a.append(int(t_a))

                if not t_mexcess.isspace():
                    mselect.append(float(t_mexcess.replace("#",""))/1000)
                else:
                    mselect.append(0.0)

    df = pd.DataFrame({"Z": z, "A": a, "Mselect": mselect})

    ## check
    if df.empty:
        raise TypeError()

    return df



def get_closest_mass(z, a):
    masslist = mass_table_df.loc[(mass_table_df.Z == z)]["A"].to_list()
    closest_mass = min(masslist, key=lambda x:abs(x-a))
    return closest_mass



mass_table_df = read_mass_table()


def mass_excess(z, a):

    try:
        mx = mass_table_df.loc[(mass_table_df.Z == z) & (mass_table_df.A == a)][
            "Mselect"
        ].to_string(index=False)
        return float(mx)
    
    except:
        closest_mass = get_closest_mass(z, a)
        mx = mass_table_df.loc[(mass_table_df.Z == z) & (mass_table_df.A == closest_mass)][
            "Mselect"
        ].to_string(index=False)

        print("# ", z, a, " not found in mass table, so chose closest data has been taken:", z, closest_mass, mx)
        
        return float(mx)
    
