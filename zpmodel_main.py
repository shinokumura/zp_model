from params import CN, ZCN, ACN
from ya import gen_massdist,read_massdist
from zp import fractionZp
from excite import gen_tke_a, fragment_excitation_energy
from func import mass_excess,n_ind_sep_en


def main(massmode="generate",filename=""):
    ## generate mass distribution
    if massmode == "generate":
        mass_df = gen_massdist(CN)
    if massmode == "read":
        if not filename:
            print("please set filename")
        else:
            mass_df = read_massdist(filename)
    # print(mass_df)

    ## generate fractional yield
    fractional_df = fractionZp(CN, mass_df, eo_model="wahl")
    # print(fractional_df)

    ## get average TKE(A)
    tke_ah_dict = gen_tke_a()
    # print(tke_ah_dict)

    fragment_dict = fragment_excitation_energy(tke_ah_dict, fractional_df)
    # print(fragment_dict)

    mass_ex_tn = mass_excess(ZCN, ACN - 1)  # only the case of the first chance fission
    mass_ex_cn = mass_excess(ZCN, ACN)
    
    # separation energy of fissioning system
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)

    ## output
    print("# Z = ", ZCN)
    print("# A = ", ACN)
    print("# Ex (MeV) = {:9.2E}".format(sep_cn))
    print("# Ntotal = ", len(fragment_dict.keys()))
    print("# Zl  Al Zh  Ah Yield      TKE[MeV]   TXE[MeV]   El[MeV]    Wl[MeV]    Eh[MeV]    Wh[MeV]")

    for k in reversed(fragment_dict.keys()):
        print("{:4d}{:4d}{:3d}{:4d}{:11.4E}{:11.4E}{:11.4E}{:11.4E}{:11.4E}{:11.4E}{:11.4E}".format(    
        int(fragment_dict[k]["zl"]),
        int(fragment_dict[k]["al"]),
        int(fragment_dict[k]["zh"]),
        int(fragment_dict[k]["ah"]),
        fragment_dict[k]["ffy"],
        fragment_dict[k]["ke"],
        fragment_dict[k]["ex"],
        fragment_dict[k]["ex_L"],
        fragment_dict[k]["ex_w_L"],
        fragment_dict[k]["ex_H"],
        fragment_dict[k]["ex_w_H"])
        )

if __name__ == "__main__":
    main()
