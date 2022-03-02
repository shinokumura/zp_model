
import sys
sys.path.append('./script')

from script.params import CN, ZCN, ACN
from script.ya import gen_massdist,read_massdist
from script.zp import fractionZp
from script.excite import gen_tke_a, fragment_excitation_energy
from script.func import mass_excess,n_ind_sep_en

if len(sys.argv) > 1:
    massmode = sys.argv[1]
    filename = sys.argv[2]
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
    ## option model: wahl_zp, fixed
    ##        eo_model: if model = fixed, choose the even-odd term model
    fractional_df = fractionZp(mass_df, model = "wahl_zp", eo_model="")

    ## get average TKE(A)
    tke_ah_dict = gen_tke_a()

    ## energy divide into two fragments
    fragment_dict = fragment_excitation_energy(tke_ah_dict, fractional_df)

    ## separation energy of fissioning system
    mass_ex_tn = mass_excess(ZCN, ACN - 1)
    mass_ex_cn = mass_excess(ZCN, ACN)
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

