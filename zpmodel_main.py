
import sys
# sys.path.append('./script')

from script.params import CN, ZCN, ACN
from script.ya import gen_massdist,read_massdist
from script.zp import fractionZp
from script.excite import gen_tke_a, fragment_excitation_energy
from script.func import mass_excess, n_ind_sep_en

if len(sys.argv) > 1:
    massmode = sys.argv[1]
    filename = sys.argv[2]


def main(massmode="generate", filename=None):

    # define output format
    output_mode = "BeoH"

    if massmode == "generate":
        ## generate mass distribution from gaussian parameters in params.py
        mass_df = gen_massdist(CN)

    if massmode == "read":
        ## read mass distribution from file
        if not filename:
            print("please set filename")
        else:
            mass_df = read_massdist(filename)

    ## generate fractional yield
    ## option model: wahl_zp, fixed
    ##     eo_model: if "fixed" model is chosen, choose the even-odd term model
    fractional_df = fractionZp(mass_df, model = "wahl_zp", eo_model="")
    # fractional_df = fractionZp(mass_df, model = "fixed", eo_model="wahl")

    ## get average TKE(A)
    tke_ah_dict = gen_tke_a()

    ## energy divide into two fragments
    fragment_dict = fragment_excitation_energy(tke_ah_dict, None, fractional_df)

    ## calculate separation energy of fissioning system
    mass_ex_tn = mass_excess(ZCN, ACN - 1)
    mass_ex_cn = mass_excess(ZCN, ACN)
    sep_cn = n_ind_sep_en(mass_ex_tn, mass_ex_cn)

    ## output fission fragment distribution TALYS mode
    if output_mode == "TALYS":
        print("# Z        =   ", ZCN)
        print("# A        =  ", ACN)
        print("# Ex (MeV) = {:9.2E}".format(sep_cn))
        print("# Ntotal   =  ", len(fragment_dict.keys()))
        print("#  Zl  Al  Zh  Ah   Yield       TKE[MeV]    TXE[MeV]    El[MeV]     Wl[MeV]     Eh[MeV]     Wh[MeV]")
        i = 0
        for k in reversed(fragment_dict.keys()):
            print("{:4d}{:4d}{:4d}{:4d}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}".format(  
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
        print("")


    if output_mode == "BeoH":
        ## output fission fragment distribution BeoH mode
        print("# Z        =   ", ZCN)
        print("# A        =  ", ACN)
        print("# Ex (MeV) = {:9.2E}".format(sep_cn))
        print("# Ntotal   =  ", len(fragment_dict.keys()))
        print("#   n    k  Zl  Al  Zh  Ah  YieldL      YieldH      TKE[MeV]    TXE[MeV]    El[MeV]     Wl[MeV]     Eh[MeV]     Wh[MeV]")
        i = 0
        for k in reversed(fragment_dict.keys()):
            i += 1
            print("{:4d}{:4d}{:4d}{:4d}{:4d}{:4d}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}{:12.4E}".format(  
            i,
            i,  
            int(fragment_dict[k]["zl"]),
            int(fragment_dict[k]["al"]),
            int(fragment_dict[k]["zh"]),
            int(fragment_dict[k]["ah"]),
            fragment_dict[k]["ffy"],
            fragment_dict[k]["ffy"],
            fragment_dict[k]["ke"],
            fragment_dict[k]["ex"],
            fragment_dict[k]["ex_L"],
            fragment_dict[k]["ex_w_L"],
            fragment_dict[k]["ex_H"],
            fragment_dict[k]["ex_w_H"])
            )
        print("")


if __name__ == "__main__":
    main()
    # main(massmode="read", filename="/Users/okumuras/Dropbox/2022/TALYS-Langevin/Pu240")

