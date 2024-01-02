## configs
# MASSFILE = "./data/mass-frdm95.dat"
MASSFILE = "/Users/okumuras/Dropbox/Development/zp_model/data/mass-frdm95.dat"
# MASSFILE = "./data/mass_1.mass20.txt"
LDFILE = "/Users/okumuras/Dropbox/Development/zp_model/data/kcksyst.dat"


## parameters
try:
    from .func import sep_nuclide
except:
    from func import sep_nuclide


CN = "U236"
ZCN, ACN = sep_nuclide(CN)

EIN = 5.0  # Incident energy in MeV
ENEUTRON = 8.071431  #  Neutron Mass Excess [MeV]

MIN_A = int(ACN / 2 - 50)  # assuming the a of FF distributes +/- 60a from symmetric
MAX_A = int(ACN / 2 + 50)

MIN_Z = int(ZCN / 2 - 30)  # assuming the z of FF distributes +/- 30z from symmetric
MAX_Z = int(ZCN / 2 + 30)


## Wahl systematics parameters
## not used
ZPFACTOR = [1.00000, 1.00000]

GAUSS1 = [0.79310, 4.82780, 22.99700]
GAUSS2 = [0.20542, 2.72770, 15.63400]
GAUSSSYM = [0.00295, 8.60000]

DELTAZ = 0.5
SIGMAZ = 0.5

## Yield cutoff
Y_CUTOFF = 1.0e-15

## Wahl ZP model parameter
VX20 = 92    # reference to U235+n reaction
VX21 = 236   # reference to U235+n reaction
VY58 = 8.0
VX58 = 20.0

## See details in 
## JNST 2018, 55, 9, 1009–1023 https://doi.org/10.1080/00223131.2018.1467288
## JNST 2021, https://doi.org/10.1080/00223131.2021.1954103
RT = 1.0

## set true if the constant TKE width is used and set TKE_WIDTH for that
TKE_WIDTH_CONST = False
TKE_WIDTH = 8.0 # in MeV if TKE_WIDTH_CONST is True

# if the width of TKE is given, set true
TKE_WIDTH_GIVEN = True

