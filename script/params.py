## configs
MASSFILE = "data/mass-frdm95.dat"
# MASSFILE = "data/mass-hfb14.dat"
LDFILE = "data/kcksyst.dat"


## parameters
from func import sep_nuclide

CN = "U236"
ZCN, ACN = sep_nuclide(CN)

EIN = 2.83e-8  # incident energy in MeV
# EEX = EIN
ENEUTRON = 8.071431  #  Neutron Mass Excess [MeV]

MIN_A = int(ACN / 2 - 50)  # assuming the a of FF distributes +/- 60a from symmetric
MAX_A = int(ACN / 2 + 50)

MIN_Z = int(ZCN / 2 - 30)  # assuming the z of FF distributes +/- 30z from symmetric
MAX_Z = int(ZCN / 2 + 30)

## not used
ZPFACTOR = [1.00000, 1.00000]

GAUSS1 = [0.79310, 4.82780, 22.99700]
GAUSS2 = [0.20542, 2.72770, 15.63400]
GAUSSSYM = [0.00295, 8.60000]

Y_CUTOFF = 1.0e-15

DELTAZ = 0.5
SIGMAZ = 0.5
    
## See details in JNST 2018, 55, 9, 1009–1023 https://doi.org/10.1080/00223131.2018.1467288
RT = 1.25
TKE_WIDTH_CONST = False
TKE_WIDTH = 8.0 # in MeV
