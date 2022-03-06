import pandas as pd
import math

def read_ld_data():
    from params import LDFILE
    from func import slices

    z = []
    a = []
    pairing = []
    eshell = []
    lda = []
    temp = []
    e0 = []

    with open(LDFILE) as f:
        lines = f.readlines()
        for line in lines[9:]:
            kz, ka, kpare, keshell, kematch, klda, kt, ke0, ktsys, kte0sys = slices(line, 6, 6, 13, 13, 10, 10, 10, 10, 10, 10)
            z.append(int(kz))
            a.append(int(ka))
            pairing.append(float(kpare))
            eshell.append(float(keshell))
            lda.append(float(klda))

            if float(kematch) > 0.0:
                temp.append(float(kt))
                e0.append(float(ke0))
            else:
                temp.append(float(ktsys))
                e0.append(float(kte0sys))

    df = pd.DataFrame({"Z": z, "A": a, "Paring": pairing, "Eshell":eshell, "LDa":lda, "Temp": temp, "E0": e0})

    ## check if data is read correctry
    if df.empty:
        raise TypeError()

    return df

ld_df = read_ld_data()

def ld_data(z, a):
    try:
        ld_data = ld_df.loc[(ld_df.Z == z) & (ld_df.A == a)].values.flatten().tolist()
    except:
        print('unknown level density nuclide: ', z, a)
    return ld_data



def gilbert_cameron_type_ld():
    pass


def ld_shell_correction(u, lda, shellcorrection, mass):
    gamma = 0.31 * (mass**(-1/3.0))
    ax = lda

    if (u > 0.0):
        ax = lda * (1.0 + (1.0-math.exp(-gamma*u))*shellcorrection/u)
    elif(u == 0.0):
        ax = lda * (1.0+gamma*shellcorrection) 

    return ax


def kckAsymptoticLevelDensity(a):
    return 0.126181 * a + 7.52191e-05 * a * a


if __name__ == "__main__":
    ld_data(30,  84)


