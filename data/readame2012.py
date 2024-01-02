

from ..script.func import slices




#57     32   89   Ge         x      -33730#    600#     8169#       7#             β−    13070#     670#    88 963790#      640#
def a():

    z = []
    n = []
    a = z + n
    mselect = []
    with open("a") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith([0-9]):
                z, n, a, el, orig, mx, Bn,     = slices(line, 5, 5, 5, 5, 10, 10, 10)


    df = pd.DataFrame({"Z": z, "A": a, "Mselect": mselect})

    ## check
    if df.empty:
        raise TypeError()

    return df


mass_table_df = read_mass_table()