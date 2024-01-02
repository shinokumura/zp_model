import numpy as np
### Added to fit fission event data


def tke_talou_polynom(x):
    ## Taken from the CGMF
    def h(x):
        return x - 135.0

    return (
        175.94
        - 0.95418 * h(x)
        - 0.064178 * h(x) ** 2
        + 0.007506 * h(x) ** 3
        - 0.0003466 * h(x) ** 4
        + 7.0311e-06 * h(x) ** 5
        - 5.264e-08 * h(x) ** 6
    )

def tke_exp_fit(x,  a, b, c, d):
    ## optimised parameters for U-235 thermal neutron case, JNST, 55:9, 1009-1023, 2018
    # a               = 335.284
    # b               = 1.17404
    # c               = 0.187596
    # d               = 69.0834
    return (a - b * x) * (1 - c * np.exp(-((x - 118) ** 2) / d))


def tke_cubic(x, a, b, c, d):
    def h(x):
        return x - 135
    
    return (
        a
        - b * h(x)
        - c * h(x) ** 2
        - d * h(x) ** 3
    )


def tke_polynom(x, a, b, c, d, e, f, g):
    def h(x):
        return x - 135.0

    return (
        a
        - b * h(x)
        - c * h(x) ** 2
        + d * h(x) ** 3
        - e * h(x) ** 4
        + f * h(x) ** 5
        - g * h(x) ** 6
    )

