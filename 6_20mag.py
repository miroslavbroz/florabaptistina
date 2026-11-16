#!/usr/bin/env python3

from florabaptistina import *

def load_data(filename):
    """Read columns a, e, i, h, p (and "orig" -> the original catalogue line no.)."""
    data = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4, 5, 6), converters=lambda s: np.nan if s.strip() == '?' else float(s))
    data = np.atleast_2d(data)

#    data = data[
#        (data[:, 0] > 2.18)    # nu6 bias
#        & (data[:, 0] < 2.33)
#        & (data[:, 1] >= 0)
#        & (data[:, 1] < 1)
#        & (data[:, 2] >= 0)
#        & (data[:, 2] < 1)
#        & (data[:, 3] < 18.5)  # observational bias
#    ]

    return data.T

def main():

    data = load_data("elmts.txt")

    model = ModelN(data, N=2)

    theta0 = np.loadtxt("theta0.tmp")

    model.fixed = theta0
    for _ in model.models: _.hmax = 18.5

    lnp = model.lnlike(theta0)
#    lnp = model.lnprob(theta0)

    print(model)
    print("theta0 = ", theta0)
    print("lnp = ", lnp)
    model.save()

if __name__ == "__main__":
    main()


