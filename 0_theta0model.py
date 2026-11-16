#!/usr/bin/env python3

from florabaptistina import *

def main():

    theta0 = np.loadtxt("theta0.tmp")

    data = load_data("elmts.txt")

    model = ModelN(data, N=2, debug=True)

    model.fixed = theta0
#    model.debug = True

    lnp = model.lnlike(theta0)
#    lnp = model.lnprob(theta0)

    print(model)
    print("theta0 = ", theta0)
    print("lnp = ", lnp)
    model.save()

if __name__ == "__main__":
    main()


