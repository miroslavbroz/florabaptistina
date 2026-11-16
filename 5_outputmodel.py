#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import corner

from florabaptistina import *

def main():

    with open("output.pkl", "rb") as f:
        model, sampler = pickle.load(f)

    log_prob = sampler.get_log_prob(flat=True)
    chain = sampler.get_chain(flat=True)

    idx = np.argmax(log_prob)
    lnp = log_prob[idx]
    theta = chain[idx]
    print("best lnp = ", lnp)
    print("best theta = ", theta)

    lnp_ = np.median(log_prob)
    theta_ = np.median(chain, axis=0)
    print("median lnp = ", lnp_)
    print("median theta = ", theta_)

    lnp__ = np.median(sampler.get_log_prob()[-1], axis=0)
    theta__ = np.median(sampler.get_chain()[-1], axis=0)
    print("last lnp = ", lnp__)
    print("last theta = ", theta__)
    print("")

    theta___ = model.expand(theta)
    np.savetxt("thetabest.tmp", theta___.T)

    model.debug = True
    lnp = model.lnlike(theta)

    print(model)
    print("theta = ", theta)
    print("lnp = ", lnp)
    model.save()

if __name__ == "__main__":
    main()


