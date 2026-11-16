#!/usr/bin/env python3

import pickle

from florabaptistina import *

def main():

    with open("output.pkl", "rb") as f:
        model, sampler = pickle.load(f)

    samples = sampler.get_chain(flat=True)
    lnp = sampler.get_log_prob(flat=True)

    steps = np.arange(0, sampler.iteration, 1)
    i = np.repeat(steps, sampler.nwalkers)
    j = np.tile(np.arange(sampler.nwalkers), len(steps))

    names = [
        "a_c^1", "w0^1", "C^1", "s_a^1", "beta_a^1", "mu_e^1", "w_e^1", "s_e^1", "beta_e^1", "mu_i^1", "w_i^1", "s_i^1", "beta_i^1", "Hb^1", "alpha1^1", "alpha2^1", "mu_p^1", "w_p^1", "s_p^1", "beta_p^1",
        "a_c^2", "w0^2", "C^2", "s_a^2", "beta_a^2", "mu_e^2", "w_e^2", "s_e^2", "beta_e^2", "mu_i^2", "w_i^2", "s_i^2", "beta_i^2", "Hb^2", "alpha1^2", "alpha2^2", "mu_p^2", "w_p^2", "s_p^2", "beta_p^2",
        "z1",
        ]

    names = np.asarray(names)
    names = names[model.free]

    np.savetxt("lnp.tmp", np.column_stack((i,j,samples,lnp)), header='i j '+' '.join(names)+' lnp', fmt=["%d","%d"]+len(names)*["%f"]+["%f"])

if __name__ == "__main__":
    main()

