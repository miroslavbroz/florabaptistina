#!/usr/bin/env python3

import pickle
import matplotlib.pyplot as plt
import corner

from florabaptistina import *

def main():

    with open("output.pkl", "rb") as f:
        model, sampler = pickle.load(f)

    names = [
        "a_c^1", "w0^1", "C^1", "s_a^1", "beta_a^1", "mu_e^1", "w_e^1", "s_e^1", "beta_e^1", "mu_i^1", "w_i^1", "s_i^1", "beta_i^1", "Hb^1", "alpha1^1", "alpha2^1", "mu_p^1", "w_p^1", "s_p^1", "beta_p^1",
        "a_c^2", "w0^2", "C^2", "s_a^2", "beta_a^2", "mu_e^2", "w_e^2", "s_e^2", "beta_e^2", "mu_i^2", "w_i^2", "s_i^2", "beta_i^2", "Hb^2", "alpha1^2", "alpha2^2", "mu_p^2", "w_p^2", "s_p^2", "beta_p^2",
        "z1",
        ]

    names = np.asarray(names)
    names = names[model.free]

    burn = 0
    burn = 500
    chain = sampler.get_chain(flat=True, discard=burn)
    log_prob = sampler.get_log_prob(flat=True, discard=burn)

    idx = np.argmax(log_prob)
    lnp = log_prob[idx]
    theta = chain[idx]

    fig = corner.corner(
        chain[:, :],
        quantiles=[0.16, 0.50, 0.84],
        truths=theta,
        truth_color='lime',
        labels=names,
        show_titles=True,
        title_fmt=".3g",
        title_kwargs={"fontsize": 10}
    )

    plt.savefig("corner_.png")

if __name__ == "__main__":
    main()


