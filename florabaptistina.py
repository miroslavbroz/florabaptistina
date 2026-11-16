#!/usr/bin/env python3

"""
florabaptistina.py
MCMC for N asteroid families.
Miroslav Broz (miroslav Broz@email.cz), Oct 22nd 2026

With help of chatgpt.com:

In Python, implement the following:
read 4-dimensional data elmts == a,e,i,H by np.loadtxt()
define generalized skewed Gaussian function G2(), w. 5 parameters (a,b,c,d,e), i.e., centre, width, normalisation, skewness, exponent
for e, i, use these functions as dN(e), dN(i)
for H, use a broken power law as dN(H)
for a, use G2() as dN(a, H), but with variable width, because asteroid families are wider in a for higher H, as H > log10(|a-a_c|/C)
define model distribution as dN(a,H) * dN(e) * dN(i) * dN(H)
define chi2() function
define lnlike() function as likelihood
define lnprior() as uninformed priors
define unimformed priors 
define run_mcmc() to fit the observed distribution of elmts

The observational data:

\vec x = (a, e, i, H)

Model:

p(\vec x) = p(a|H) p(e) p(i) p(H)

Model2:

p(\vec x) = f p1 + (1-f) p2

Hence, the family membership:

                   f p1
p(1|\vec x) = ---------------
              f p1 + (1-f) p2

                 (1-f) p2
p(2|\vec x) = ---------------
              f p1 + (1-f) p2

"""

import numpy as np
from scipy.special import erf, gammaln, log_ndtr, logsumexp, softmax
from scipy.integrate import cumulative_trapezoid, quad
import emcee

# ============================================================
# Data
# ============================================================

def load_data(filename):
    """Read columns a, e, i, h, p (and "orig" -> the original catalogue line no.)."""
    data = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4, 5, 6), converters=lambda s: np.nan if s.strip() == '?' else float(s))
    data = np.atleast_2d(data)

    data = data[
        (data[:, 0] > 2.18)    # nu6 bias
        & (data[:, 0] < 2.33)
        & (data[:, 1] >= 0)
        & (data[:, 1] < 1)
        & (data[:, 2] >= 0)
        & (data[:, 2] < 1)
        & (data[:, 3] < 18.5)  # observational bias
    ]

    if len(data) < 20:
        raise ValueError("Too few valid asteroids.")

    return data.T

# ============================================================
# Special functions.
# ============================================================

def logG(x, a, b, c, d):
    """
    Log of normalized skewed Gaussian.

    a .. centre
    b .. width
    c .. skewness
    d .. exponent

        x-a
    z = ---
         b

            1            2
    G = ---------- exp(-z )
        b sqrt(pi)

             d               d
    G = ------------ exp(-|z| ) Phi(sqrt(2) c z)
        b Gamma(1/d)      

    """

    if b.any() <= 0 or d <= 0:
        return np.full_like(np.asarray(x, dtype=float), -np.inf)

    z = (x-a)/b

    return (np.log(d) - np.log(b) - gammaln(1.0/d) - np.abs(z)**d + log_ndtr(np.sqrt(2.0)*c*z))

def logG_(x, a, b, c, d, lo, hi):
    """Normalisation on <lo, hi>."""

    # Use only ~10 normalisations and then an interpolation => 2.2 iterations s^-1.
    b = np.asarray(b)
    b10 = np.quantile(b, np.linspace(0.0, 1.0, 50))
    Z10 = np.asarray([quad(lambda x: np.exp(logG(x, a, b1, c, d)), lo, hi)[0] for b1 in b10])
    idx = (Z10 > 0.0)
    Z = np.interp(b, b10[idx], Z10[idx])
    Z = np.asarray(Z)

    if not np.isfinite(Z.any()) or Z.any() <= 0:
        return np.full_like(np.asarray(x, dtype=float), -np.inf)

    return logG(x, a, b, c, d) - np.log(Z)

def int_exp_kx(k, lo, hi):
    """Integral of exp(k x) between lo and hi; stable for k -> 0."""

    if hi <= lo:
        return 0.0

    if abs(k) < 1e-10:
        return hi - lo

    return np.exp(k*lo) * np.expm1(k*(hi - lo)) / k

def logH_normalization(Hb, alpha1, alpha2, lo, hi):
    """Normalization on <lo, hi>."""

    k1 = np.log(10.0)*alpha1
    k2 = np.log(10.0)*alpha2

    Z1 = int_exp_kx(k1, lo - Hb, min(hi, Hb) - Hb)
    Z2 = int_exp_kx(k2, max(lo, Hb) - Hb, hi - Hb)

    Z = Z1 + Z2

    if not np.isfinite(Z) or Z <= 0:
        return np.inf

    return np.log(Z)

def logH(H, Hb, alpha1, alpha2, lo, hi):
    """Broken power law."""
    H = np.asarray(H)

    logZ = logH_normalization(Hb, alpha1, alpha2, lo, hi)

    if not np.isfinite(logZ):
        return np.full_like(H, -np.inf, dtype=float)

    alpha = np.where(H <= Hb, alpha1, alpha2)

    return (np.log(10.0)*alpha*(H - Hb) - logZ)

def _sample(logpdf, lo, hi, n, rng):
    if n == 0:
        return np.empty(0)

    x = np.linspace(lo, hi, n)
    lp = np.asarray(logpdf(x), dtype=float)
    p = np.exp(lp - np.max(lp))

    cdf = cumulative_trapezoid(p, x, initial=0.0)

    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Invalid PDF or sampling interval.")

    cdf /= cdf[-1]
    cdf, idx = np.unique(cdf, return_index=True)

    return np.interp(rng.random(n), cdf, x[idx])

# ============================================================
# Model
# ============================================================

# Parameter order:
#
#  0 ac       centre of a distribution
#  1 w0       a minimum width
#  2 C        magnitude-dependent width scale
#  3 sa       a skewness
#  4 beta_a   a exponent
#  5 mu_e     e centre
#  6 w_e      e width
#  7 s_e      e skewness
#  8 beta_e   e exponent
#  9 mu_i     i centre
# 10 w_i      i width
# 11 s_i      i skewness
# 12 beta_i   i exponent
# 13 Hb       H break
# 14 alpha1   bright-side slope
# 15 alpha2   faint-side slope
# 16 mu_p     pv centre
# 17 w_p      pv width
# 18 s_p      pv skewness
# 19 beta_p   pv exponent

class Model:

    def __init__(self, data, tmp='lnlike.tmp', debug=False):
        self.a = np.asarray(data[0])
        self.e = np.asarray(data[1])
        self.i = np.asarray(data[2])
        self.h = np.asarray(data[3])
        self.p = np.asarray(data[4])
        self.psigma = np.asarray(data[5])
        self.tmp = tmp
        self.debug = debug

        self.amin, self.amax = np.min(self.a), np.max(self.a)
        self.emin, self.emax = np.min(self.e), np.max(self.e)
        self.imin, self.imax = np.min(self.i), np.max(self.i)
        self.hmin, self.hmax = np.min(self.h), np.max(self.h)
        self.pmin, self.pmax = np.min(self.p), np.max(self.p)

        self.arange = max(self.amax - self.amin, 1e-8)
        self.erange = max(self.emax - self.emin, 1e-8)
        self.irange = max(self.imax - self.imin, 1e-8)

        self.bounds = [
            (self.amin, self.amax),	# a_c
            (0.0, 0.200),		# w0
            (0.0, 3.0e-4),		# C
            (-10.0, 10.0),		# s_a
            (1.1, 6.5),			# beta_a
            (0.0, 1.0),			# mu_e
            (0.001, 0.100),		# w_e
            (-10.0, 10.0),		# s_e
            (1.1, 6.5),			# beta_e
            (0.0, 1.0),			# mu_i
            (0.001, 0.100),		# w_i
            (-10.0, 10.0),		# s_i
            (1.1, 6.5),			# beta_i
            (13.5,17.5),		# Hb
            (0.35, 0.7),		# alpha1
            (0.0, 0.7),			# alpha2
            (0.0, 1.0),			# mu_p
            (0.001, 1.0),		# w_p
            (-10.0, 10.0),		# s_p
            (1.1, 6.5),			# beta_p
        ]

        self.names = ["a_c", "w0", "C", "s_a", "beta_a", "mu_e", "w_e", "s_e", "beta_e", "mu_i", "w_i", "s_i", "beta_i", "Hb", "alpha1", "alpha2", "mu_p", "w_p", "s_p", "beta_p"]

        self.ndim = len(self.bounds)

    def lnlike(self, theta):
        """
        p(\vec x|theta)

        """
        (a_c, w0, C, s_a, beta_a, mu_e, w_e, s_e, beta_e, mu_i, w_i, s_i, beta_i, Hb, alpha1, alpha2, mu_p, w_p, s_p, beta_p) = theta
        self.theta = theta

        # Magnitude-dependent semimajor-axis width.
        w_a = np.sqrt(w0**2 + (C * 10.0**(0.2*(self.h-0.0)))**2)

        # Conditional probabilty p(a|H).
        log_a = logG_(self.a, a_c, w_a, s_a, beta_a, self.amin, self.amax)

        # Marginal eccentricity and inclination probabilities.
        log_e = logG(self.e, mu_e, w_e, s_e, beta_e)
        log_i = logG(self.i, mu_i, w_i, s_i, beta_i)

        # Magnitude distribution, normalized on observed H range.
        log_h = logH(self.h, Hb, alpha1, alpha2, self.hmin, self.hmax)

        # Geometric albedo (might be undefined, is uncertain)
        mask = np.isfinite(self.p)
        log_p = np.zeros_like(self.p)
        w_ = np.sqrt(w_p**2 + self.psigma[mask]**2)
        log_p[mask] = logG(self.p[mask], mu_p, w_, s_p, beta_p)

        lnp = log_a + log_e + log_i + log_h + log_p

        self.log_a = log_a
        self.log_e = log_e
        self.log_i = log_i
        self.log_h = log_h
        self.log_p = log_p

        if not np.all(np.isfinite(lnp)):
# Note: Ditto.
#            return -np.inf
            return np.full_like(lnp, -np.inf, dtype=float)

# Note: The sum will be done in ModelN.lnlike().
#        return np.sum(lnp)
        return lnp

    def lnprior(self, theta):
        """
        p(theta)

        """
        for val, (lo, hi) in zip(theta, self.bounds):
            if not np.isfinite(val) or not lo <= val <= hi:
                print("lnprior(): val = ", val, " lo = ", lo, " hi = ", hi)
                return -np.inf

        return 0.0

    def lnprob(self, theta):
        """
        p(theta|\vec x)

        """
        lp = self.lnprior(theta)

        if not np.isfinite(lp):
            return -np.inf

        return lp + self.lnlike(theta)

    def minuslnprob(self, theta):
        lp = self.lnprob(theta)

        return -lp if np.isfinite(lp) else 1.e100

    def __str__(self):
        s = "Model:\n"
        for i in range(len(self.names)):
            s = s + f"{i:2d} .. {self.names[i]:10s} = {self.theta[i]}\n"
        return s[:-1]

    def sample(self, theta, n=1, rng=None):

        rng = np.random.default_rng() if rng is None else rng

        (a_c, w0, C, s_a, beta_a, mu_e, w_e, s_e, beta_e, mu_i, w_i, s_i, beta_i, Hb, alpha1, alpha2) = theta

        e = _sample(lambda x: logG(x, mu_e, w_e, s_e, beta_e), 0.0, 1.0, n, rng)
        i = _sample(lambda x: logG(x, mu_i, w_i, s_i, beta_i), 0.0, 1.0, n, rng)
        p = _sample(lambda x: logG(x, mu_p, w_p, s_p, beta_p), 0.0, 1.0, n, rng)
        h = _sample(lambda x: logH(x, Hb, alpha1, alpha2, self.hmin, self.hmax), self.hmin, self.hmax, n, rng)

        w_a = np.sqrt(w0**2 + (C * 10.0**(0.2 * h))**2)

        a = _sample(lambda x: logG_(x, a_c, w_a, s_a, beta_a, self.amin, self.amax), self.amin, self.amax, n, rng)

        return a, e, i, h, p

# ============================================================
# N Models
# ============================================================

class ModelN(Model):

    def __init__(self, data, N=2, debug=False):
        self.orig = data[6]
        self.N = N
        self.debug = debug

        self.models = [Model(data) for _ in range(N)]

        self.bounds = sum((model.bounds for model in self.models), [])
        self.bounds += (N-1)*[(-5, 5)]
        self.bounds = np.asarray(self.bounds)

        self.ndim = len(self.bounds)
        self.ndims = [model.ndim for model in self.models]

        self.free = np.ones(self.ndim, dtype=bool)
        self.fixed = None
        self.w = np.full(N, 1/N)

    def expand(self, theta):
        theta_ = self.fixed.copy()
        theta_[self.free] = theta
        return theta_

    def expandN(self, theta):
        theta_ = self.expand(theta)
        N = self.N

        # Per-model parameters
        offsets = np.cumsum([0] + self.ndims)
        thetas = [theta_[offsets[i]:offsets[i+1]] for i in range(N) ]

        # N-1 logits; z_N = 0.0
        z = theta_[offsets[N]:]
        z = np.r_[z, 0.0]
        self.z = z.copy()
        z -= np.max(z)
        w = np.exp(z)
        w /= w.sum()
        self.w = w

        return thetas, w

    def lnlike(self, theta):
        """
        p(\vec x|theta)

        """
        thetas, w = self.expandN(theta)

        lnlikes = np.array([
            np.log(w[i]) + model.lnlike(thetas[i])
            for i, model in enumerate(self.models)
        ])

        lp = logsumexp(lnlikes, axis=0)

        self._lnlikes = lnlikes

        return np.sum(lp)

    def save(self):
        """
                           f p1
        p(1|\vec x) = ---------------
                      f p1 + (1-f) p2
        """
        lps = self._lnlikes
        lp = logsumexp(lps, axis=0)

        p_of_families = np.exp(lps - lp)
        i = range(len(p_of_families[0]))
        j = np.argmax(p_of_families.T, axis=1)
        k = np.array([ np.random.choice(len(p), p=p) for p in p_of_families.T ])
        model = self.models[0]
        elmts = np.array([model.a, model.e, model.i, model.h, model.p])
        np.savetxt("families.tmp",
            np.column_stack((self.orig, j, k, p_of_families.T, elmts.T)),
            header='orig j k_random p_of_families .. a e inc H p',
            fmt=3*['%d'] + self.N*['%f'] + 5*['%f'])

        for attr in ["a", "e", "i", "h", "p"]:
            x = getattr(self.models[0], attr)
            lps = [np.log(self.w[i]) + getattr(model, "log_" + attr) for i, model in enumerate(self.models)]
            lp = logsumexp(lps, axis=0)
            print("log_" + attr + " = ", np.sum(lp))
            np.savetxt(f"lnlike{attr}.tmp", np.column_stack([x] + lps + [lp]))

    def __str__(self):
        models = "\n".join(str(model) for model in self.models)
        return f"""
==
ModelN:
N = {self.N}
{models}
logits z = {self.z}
weights w = {self.w}
parameters free = {self.free}
dimensions ndim = {self.ndim}
==
"""

    def sample(self, theta, n=10000, seed=None):
        """Generate n samples from the N-component mixture."""
        rng = np.random.default_rng(seed)
    
        thetas, w = self.expandN(theta)
    
        family = rng.choice(self.N, size=n, p=w)
    
        a = np.empty(n)
        e = np.empty(n)
        i = np.empty(n)
        h = np.empty(n)
        p = np.empty(n)
    
        for k, model in enumerate(self.models):
            idx = np.flatnonzero(family==k)
            if len(idx) == 0:
                continue
    
            samples = model.sample(thetas[k], n=len(idx), rng=rng)
            a[idx], e[idx], i[idx], h[idx], p[idx] = samples
    
        return a, e, i, h, p, family

# ============================================================
# MCMC
# ============================================================

def run_mcmc(model, theta0, nwalkers=72, nsteps=10000, seed=1234):
    """
    Run the MCMC.

    """
    ndim = np.sum(model.free)
    bounds = np.asarray(model.bounds)[model.free]
    theta0 = theta0[model.free]
    rng = np.random.default_rng(seed)

    # Initialize walkers uniformly inside the prior bounds.
    #p0 = np.array([ rng.uniform(lo, hi, size=nwalkers) for lo, hi in model.bounds ]).T

    # Initial walker dispersion: 0.1% of each prior range
    scale = bounds[:, 1] - bounds[:, 0]
    p0 = theta0 + rng.normal(0, 1e-3 * scale, size=(nwalkers, ndim))

    for p in p0:
        if not np.isfinite(model.lnprob(p)):
            raise ValueError("Invalid initial walker. Check theta0 and its dispersion.")

    sampler = emcee.EnsembleSampler(nwalkers, ndim, model.lnprob)

    sampler.run_mcmc(p0, nsteps, progress=True)

    return sampler

# ============================================================
# Run
# ============================================================

def main():

    data = load_data("elmts.txt")

    model = ModelN(data, N=2)

# 0  .. 20 .. a_c   
# 1  .. 21 .. w0    
# 2  .. 22 .. C     
# 3  .. 23 .. s_a   
# 4  .. 24 .. beta_a
# 5  .. 25 .. mu_e  
# 6  .. 26 .. w_e   
# 7  .. 27 .. s_e   
# 8  .. 28 .. beta_e
# 9  .. 29 .. mu_i  
# 10 .. 30 .. w_i   
# 11 .. 31 .. s_i   
# 12 .. 32 .. beta_i
# 13 .. 33 .. Hb    
# 14 .. 34 .. alpha1
# 15 .. 35 .. alpha2
# 16 .. 36 .. mu_p  
# 17 .. 37 .. w_p   
# 18 .. 38 .. s_p   
# 19 .. 39 .. beta_p

# 40 .. z1

    theta0 = np.loadtxt("theta0.tmp")
#    theta0 = np.loadtxt("thetabest.tmp1__103978")

#    model.free[:] = False
#    model.free[[64,65,66]] = True  # only weights
#    model.free[[0,5,9,  16,21,25, 32,37,41, 48,53,57,  64,65,66]] = True  # only centres+weights
#    model.free[[0,1,5,6,9,10,  16,17,21,22,25,26, 32,33,37,38,41,42, 48,49,53,54,57,58,  64,65,66]] = True  # centres+widths+weights
#    model.free[[0,1,5,6,9,10,13,14,15,  16,17,21,22,25,26,29,30,31, 32,33,37,38,41,42,45,46,47, 48,49,53,54,57,58,  64,65,66]] = True  # centres+widths+magnitudes+weights
#    model.free[[0,20, 1,21, 2,22, 36,37,38,39, 40]] = False  # a_c, w0, C, albedo of B., z1
    model.fixed = theta0

    # per-model bounds, to prevent <->
    # a_c
    model.bounds[0] = (2.20-0.005, 2.20+0.005)
    model.bounds[20] = (2.26-0.005, 2.26+0.005)

    model.bounds[1] = (0.0, 0.005)
    model.bounds[21] = (0.0, 0.005)

    model.bounds[2] = (2.00e-4, 2.20e-4)
    model.bounds[22] = (0.19e-4, 0.21e-4)

    model.bounds[3] = (-0.5, 0.5)
    model.bounds[23] = (-0.5, 0.5)

    model.bounds[4] = (1.5, 2.5)
    model.bounds[24] = (5.5, 6.5)

    # mu_e
    model.bounds[5] = (0.12, 0.17)
    model.bounds[25] = (0.145, 0.150)

    model.bounds[6] = (0.020, 0.040)
    model.bounds[26] = (0.008, 0.012)

    model.bounds[7] = (-1.5, 0.5)
    model.bounds[27] = (-1.0, 0.0)

    model.bounds[8] = (3.0, 5.0)
    model.bounds[28] = (2.5, 3.5)

    # mu_i
    model.bounds[9] = (0.06,0.11)
    model.bounds[29] = (0.09,0.11)

    model.bounds[10] = (0.020,0.040)
    model.bounds[30] = (0.010,0.015)

    model.bounds[11] = (-0.5, 0.5)
    model.bounds[31] = (0.0, 1.0)

    model.bounds[12] = (1.5, 2.5)
    model.bounds[32] = (2.5, 3.5)

    # Hb
#    model.bounds[13] = (13.5,17.5)
#    model.bounds[33] = (13.5,17.5)

#    model.bounds[14] = (0.35, 0.7)
#    model.bounds[34] = (0.35, 0.7)

#    model.bounds[15] = (0.0, 0.7)
#    model.bounds[35] = (0.0, 0.7)

    # mu_A
    model.bounds[16] = (0.20,0.25)
    model.bounds[36] = (0.13,0.17)

    model.bounds[17] = (0.15,0.25)
    model.bounds[37] = (0.02,0.05)

    model.bounds[18] = (4.0,5.0)
    model.bounds[38] = (-0.1,0.1)

    model.bounds[19] = (1.5,2.5)
    model.bounds[39] = (2.5,3.5)

    # Z_1
    model.bounds[40] = (1.2,1.4)

    model.lnlike(theta0[model.free])
    print(model)

    sampler = run_mcmc(model, theta0, nwalkers=41*2, nsteps=1500)

    print(model)

    import pickle
    with open("output.pkl", "wb") as f:
        pickle.dump((model, sampler), f)

if __name__ == "__main__":
    main()


