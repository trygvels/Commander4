import numpy as np
COEF = (2.0 * np.pi) ** 0.25 / 2.0
def run(k, numstep, tau, nchain=4000, burn=300, keep=600, seed=5):
    tau = np.asarray(tau, float); dim = tau.size
    rng = np.random.default_rng([seed, k, numstep])
    beta = np.zeros((nchain, dim)); out = np.empty((keep, nchain, dim)); nacc = 0
    for it in range(burn + keep):
        X = rng.chisquare(k, nchain)
        for _ in range(numstep):
            prop = beta + tau * rng.standard_normal((nchain, dim))
            Xp = rng.chisquare(k, nchain)
            loga = -0.5 * ((prop ** 2).sum(1) - (beta ** 2).sum(1)) - 0.5 * (Xp - X)
            acc = np.log(rng.random(nchain)) < loga
            beta[acc] = prop[acc]; X[acc] = Xp[acc]
            if it >= burn: nacc += acc.sum()
        if it >= burn: out[it - burn] = beta
    h = keep // 2
    return out.std((0, 1)), out[:h].std((0, 1)), out[h:].std((0, 1)), nacc / (keep * numstep * nchain)
k = 190825; law = COEF * k ** 0.25
print(f"# E014 settings with the recorded proposal widths: 0.002 in beta_s (ridge sd 0.00036) and 0.004 in beta_d (ridge sd 0.000043),")
print(f"# i.e. 5.6 and 93 ridge sd; k = {k}, 20 steps per refresh; law (one step per refresh) {law:.3f}; measured 17.1 (beta_s), 16.5 (beta_d), acceptance 0.073")
sd, s1, s2, acc = run(k, 20, [0.002 / 0.00036, 0.004 / 0.000043])
print(f"toy sd beta_s {sd[0]:.3f} (halves {s1[0]:.3f}, {s2[0]:.3f}), beta_d {sd[1]:.3f} (halves {s1[1]:.3f}, {s2[1]:.3f}), acceptance {acc:.4f}")
sd, s1, s2, acc = run(k, 20, [0.002 / 0.00036, 0.004 / 0.000043], seed=6)
print(f"second seed: beta_s {sd[0]:.3f}, beta_d {sd[1]:.3f}, acceptance {acc:.4f}")
