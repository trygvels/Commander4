"""Reproduce E005: per-pixel toy, scalar index, chain width and mean shift vs the exact posterior.
Uses the toy and chain helpers from tests/test_mcmc_likelihood.py unchanged."""
import sys, time
sys.path.insert(0, "tests")  # run from the repository root
import numpy as np
from test_mcmc_likelihood import PerPixelToy, _run_chain

print(f"{'npix':>5s} {'likelihood':14s} {'n_iter':>6s} {'sd ratio':>9s} {'mean shift/sd':>14s} {'time s':>7s}")
for npix, modes in [(150, ["marginal", "ridge", "conditional", "resolve_chisq"]),
                    (600, ["marginal", "ridge", "resolve_chisq"])]:
    for mode in modes:
        toy = PerPixelToy(npix=npix)
        mean, sd = toy.exact_posterior(use_logdet=(mode != "ridge"))
        n_iter = 2500 if mode in ("conditional", "resolve_chisq") else 1500
        t0 = time.time()
        samples = _run_chain(toy, mode, n_iter, 4, sigma=1.5 * sd, seed=7)
        kept = samples[n_iter // 5:]
        print(f"{npix:5d} {mode:14s} {n_iter:6d} {kept.std() / sd:9.2f} "
              f"{(kept.mean() - mean) / sd:+14.2f} {time.time() - t0:7.1f}", flush=True)
