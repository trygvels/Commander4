import numpy as np, sys
sys.argv = ["x"]
exec(open("/private/tmp/claude-501/-Users-trygvels-code-Commander4/e8b23887-d1a3-4d04-b2ac-228d7b739a81/scratchpad/f26_check.py").read().split("deltas = np.arange")[0])
deltas = np.arange(-0.10, 0.1001, 0.02)
rng = np.random.default_rng(1)
base = sky - cmb[:, None]*cm[None]
for label, sel, comp in (("d_Q2", qd == 1, "d"), ("d_Q3", qd == 2, "d"), ("s_Q4", qs == 3, "s")):
    for k in range(4):
        c2 = rng.permutation(cmb)
        scan(f"{label} cmb shuffled {k}", sel, base + c2[:, None]*cm[None], rms, 100., comp, deltas, THETA[100.])
    for off in (6000, 12000, 30000, 40000):
        c2 = np.roll(cmb, off)
        scan(f"{label} cmb rolled {off}", sel, base + c2[:, None]*cm[None], rms, 100., comp, deltas, THETA[100.])
