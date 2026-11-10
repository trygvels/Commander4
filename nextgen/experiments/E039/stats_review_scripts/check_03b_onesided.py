"""Check 3b: one-sided miss rates of the 1 sd and 2 sd intervals, from the rows saved by check_03.
Nominal: 0.1587 on each side at 1 sd, 0.0228 at 2 sd. 'low' = truth below the interval (beta_hat too flat)."""
import glob
import os

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
print(f"{'config':6s} {'interval':16s} {'low 1sd':>9s} {'high 1sd':>9s} {'low 2sd':>9s} {'high 2sd':>9s} {'se 1sd':>8s} {'se 2sd':>8s}")
for fn in sorted(glob.glob(os.path.join(here, "out_03_rows_*.npz"))):
    r = np.load(fn)
    name = os.path.basename(fn)[12:-4]
    for label, sd, sel in (("truth-sandwich", r["sd_sw"][None, :] * np.ones_like(r["err"]), np.ones_like(r["ok"])),
                           ("plugin", r["sd_plugin"], r["ok"])):
        z = (r["err"] / sd)[sel]
        m = z.size
        print(f"{name:6s} {label:16s} {np.mean(z > 1):9.3f} {np.mean(z < -1):9.3f} {np.mean(z > 2):9.4f} {np.mean(z < -2):9.4f} "
              f"{np.sqrt(0.1587 * 0.8413 / m):8.3f} {np.sqrt(0.0228 * 0.9772 / m):8.4f}")
