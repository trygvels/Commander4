# Spurious-map chi-square: compare code's sum_j S_j^2 / Cov_jj against the exact
# likelihood-ratio statistic S^T Cov_SS^{-1} S for one pixel, 4 detectors (3 spurious cols).
import numpy as np
rng = np.random.default_rng(2)
ndet, nhit = 4, 30
for trial in range(3):
    rows, w, y = [], [], []
    truth = np.r_[10., 2., -1., 0.3, -0.2, 0.5, 0.]   # T Q U S1 S2 S3 (S4=0 by construction)
    for det in range(ndet):
        for t in range(nhit):
            psi = rng.uniform(0, np.pi)
            r = np.zeros(ndet-1+3); r[:3] = [1, np.cos(2*psi), np.sin(2*psi)]
            if det < ndet-1: r[3+det] = 1
            rows.append(r); w.append(1/rng.uniform(0.5,1.5)**2)
    P = np.array(rows); W = np.diag(w)
    d = P@truth[:6] if False else P@np.r_[truth[:3], truth[3:6]] + rng.standard_normal(len(w))/np.sqrt(w)
    A = P.T@W@P; Cov = np.linalg.inv(A); xh = Cov@P.T@W@d
    S = xh[3:]; CSS = Cov[3:,3:]
    code = np.sum(S**2/np.diag(CSS))
    exact = S@np.linalg.solve(CSS, S)
    # chi2(S=0 fit) - chi2(free fit)
    P0 = P[:, :3]; x0 = np.linalg.solve(P0.T@W@P0, P0.T@W@d)
    lr = (d-P0@x0)@W@(d-P0@x0) - (d-P@xh)@W@(d-P@xh)
    print(f"trial {trial}: code(diag)={code:.3f}  exact S^T C^-1 S={exact:.3f}  LR stat={lr:.3f}")
