import sys
sys.path.insert(0, "nextgen/scripts")
from e035_legacy_tempering import run, COEF
print("# P5 (added after P1 to P4 were read; not pre-registered): the toy at the settings of the earlier measurements")
print(f"{'case':42s} {'k':>7s} {'dim':>3s} {'nstep':>5s} {'tau':>7s} {'law':>7s} {'toy sd':>7s} {'+-':>6s} {'measured':>9s} {'acc':>6s}")
for label, k, dim, ns, tau, meas in [
    ("E002 prototype, 256 px x 3 comps", 768, 2, 1, 1.0 * COEF * 768 ** 0.25, "3.9"),
    ("E002 prototype, 256 px, wider step", 768, 2, 1, 2.0 * COEF * 768 ** 0.25, "3.9"),
    ("E002 prototype, 1024 px x 3 comps", 3072, 2, 1, 1.0 * COEF * 3072 ** 0.25, "5.7"),
    ("E002 prototype, 1024 px, wider step", 3072, 2, 1, 2.0 * COEF * 3072 ** 0.25, "5.7"),
    ("E005 unit test, 150 px x 3 comps", 450, 1, 4, 1.5, "3.2"),
    ("E005 unit test, 600 px x 3 comps", 1800, 1, 4, 1.5, "4.4"),
    ("E014 real maps, 38165 px x 5 comps", 190825, 2, 20, 2.4 * 16.8, "16.5/17.1"),
]:
    sd, se, s1, s2, ku, acc, _ = run(k, ns, tau, "chi2", dim=dim, burn=300, keep=600)
    print(f"{label:42s} {k:7d} {dim:3d} {ns:5d} {tau:7.2f} {COEF * k ** 0.25:7.3f} {sd:7.3f} {se:6.3f} {meas:>9s} {acc:6.3f}", flush=True)
