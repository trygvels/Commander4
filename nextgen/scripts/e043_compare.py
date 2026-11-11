"""E043: compare a recorded simulation output with its rerun on the sound Cholesky path, line by line.

Every line that starts with "real" is split into its numbers (the elapsed time in brackets is dropped);
lines are matched by realization index. Prints the largest absolute difference per kind of number, the
realizations that differ by more than the tolerance, and the summary lines of both files side by side.

    python -I nextgen/scripts/e043_compare.py ORIGINAL RERUN [--tol 2e-4]
"""
import argparse
import re

NUM = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")


def parse(path):
    rows, other = {}, []
    for ln in open(path):
        if ln.startswith("real "):
            body = re.sub(r"\[\s*\d+ s\]", "", ln)
            idx = int(body.split()[1])
            rows[idx] = ([float(x) for x in NUM.findall(body)], body.strip())
        elif ln.strip() and not ln.startswith("#"):
            other.append(ln.rstrip())
    return rows, other


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("original")
    ap.add_argument("rerun")
    ap.add_argument("--tol", type=float, default=2e-4)
    a = ap.parse_args()
    A, oa = parse(a.original)
    B, ob = parse(a.rerun)
    common = sorted(set(A) & set(B))
    print(f"original: {len(A)} realizations; rerun: {len(B)}; compared: {len(common)}; tolerance {a.tol:g}")
    worst, bad, nnum = 0.0, [], 0
    for i in common:
        xa, xb = A[i][0], B[i][0]
        if len(xa) != len(xb):
            bad.append((i, float("inf"), "different number of fields"))
            continue
        d = max((abs(p - q) for p, q in zip(xa, xb)), default=0.0)
        nnum += len(xa)
        worst = max(worst, d)
        if d > a.tol:
            bad.append((i, d, ""))
    print(f"numbers compared: {nnum}; largest absolute difference: {worst:.2e}; realizations beyond tolerance: {len(bad)}")
    for i, d, why in bad[:40]:
        print(f"  realization {i}: max difference {d:.3e} {why}")
        print("    original: " + A[i][1][:300])
        print("    rerun:    " + B[i][1][:300])
    ident = sum(A[i][1] == B[i][1] for i in common)
    print(f"lines identical as printed: {ident} of {len(common)}")
    if len(B) >= len(A):
        print("summary lines, original | rerun:")
        for la, lb in zip(oa, ob):
            flag = "" if la == lb else "   <-- differs"
            print("  " + la[:150])
            if la != lb:
                print("  " + lb[:150] + flag)


if __name__ == "__main__":
    main()
