"""Referee check (draft 3): independent comparison of the exposed E039 outputs with their E043 reruns.
Strips the elapsed-time bracket and compares every remaining token, as text and as numbers."""
import re, sys
pairs = [
    ("experiments/E039/sim_n8_bpxv.txt", "experiments/E043/sim_n8_bpxv_rerun.txt"),
    ("experiments/E039/plugin_n8_bpxv.txt", "experiments/E043/plugin_n8_bpxv_rerun.txt"),
    ("experiments/E040/bpxv_n8_wm_all.txt", "experiments/E043/bpxv_n8_wm_all_rerun.txt"),
    ("experiments/E040/bpxv_n16_wm.txt", "experiments/E043/bpxv_n16_wm_oldpath.txt"),
    ("experiments/E040/bpxv_n16_wm.txt", "experiments/E043/bpxv_n16_wm_newpath.txt"),
    ("experiments/E043/bpxv_n16_wm_oldpath.txt", "experiments/E043/bpxv_n16_wm_newpath.txt"),
    ("experiments/E039/sim_n8_hp2_part1.txt", "experiments/E043/sim_n8_hp2_part1_rerun.txt"),
]
NUM = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")
def clean(l):
    l = re.sub(r"\[\s*\d+ s\]", "", l)
    l = re.sub(r";?\s*\d+ s\b", "", l)
    return l.rstrip()
for a, b in pairs:
    try:
        A = [clean(l) for l in open(a)]; B = [clean(l) for l in open(b)]
    except FileNotFoundError as e:
        print("missing", e); continue
    n = min(len(A), len(B))
    ident = sum(x == y for x, y in zip(A, B))
    worst = 0.0; nnum = 0; ndiff_lines = []
    for i, (x, y) in enumerate(zip(A, B)):
        if x != y:
            fx = [float(t) for t in NUM.findall(x)]; fy = [float(t) for t in NUM.findall(y)]
            if len(fx) == len(fy):
                d = max((abs(p - q) for p, q in zip(fx, fy)), default=0.0)
            else:
                d = float("inf")
            worst = max(worst, d if d != float("inf") else 0)
            ndiff_lines.append((i + 1, d))
        nnum += len(NUM.findall(x))
    print(f"{a}\n  vs {b}\n  lines {len(A)} / {len(B)}; identical after removing timings: {ident} of {n}; numbers on compared lines: {nnum}; largest numeric difference on differing lines: {worst:.2e}")
    for i, d in ndiff_lines[:12]:
        print(f"    line {i}: max diff {d:.2e}\n      A: {A[i-1][:230]}\n      B: {B[i-1][:230]}")
