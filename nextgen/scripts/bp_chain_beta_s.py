"""Read the per-region synchrotron index values of every sample of a remote BeyondPlanck chain
file over HTTP range requests (E022: the BP XV beta_s chain itself). Usage:
  python -I bp_chain_beta_s.py URL OUT.npz [max_bytes] [max_samples]
"""
import sys, io, time, urllib.request, h5py, numpy as np

URL, OUT = sys.argv[1], sys.argv[2]
MAXBYTES = int(sys.argv[3]) if len(sys.argv) > 3 else 3_000_000_000
MAXS = int(sys.argv[4]) if len(sys.argv) > 4 else 100000


class HttpFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url; self.pos = 0; self.fetched = 0; self.ncalls = 0
        r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60)
        self.size = int(r.headers["Content-Length"])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]; return self.pos
    def readinto(self, b):
        n = len(b)
        if self.pos >= self.size: return 0
        end = min(self.pos + n, self.size) - 1
        if self.fetched + (end - self.pos + 1) > MAXBYTES:
            raise RuntimeError(f"byte budget exceeded ({self.fetched} fetched, {self.ncalls} requests)")
        for attempt in range(5):
            try:
                req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
                data = urllib.request.urlopen(req, timeout=120).read(); break
            except Exception as e:
                if attempt == 4: raise
                time.sleep(5 * (attempt + 1))
        self.fetched += len(data); self.ncalls += 1
        b[:len(data)] = data; self.pos += len(data); return len(data)


hf = HttpFile(URL)
print("file size bytes:", hf.size, flush=True)
f = h5py.File(io.BufferedReader(hf, buffer_size=1 << 16), "r")
t0 = time.time()
keys = list(f.keys())
samples = sorted(k for k in keys if k != "parameters")
print(f"{len(samples)} samples, first {samples[0]} last {samples[-1]}; fetched {hf.fetched/1e6:.1f} MB in {hf.ncalls} requests, {time.time()-t0:.0f}s", flush=True)
g = f[samples[0] + "/synch"]
print("synch datasets in", samples[0], ":", list(g.keys()), flush=True)
for k in g.keys():
    if "pixreg" in k or k.endswith("lmax") or k.endswith("nmaps"):
        print("  ", k, g[k].shape, g[k][()].ravel()[:12], flush=True)
vals = {}; names = [k for k in g.keys() if k.endswith("pixreg_val")]
for i, s in enumerate(samples[:MAXS]):
    try:
        for k in names: vals.setdefault(k, []).append(np.array(f[s + "/synch/" + k][()]))
        vals.setdefault("sample", []).append(int(s))
    except Exception as e:
        print("  sample", s, "failed:", e, flush=True); continue
    if i % 50 == 0:
        print(f"  {i+1} samples read, fetched {hf.fetched/1e6:.1f} MB in {hf.ncalls} requests, {time.time()-t0:.0f}s", flush=True)
        np.savez(OUT, **{k: np.array(v) for k, v in vals.items()})
np.savez(OUT, **{k: np.array(v) for k, v in vals.items()})
print("done:", {k: np.array(v).shape for k, v in vals.items()}, "fetched", hf.fetched / 1e6, "MB", flush=True)
