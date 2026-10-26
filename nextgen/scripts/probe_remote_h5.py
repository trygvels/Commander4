import sys, io, urllib.request, h5py, numpy as np
URL = sys.argv[1]
MAXBYTES = int(sys.argv[2]) if len(sys.argv) > 2 else 40_000_000

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
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(req, timeout=120).read()
        self.fetched += len(data); self.ncalls += 1
        b[:len(data)] = data; self.pos += len(data); return len(data)

hf = HttpFile(URL)
print("file size bytes:", hf.size)
f = h5py.File(io.BufferedReader(hf, buffer_size=1 << 20), "r")
top = list(f.keys())
print("top-level keys:", len(top), top[:5], "...", top[-5:])
samples = [k for k in top if k != "parameters"]
print("nsamples:", len(samples), "first:", samples[0], "last:", samples[-1])
if "parameters" in f:
    print("=== parameters group ===")
    for comp, g in f["parameters"].items():
        d = {}
        for k, v in g.items():
            val = v[()]
            if isinstance(val, bytes): val = val.decode()
            d[k] = val
        print(" ", comp, d)
s0 = f[samples[0]]
print("=== groups in sample", samples[0], ":", list(s0.keys()))
def walk(name, obj):
    if isinstance(obj, h5py.Dataset):
        print("   ", name, obj.shape, obj.dtype, end="")
        if obj.size <= 4:
            print("  value=", obj[()], end="")
        print()
for comp in s0.keys():
    print(" --", comp)
    s0[comp].visititems(walk)
print("bytes fetched:", hf.fetched, "requests:", hf.ncalls)
