"""Read-only: load the E014 bands through Commander4's own reader (same degrade, units, smoothing
as the chain saw) and cache full-sky nside-64 maps and inverse variances in the scratchpad.
Stops after a time budget; rerun until all bands are cached."""
import os, sys, time
import numpy as np
import healpy as hp
from commander4.file_io.map_reader import read_data_map_from_file
from commander4.parameters.parse import load_params
from commander4.compsep.mcmc import resolve_chisq_mask

pf, cache, budget = sys.argv[1], sys.argv[2], float(sys.argv[3])
t0 = time.time()
params, _, _ = load_params(pf)
bands = params.compsep.bands
names = [n for n in bands if bands[n].enabled]
todo = 0
for name in names:
    out = os.path.join(cache, f"band_{name}.npz")
    if os.path.exists(out):
        continue
    if time.time() - t0 > budget:
        todo += 1
        continue
    b = bands[name]
    if "_name" not in b:
        b._name = name
    t1 = time.time()
    dm = read_data_map_from_file(b, params)
    np.savez(out, d=np.asarray(dm.map_sky[0], dtype=np.float64), ivar=np.asarray(dm.map_inv_var[0], dtype=np.float64),
             freq=float(b.freq), fwhm=float(dm.fwhm), d_dtype=str(dm.map_sky.dtype), ivar_dtype=str(dm.map_inv_var.dtype))
    print(f"{name}: freq {float(b.freq)} fwhm {dm.fwhm} dtype {dm.map_sky.dtype}/{dm.map_inv_var.dtype} npix {dm.map_sky.shape} "
          f"{time.time()-t1:.1f}s", flush=True)
group = params.compsep.mcmc_sampling_groups.sample_spectral_indices
mfile = os.path.join(cache, "mask.npy")
if not os.path.exists(mfile):
    keep = resolve_chisq_mask(hp.read_map(group.chisq_mask, field=0, dtype=np.float64), 64)
    np.save(mfile, keep)
    print("mask kept pixels", int(keep.sum()))
print(f"bands {len(names)}, still to do {todo}, elapsed {time.time()-t0:.1f}s")
