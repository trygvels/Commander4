import numpy as np, time
npix=12*512**2; nb=8
f=np.array([23,30,44,70,100,143,217,353.]); hk=0.04799243
x=hk*f/2.7255; g=x**2*np.exp(x)/np.expm1(x)**2
def mix(bs,bd): return np.stack([g,(f/30)**bs,(f/353)**(bd+1)*np.expm1(hk*353/20)/np.expm1(hk*f/20)],1)  # (8,3)
rng=np.random.default_rng(1); w=1/np.array([4,3,3,3.5,1.2,.8,1.2,4.])**2
a=np.stack([70*rng.standard_normal(npix),10**rng.uniform(0,2.5,npix),10**rng.uniform(.5,3,npix)],1)
d=a@mix(-3.05,1.58).T+rng.standard_normal((npix,nb))/np.sqrt(w)
W=np.broadcast_to(w,(npix,nb)).copy()   # per-pixel inverse variance, as in real maps
def marg(bs,bd):
    A=mix(bs,bd)                                     # (8,3)
    M=np.einsum('pb,bi,bj->pij',W,A,A)               # (npix,3,3)
    b=np.einsum('pb,bi->pi',W*d,A)
    sol=np.linalg.solve(M,b[...,None])[...,0]
    sign,ld=np.linalg.slogdet(M)
    return 0.5*np.sum(b*sol)-0.5*np.sum(ld)
marg(-3.0,1.6)
t=time.perf_counter(); n=5
for i in range(n): v=marg(-3.0+0.01*i,1.6)
print(f"numpy: {(time.perf_counter()-t)/n:.3f} s per full-sky marginal lnL (npix={npix})")

# --- optimized NumPy: closed-form 3x3 symmetric solve/det ---
def marg_fast(bs,bd):
    A=mix(bs,bd)
    WA=[W*A[:,i] for i in range(3)]
    m00=WA[0]@np.ones(nb)*0; # placeholder to keep shapes simple
    m=lambda i,j: (W*(A[:,i]*A[:,j])).sum(1)
    M00,M01,M02,M11,M12,M22=m(0,0),m(0,1),m(0,2),m(1,1),m(1,2),m(2,2)
    Wd=W*d
    b0,b1,b2=Wd@A[:,0],Wd@A[:,1],Wd@A[:,2]
    c00=M11*M22-M12*M12; c01=M02*M12-M01*M22; c02=M01*M12-M02*M11
    c11=M00*M22-M02*M02; c12=M01*M02-M00*M12; c22=M00*M11-M01*M01
    det=M00*c00+M01*c01+M02*c02
    q=(b0*(c00*b0+c01*b1+c02*b2)+b1*(c01*b0+c11*b1+c12*b2)+b2*(c02*b0+c12*b1+c22*b2))/det
    return 0.5*q.sum()-0.5*np.log(det).sum()
print("check", marg(-3.0,1.6), marg_fast(-3.0,1.6))
t=time.perf_counter()
for i in range(n): v=marg_fast(-3.0+0.01*i,1.6)
print(f"numpy closed-form: {(time.perf_counter()-t)/n:.3f} s")

import numba as nb_
@nb_.njit(parallel=True, fastmath=True)
def marg_nb(A, W, d):
    npix=d.shape[0]; out=np.zeros(npix)
    for p in nb_.prange(npix):
        M00=M01=M02=M11=M12=M22=b0=b1=b2=0.0
        for k in range(A.shape[0]):
            wk=W[p,k]; x0=A[k,0]; x1=A[k,1]; x2=A[k,2]; wd=wk*d[p,k]
            M00+=wk*x0*x0; M01+=wk*x0*x1; M02+=wk*x0*x2; M11+=wk*x1*x1; M12+=wk*x1*x2; M22+=wk*x2*x2
            b0+=wd*x0; b1+=wd*x1; b2+=wd*x2
        c00=M11*M22-M12*M12; c01=M02*M12-M01*M22; c02=M01*M12-M02*M11
        c11=M00*M22-M02*M02; c12=M01*M02-M00*M12; c22=M00*M11-M01*M01
        det=M00*c00+M01*c01+M02*c02
        q=(b0*(c00*b0+c01*b1+c02*b2)+b1*(c01*b0+c11*b1+c12*b2)+b2*(c02*b0+c12*b1+c22*b2))/det
        out[p]=0.5*q-0.5*np.log(det)
    return out.sum()
A0=mix(-3.0,1.6); print("numba check", marg_nb(A0,W,d))
t=time.perf_counter()
for i in range(n): v=marg_nb(mix(-3.0+0.01*i,1.6),W,d)
print(f"numba (2 threads): {(time.perf_counter()-t)/n:.3f} s")
