# Toy per-pixel foreground model mirroring Commander3's intensity component separation:
#   d_nu = g(nu) * a_cmb + A_s (nu/30)^beta_s + A_d MBB(nu; beta_d, T_d=20K) + n_nu
# Units: uK_RJ. Amplitudes have flat priors, spectral indices Gaussian priors.

using StaticArrays, LinearAlgebra, Random

const H_OVER_K = 0.04799243      # K per GHz
const TCMB     = 2.7255          # K
const TDUST    = 20.0            # K, fixed
const NB       = 8
const FREQS    = SVector{NB,Float64}(23.0, 30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0)
# Toy white-noise rms per nside-16 pixel, uK_RJ (WMAP K + Planck-like bands)
const SIGMA0   = SVector{NB,Float64}(4.0, 3.0, 3.0, 3.5, 1.2, 0.8, 1.2, 4.0)

const PRIOR_S  = (-3.0, 0.2)     # beta_s ~ N(mean, sd)
const PRIOR_D  = (1.6, 0.1)      # beta_d ~ N(mean, sd)

cmb_rj(nu)            = (x = H_OVER_K * nu / TCMB; x^2 * exp(x) / expm1(x)^2)
synch_sed(nu, beta)   = (nu / 30.0)^beta
dust_sed(nu, beta)    = (nu / 353.0)^(beta + 1) * expm1(H_OVER_K * 353.0 / TDUST) / expm1(H_OVER_K * nu / TDUST)

const GCMB = SVector{NB,Float64}(cmb_rj.(FREQS))

# 8x3 mixing matrix: columns CMB, synchrotron, dust
@inline mixing(bs, bd) = hcat(GCMB, synch_sed.(FREQS, bs), dust_sed.(FREQS, bd))

@inline logprior(bs, bd) = -0.5 * ((bs - PRIOR_S[1]) / PRIOR_S[2])^2 - 0.5 * ((bd - PRIOR_D[1]) / PRIOR_D[2])^2

struct Pixel
    d::SVector{NB,Float64}       # data
    w::SVector{NB,Float64}       # inverse noise variance
    truth::SVector{5,Float64}    # a_cmb, A_s, A_d, beta_s, beta_d
    t::Float64                   # 0 = high latitude (faint), 1 = plane (bright)
end

"""
Simulate `npix` independent pixels. Foreground amplitudes rise with t in [0,1],
so the set spans the full range from noise-dominated to bright pixels.
Spectral indices are drawn from their priors so coverage is testable.
"""
function simulate(npix; seed = 1, noise_scale = 1.0)
    rng = Xoshiro(seed)
    sig = SIGMA0 .* noise_scale
    w = 1.0 ./ sig .^ 2
    pix = Vector{Pixel}(undef, npix)
    for i in 1:npix
        t   = (i - 0.5) / npix
        As  = 10^(0.0 + 2.5t + 0.1randn(rng))    # ~1 to ~300 uK_RJ at 30 GHz
        Ad  = 10^(0.5 + 2.5t + 0.1randn(rng))    # ~3 to ~1000 uK_RJ at 353 GHz
        cmb = 70.0randn(rng)
        bs  = PRIOR_S[1] + PRIOR_S[2] * randn(rng)
        bd  = PRIOR_D[1] + PRIOR_D[2] * randn(rng)
        a   = SVector(cmb, As, Ad)
        d   = mixing(bs, bd) * a + sig .* randn(rng, SVector{NB,Float64})
        pix[i] = Pixel(d, w, SVector(cmb, As, Ad, bs, bd), t)
    end
    pix
end

# ---------- log densities ----------

@inline function loglike(p::Pixel, a, bs, bd)
    r = p.d - mixing(bs, bd) * a
    -0.5 * sum(p.w .* r .^ 2)
end

# Commander 'chisq' type: beta | amplitudes
@inline logp_cond(p::Pixel, a, bs, bd) = loglike(p, a, bs, bd) + logprior(bs, bd)

# Joint over theta = (a_cmb, A_s, A_d, beta_s, beta_d), flat amplitude priors
@inline logp_joint(p::Pixel, th) = loglike(p, SVector(th[1], th[2], th[3]), th[4], th[5]) + logprior(th[4], th[5])

# Commander 'marginal' type: marginalize the amplitude of component j (2 = synch, 3 = dust)
# conditional on the other amplitudes. Same form as comp_lnL_max_chisq_diagonal with use_det.
@inline function logp_marg1(p::Pixel, a, bs, bd, j)
    A = mixing(bs, bd)
    m = A[:, j]
    r = p.d - A * a + m * a[j]
    F = sum(p.w .* m .^ 2)
    g = sum(p.w .* m .* r)
    0.5 * g^2 / F - 0.5 * log(F) + logprior(bs, bd)
end

# Fully collapsed: all three amplitudes marginalized analytically
@inline function logp_collapsed(p::Pixel, bs, bd)
    A  = mixing(bs, bd)
    WA = p.w .* A
    M  = A' * WA
    b  = WA' * p.d
    dM = det(M)
    dM > 0 || return convert(typeof(dM), -Inf)   # numerically singular far out in the tails
    0.5 * dot(b, M \ b) - 0.5 * log(dM) + logprior(bs, bd)
end

# Exact amplitude draw a | beta, d  (the 'amplitude solve'; in Commander3 this is the global CG step)
function draw_amps(rng, p::Pixel, bs, bd)
    A  = mixing(bs, bd)
    WA = p.w .* A
    M  = Symmetric(A' * WA)
    b  = WA' * p.d
    C  = cholesky(M)
    mu = C \ b
    mu + C.U \ randn(rng, SVector{3,Float64})
end
