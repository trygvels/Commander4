# Samplers. Every sampler returns (draws::Matrix n x 5 [a_cmb A_s A_d beta_s beta_d], cost::NamedTuple)
# cost fields: nsolve = exact amplitude solves, neval = likelihood evaluations,
#              ngrad = gradient evaluations (leapfrog steps), all including burn-in/warmup.

using AdvancedHMC, LogDensityProblems, LogDensityProblemsAD, ForwardDiff

# ------------------------------------------------------------------
# 1) Commander3-style Gibbs: exact amplitude draw, then per-index random-walk MH.
#    Proposal length tuned during burn-in with Commander's rule
#    (acceptance > 0.7 -> x1.5, < 0.3 -> x0.5), k proposals per index per Gibbs iteration.
#    lnl = :chisq (beta | amplitudes) or :marginal (marginalize own amplitude, then resample amplitudes).
# ------------------------------------------------------------------
function gibbs_mh(rng, p::Pixel; niter = 4000, nburn = 1000, k = 10, lnl = :chisq)
    bs, bd = PRIOR_S[1], PRIOR_D[1]
    a = draw_amps(rng, p, bs, bd)
    nsolve, neval = 1, 0
    step = MVector(0.05, 0.02)
    nacc = MVector(0, 0); nprop = MVector(0, 0)
    out = Matrix{Float64}(undef, niter, 5)
    for it in 1:(nburn + niter)
        if lnl === :chisq
            a = draw_amps(rng, p, bs, bd); nsolve += 1
        end
        for j in 1:2
            lp(x, y) = lnl === :chisq ? logp_cond(p, a, x, y) : logp_marg1(p, a, x, y, j + 1)
            cur = lp(bs, bd); neval += 1
            for _ in 1:k
                bs2 = j == 1 ? bs + step[1] * randn(rng) : bs
                bd2 = j == 2 ? bd + step[2] * randn(rng) : bd
                new = lp(bs2, bd2); neval += 1
                if log(rand(rng)) < new - cur
                    bs, bd, cur = bs2, bd2, new
                    nacc[j] += 1
                end
                nprop[j] += 1
            end
            if lnl === :marginal
                a = draw_amps(rng, p, bs, bd); nsolve += 1   # keeps the partially collapsed Gibbs sampler valid
            end
            if it <= nburn && nprop[j] >= 100
                acc = nacc[j] / nprop[j]
                step[j] *= acc > 0.7 ? 1.5 : (acc < 0.3 ? 0.5 : 1.0)
                nacc[j] = 0; nprop[j] = 0
            end
        end
        if it > nburn
            out[it - nburn, :] .= (a[1], a[2], a[3], bs, bd)
        end
    end
    out, (nsolve = nsolve, neval = neval, ngrad = 0)
end

# ------------------------------------------------------------------
# LogDensityProblems wrappers (ForwardDiff gradients)
# ------------------------------------------------------------------
mutable struct CondTarget
    p::Pixel
    a::SVector{3,Float64}
end
LogDensityProblems.logdensity(t::CondTarget, x) = logp_cond(t.p, t.a, x[1], x[2])
LogDensityProblems.dimension(::CondTarget) = 2
LogDensityProblems.capabilities(::Type{CondTarget}) = LogDensityProblems.LogDensityOrder{0}()

struct JointTarget; p::Pixel; end
LogDensityProblems.logdensity(t::JointTarget, x) = logp_joint(t.p, x)
LogDensityProblems.dimension(::JointTarget) = 5
LogDensityProblems.capabilities(::Type{JointTarget}) = LogDensityProblems.LogDensityOrder{0}()

struct CollapsedTarget; p::Pixel; end
LogDensityProblems.logdensity(t::CollapsedTarget, x) = logp_collapsed(t.p, x[1], x[2])
LogDensityProblems.dimension(::CollapsedTarget) = 2
LogDensityProblems.capabilities(::Type{CollapsedTarget}) = LogDensityProblems.LogDensityOrder{0}()

function run_nuts(rng, target, x0; nadapt = 1000, nsamp = 2000)
    ell = ADgradient(:ForwardDiff, target)
    D = LogDensityProblems.dimension(target)
    metric = DiagEuclideanMetric(D)
    ham = Hamiltonian(metric, ell)
    eps0 = find_good_stepsize(rng, ham, x0)
    integ = Leapfrog(eps0)
    kernel = HMCKernel(Trajectory{MultinomialTS}(integ, GeneralisedNoUTurn()))
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integ))
    samples, stats = sample(rng, ham, kernel, x0, nsamp + nadapt, adaptor, nadapt;
                            progress = false, verbose = false, drop_warmup = false)
    ngrad = sum(s.n_steps for s in stats)
    samples[(nadapt + 1):end], ngrad, stats
end

# ------------------------------------------------------------------
# 2) Gibbs with NUTS replacing MH: exact amplitude draw, then one NUTS transition on (beta_s, beta_d) | a.
#    Step size and diagonal metric from a short pilot adaptation at the initial amplitudes.
# ------------------------------------------------------------------
function gibbs_nuts(rng, p::Pixel; niter = 4000, nburn = 1000)
    bs, bd = PRIOR_S[1], PRIOR_D[1]
    a = draw_amps(rng, p, bs, bd)
    tgt = CondTarget(p, a)
    ell = ADgradient(:ForwardDiff, tgt)
    metric = DiagEuclideanMetric(2)
    ham = Hamiltonian(metric, ell)
    x = [bs, bd]
    eps0 = find_good_stepsize(rng, ham, x)
    integ = Leapfrog(eps0)
    kernel = HMCKernel(Trajectory{MultinomialTS}(integ, GeneralisedNoUTurn()))
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integ))
    ham, _ = AdvancedHMC.sample_init(rng, ham, x)
    nsolve, ngrad = 1, 0
    out = Matrix{Float64}(undef, niter, 5)
    for it in 1:(nburn + niter)
        a = draw_amps(rng, p, x[1], x[2]); nsolve += 1
        tgt.a = a                                   # conditional target changes every Gibbs iteration
        z = AdvancedHMC.phasepoint(rng, x, ham); ngrad += 1
        tr = AdvancedHMC.transition(rng, ham, kernel, z)
        st = AdvancedHMC.stat(tr)
        ngrad += st.n_steps
        # adapt step size and diagonal metric during burn-in (Stan windows), frozen afterwards
        ham, kernel, _ = AdvancedHMC.adapt!(ham, kernel, adaptor, it, nburn, tr.z, st.acceptance_rate)
        x = tr.z.θ
        if it > nburn
            out[it - nburn, :] .= (a[1], a[2], a[3], x[1], x[2])
        end
    end
    out, (nsolve = nsolve, neval = 0, ngrad = ngrad)
end

# ------------------------------------------------------------------
# 3) Joint NUTS over amplitudes and indices: no Gibbs alternation, no amplitude solves.
# ------------------------------------------------------------------
function joint_nuts(rng, p::Pixel; nadapt = 1000, nsamp = 2000)
    a0 = draw_amps(rng, p, PRIOR_S[1], PRIOR_D[1])
    x0 = [a0[1], a0[2], a0[3], PRIOR_S[1], PRIOR_D[1]]
    s, ngrad, _ = run_nuts(rng, JointTarget(p), x0; nadapt, nsamp)
    out = reduce(vcat, permutedims.(s))
    out, (nsolve = 1, neval = 0, ngrad = ngrad)
end

# ------------------------------------------------------------------
# 4) Collapsed NUTS: amplitudes marginalized analytically, NUTS on (beta_s, beta_d),
#    then one exact amplitude draw per retained sample.
# ------------------------------------------------------------------
function collapsed_nuts(rng, p::Pixel; nadapt = 1000, nsamp = 2000)
    s, ngrad, _ = run_nuts(rng, CollapsedTarget(p), [PRIOR_S[1], PRIOR_D[1]]; nadapt, nsamp)
    out = Matrix{Float64}(undef, nsamp, 5)
    for (i, x) in enumerate(s)
        a = draw_amps(rng, p, x[1], x[2])
        out[i, :] .= (a[1], a[2], a[3], x[1], x[2])
    end
    out, (nsolve = nsamp, neval = 0, ngrad = ngrad)
end
