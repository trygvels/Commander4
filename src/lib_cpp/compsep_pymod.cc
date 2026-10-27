#include "ducc0/bindings/pybind_utils.h"  // must be the first include

#include <algorithm>
#include <cmath>
#include <vector>

#include "ducc0/infra/threading.h"

namespace cmdr4 {

namespace detail_pymodule_compsep {

using namespace std;
using namespace ducc0;

/** Solves A x = b for a small symmetric positive definite A by Cholesky decomposition A = L L^T.
 *
 * A is overwritten: its lower triangle holds L afterwards. Returns false, leaving x unset, if A is
 * not positive definite.
 */
bool solve_cholesky(vmav<double,2> &A, const vmav<double,1> &b, vmav<double,1> &x){
    const size_t n = A.shape(0);
    // Cholesky decomposition, column by column. Entry (i,j) of L only needs the entries of L left
    // of column j, which are already stored in A, so A can be overwritten as we go.
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j <= i; j++){
            double sum = 0.0;
            for (size_t k = 0; k < j; k++)
                sum += A(i, k) * A(j, k);
            if (i == j){
                const double val = A(i, i) - sum;
                if (val <= 0.0) return false;  // Not positive definite.
                A(i, i) = std::sqrt(val);
            } else {
                A(i, j) = (1.0 / A(j, j)) * (A(i, j) - sum);
            }
        }
    }
    // Forward substitution L y = b, storing y in x.
    for (size_t i = 0; i < n; i++){
        double sum = 0.0;
        for (size_t j = 0; j < i; j++)
            sum += A(i, j) * x(j);
        x(i) = (b(i) - sum) / A(i, i);
    }
    // Backward substitution L^T x = y, in place, for i = n-1, ..., 0.
    for (size_t i = n; i-- > 0;){
        double sum = 0.0;
        for (size_t j = i + 1; j < n; j++)
            sum += A(j, i) * x(j);
        x(i) = (x(i) - sum) / A(i, i);
    }
    return true;
}


/** Solves L L^T x = b for x, given the Cholesky factor L in the lower triangle of `L` (as left by
 * solve_cholesky). */
void cholesky_solve(const vmav<double,2> &L, const vmav<double,1> &b, vmav<double,1> &x){
    const size_t n = L.shape(0);
    for (size_t i = 0; i < n; i++){
        double sum = 0.0;
        for (size_t j = 0; j < i; j++)
            sum += L(i, j) * x(j);
        x(i) = (b(i) - sum) / L(i, i);
    }
    for (size_t i = n; i-- > 0;){
        double sum = 0.0;
        for (size_t j = i + 1; j < n; j++)
            sum += L(j, i) * x(j);
        x(i) = (x(i) - sum) / L(i, i);
    }
}


/** Shared body of the per-pixel amplitude draw (Py_solve_perpix and Py_solve_perpix_regions).
 *
 * Solves (M_p^T N^-1 M_p) x = M_p^T (N^-1 d + N^-1/2 eta) per pixel, where N^-1 is diagonal over
 * the bands and holds the inverse noise variance, d is the sky map, eta are standard normal random
 * numbers and M_p is the (nband, ncomp) mixing matrix of pixel p, given by `region_of(p)` (an index
 * r >= 0, or r < 0 for a pixel to leave at zero) and `M(r, b, c)`. The arithmetic is the same
 * for both callers, so a single region reproduces the scalar kernel bit for bit. Pixels where
 * M_p^T N^-1 M_p is not positive definite (e.g. zero inverse variance in every band) get zero
 * amplitudes.
 */
template<typename RegionOf, typename Mixing>
void solve_perpix_body(const cmav<double,2> &map_sky, const cmav<double,2> &map_inv_var,
                       const cmav<double,2> &rand, size_t ncomp, size_t nthreads,
                       vmav<double,2> &comp_maps, RegionOf &&region_of, Mixing &&M){
    const size_t nband = map_sky.shape(0), npix = map_sky.shape(1);
    execStatic(npix, nthreads, 0, [&](Scheduler &sched){
        // Per-thread workspace, reused for every pixel.
        vmav<double,2> A({ncomp, ncomp});
        vmav<double,1> rhs({ncomp}), x({ncomp});
        while (auto rng = sched.getNext())
            for (size_t ipix = rng.lo; ipix < rng.hi; ipix++){
                const int64_t r = region_of(ipix);
                if (r < 0){
                    for (size_t c = 0; c < ncomp; c++)
                        comp_maps(c, ipix) = 0.0;
                    continue;
                }
                // A = M^T N^-1 M and rhs = M^T (N^-1 d + N^-1/2 eta).
                for (size_t c1 = 0; c1 < ncomp; c1++){
                    rhs(c1) = 0.0;
                    for (size_t c2 = 0; c2 < ncomp; c2++)
                        A(c1, c2) = 0.0;
                }
                for (size_t b = 0; b < nband; b++){
                    const double inv_var = map_inv_var(b, ipix);
                    const double weighted_d = inv_var * map_sky(b, ipix)
                                            + rand(ipix, b) * std::sqrt(inv_var);
                    for (size_t c1 = 0; c1 < ncomp; c1++){
                        rhs(c1) += M(r, b, c1) * weighted_d;
                        for (size_t c2 = 0; c2 < ncomp; c2++)
                            A(c1, c2) += M(r, b, c1) * inv_var * M(r, b, c2);
                    }
                }
                const bool success = solve_cholesky(A, rhs, x);
                for (size_t c = 0; c < ncomp; c++)
                    comp_maps(c, ipix) = success ? x(c) : 0.0;
            }
    });
}


/** Draws one sample of the component amplitudes per pixel, independently in each pixel.
 *
 * Solves (M^T N^-1 M) x = M^T (N^-1 d + N^-1/2 eta) per pixel, where N^-1 is diagonal over the
 * bands and holds the inverse noise variance, d is the sky map and eta are standard normal random
 * numbers. Pixels where M^T N^-1 M is not positive definite (e.g. zero inverse variance in every
 * band) get zero amplitudes.
 *
 * Args:
 *   map_sky -- (nband, npix) sky maps.
 *   map_inv_var -- (nband, npix) inverse noise variance maps (0 where unobserved).
 *   M -- (nband, ncomp) mixing matrix.
 *   rand -- (npix, nband) standard normal random numbers (pixel-major, so that the numbers of
 *           one pixel are adjacent in memory).
 *   nthreads -- number of threads.
 * Returns:
 *   (ncomp, npix) component amplitude maps.
 */
NpArr Py_solve_perpix(const CNpArr &map_sky_, const CNpArr &map_inv_var_, const CNpArr &M_,
                      const CNpArr &rand_, size_t nthreads){
    auto map_sky = to_cmav<double,2>(map_sky_, "map_sky");
    auto map_inv_var = to_cmav<double,2>(map_inv_var_, "map_inv_var");
    auto M = to_cmav<double,2>(M_, "M");
    auto rand = to_cmav<double,2>(rand_, "rand");
    const size_t nband = map_sky.shape(0), npix = map_sky.shape(1), ncomp = M.shape(1);
    MR_assert(map_inv_var.shape(0)==nband && map_inv_var.shape(1)==npix,
              "map_inv_var and map_sky shapes differ");
    MR_assert(M.shape(0)==nband, "M must have shape (nband, ncomp)");
    MR_assert(rand.shape(0)==npix && rand.shape(1)==nband, "rand must have shape (npix, nband)");

    auto comp_maps_ = make_Pyarr<double>({ncomp, npix});
    auto comp_maps = to_vmav<double,2>(comp_maps_);
    solve_perpix_body(map_sky, map_inv_var, rand, ncomp, nthreads, comp_maps,
                      [](size_t) -> int64_t { return 0; },
                      [&M](int64_t, size_t b, size_t c) { return M(b, c); });
    return comp_maps_;
}


/** Per-pixel amplitude draw with a mixing matrix that is constant within regions.
 *
 * As Py_solve_perpix, with M_reg[region[p]] as the mixing matrix of pixel p: the amplitude draw
 * conditional on spatially varying spectral indices (one index value per region). Pixels with
 * region < 0 get zero amplitudes. With a single region (region = 0 everywhere, M_reg of shape
 * (1, nband, ncomp)) the result equals Py_solve_perpix bit for bit, since the two share
 * solve_perpix_body. Prior terms enter through pseudo-bands appended to the stream (see
 * compsep/index_field.py): a Gaussian prior a_c ~ N(m, t^2) is a band with map m, inverse
 * variance 1/t^2 and mixing row e_c, so the draw is N(A'^-1 b', A'^-1) with A' = A + T^-1 and
 * b' = b + T^-1 m, the fluctuation of the pseudo-band being T^-1/2 eta.
 *
 * Args:
 *   map_sky -- (nband, npix) sky maps.
 *   map_inv_var -- (nband, npix) inverse noise variance maps (0 where unobserved).
 *   region -- (npix,) int64 region index per pixel, -1 to leave the pixel at zero.
 *   M_reg -- (nreg, nband, ncomp) mixing matrix per region.
 *   rand -- (npix, nband) standard normal random numbers (pixel-major).
 *   nthreads -- number of threads.
 * Returns:
 *   (ncomp, npix) component amplitude maps.
 */
NpArr Py_solve_perpix_regions(const CNpArr &map_sky_, const CNpArr &map_inv_var_,
                              const CNpArr &region_, const CNpArr &M_reg_, const CNpArr &rand_,
                              size_t nthreads){
    auto map_sky = to_cmav<double,2>(map_sky_, "map_sky");
    auto map_inv_var = to_cmav<double,2>(map_inv_var_, "map_inv_var");
    auto region = to_cmav<int64_t,1>(region_, "region");
    auto M_reg = to_cmav<double,3>(M_reg_, "M_reg");
    auto rand = to_cmav<double,2>(rand_, "rand");
    const size_t nband = map_sky.shape(0), npix = map_sky.shape(1);
    const size_t nreg = M_reg.shape(0), ncomp = M_reg.shape(2);
    MR_assert(map_inv_var.shape(0)==nband && map_inv_var.shape(1)==npix,
              "map_inv_var and map_sky shapes differ");
    MR_assert(region.shape(0)==npix, "region must have one entry per pixel");
    MR_assert(M_reg.shape(1)==nband, "M_reg must have shape (nreg, nband, ncomp)");
    MR_assert(rand.shape(0)==npix && rand.shape(1)==nband, "rand must have shape (npix, nband)");
    for (size_t ipix = 0; ipix < npix; ipix++)
        MR_assert(region(ipix) < int64_t(nreg), "region index out of range");

    auto comp_maps_ = make_Pyarr<double>({ncomp, npix});
    auto comp_maps = to_vmav<double,2>(comp_maps_);
    solve_perpix_body(map_sky, map_inv_var, rand, ncomp, nthreads, comp_maps,
                      [&region](size_t ipix) -> int64_t { return region(ipix); },
                      [&M_reg](int64_t r, size_t b, size_t c) { return M_reg(r, b, c); });
    return comp_maps_;
}


/** Log marginal likelihood of the mixing matrix, with the per-pixel amplitudes integrated out.
 *
 * For the per-pixel data model d = M a + n with diagonal noise and a flat amplitude prior, the
 * amplitudes can be marginalized analytically in each pixel:
 *
 *     ln p(d | M) = sum_pix [ 1/2 b^T A^-1 b - 1/2 ln det A ] + const,
 *     A = M^T N^-1 M,   b = M^T N^-1 d,
 *
 * where the constant (-1/2 d^T N^-1 d and the 2 pi factors) does not depend on M. This is the
 * acceptance quantity for a Metropolis step on parameters of M (e.g. spectral indices) when the
 * amplitudes are drawn from their conditional afterwards; it is the same system Py_solve_perpix
 * solves, without the fluctuation term.
 *
 * With use_logdet false the ln det term is dropped, giving the "ridge" (profile) likelihood
 * 1/2 b^T A^-1 b = max_a ln p(d | a, M) + const. As an MH target this is the marginal above times
 * det(A)^(1/2), i.e. the flat-amplitude-prior posterior with an extra Jeffreys-type factor on the
 * parameters of M; it removes the prior-volume term that a flat amplitude prior, defined at an
 * SED reference frequency, imposes on those parameters (C3's 'ridge' lnL type).
 *
 * A pixel enters the sum only if at least ncomp bands observe it (inverse variance > 0) and A is
 * positive definite. The first condition depends only on the inverse-variance maps, never on M,
 * so the set of pixels does not change between the two states of an MH step. Pixels that fail the
 * second condition are counted in the returned `nskipped`; callers should treat a nonzero count
 * for an observed pixel as a numerical problem rather than silently ignoring it.
 *
 * Args:
 *   map_sky -- (nband, npix) sky (residual) maps.
 *   map_inv_var -- (nband, npix) inverse noise variance maps (0 where unobserved or masked).
 *   M -- (nband, ncomp) mixing matrix.
 *   nthreads -- number of threads.
 *   use_logdet -- include the -1/2 ln det A term (marginal) or not (ridge).
 * Returns:
 *   (loglike, nused, nskipped): the summed log marginal likelihood, the number of pixels that
 *   entered it, and the number of observed pixels skipped because A was not positive definite.
 */
py::tuple Py_marginal_loglike_perpix(const CNpArr &map_sky_, const CNpArr &map_inv_var_,
                                     const CNpArr &M_, size_t nthreads, bool use_logdet){
    auto map_sky = to_cmav<double,2>(map_sky_, "map_sky");
    auto map_inv_var = to_cmav<double,2>(map_inv_var_, "map_inv_var");
    auto M = to_cmav<double,2>(M_, "M");
    const size_t nband = map_sky.shape(0), npix = map_sky.shape(1), ncomp = M.shape(1);
    MR_assert(map_inv_var.shape(0)==nband && map_inv_var.shape(1)==npix,
              "map_inv_var and map_sky shapes differ");
    MR_assert(M.shape(0)==nband, "M must have shape (nband, ncomp)");

    const size_t nthr = std::max<size_t>(1, nthreads);
    std::vector<double> partial(nthr, 0.0);
    std::vector<size_t> used(nthr, 0), skipped(nthr, 0);
    execStatic(npix, nthr, 0, [&](Scheduler &sched){
        const size_t ithr = sched.thread_num();
        vmav<double,2> A({ncomp, ncomp});
        vmav<double,1> rhs({ncomp}), x({ncomp});
        double sum = 0.0;
        size_t nused = 0, nskip = 0;
        while (auto rng = sched.getNext())
            for (size_t ipix = rng.lo; ipix < rng.hi; ipix++){
                size_t nobs = 0;
                for (size_t b = 0; b < nband; b++)
                    if (map_inv_var(b, ipix) > 0.0) nobs++;
                if (nobs < ncomp) continue;  // Underdetermined for every M: never in the sum.
                for (size_t c1 = 0; c1 < ncomp; c1++){
                    rhs(c1) = 0.0;
                    for (size_t c2 = 0; c2 < ncomp; c2++)
                        A(c1, c2) = 0.0;
                }
                for (size_t b = 0; b < nband; b++){
                    const double inv_var = map_inv_var(b, ipix);
                    if (inv_var <= 0.0) continue;
                    const double weighted_d = inv_var * map_sky(b, ipix);
                    for (size_t c1 = 0; c1 < ncomp; c1++){
                        rhs(c1) += M(b, c1) * weighted_d;
                        for (size_t c2 = 0; c2 <= c1; c2++)
                            A(c1, c2) += M(b, c1) * inv_var * M(b, c2);
                    }
                }
                for (size_t c1 = 0; c1 < ncomp; c1++)
                    for (size_t c2 = c1 + 1; c2 < ncomp; c2++)
                        A(c1, c2) = A(c2, c1);
                if (!solve_cholesky(A, rhs, x)){
                    nskip++;
                    continue;
                }
                // A now holds L (A = L L^T), x = A^-1 b. ln det A = 2 sum ln L_ii.
                double quad = 0.0, logdet = 0.0;
                for (size_t c = 0; c < ncomp; c++){
                    quad += rhs(c) * x(c);
                    logdet += 2.0 * std::log(A(c, c));
                }
                sum += use_logdet ? 0.5 * quad - 0.5 * logdet : 0.5 * quad;
                nused++;
            }
        partial[ithr] = sum;
        used[ithr] = nused;
        skipped[ithr] = nskip;
    });
    double total = 0.0;
    size_t nused = 0, nskip = 0;
    for (size_t i = 0; i < nthr; i++){
        total += partial[i];
        nused += used[i];
        nskip += skipped[i];
    }
    return py::make_tuple(total, nused, nskip);
}


/** Per-region marginal (or ridge) log likelihood and its gradient, for a mixing matrix that is
 * constant within each region but differs between regions (spatially varying spectral indices).
 *
 * Pixel p belongs to region r = region[p] (pixels with r < 0 are skipped) and uses mixing matrix
 * M_reg[r]. The per-pixel term is the one of marginal_loglike_perpix. Each component c carries one
 * parameter per region, theta_rc, with d M[:, c] / d theta_rc = dM_reg[r, :, c] and no other
 * column depending on it (for power-law and modified-blackbody SEDs, dM = M * ln(nu / nu_ref)).
 * With x = A^-1 b and g = dM[:, c], the gradient of one pixel's term is
 *
 *     d/dtheta [ b^T A^-1 b / 2 ] = x_c g^T N^-1 (d - M x),
 *     d/dtheta [ -ln det A / 2 ]  = -(A^-1 M^T N^-1 g)_c,
 *
 * summed over the region's pixels. Components without a parameter should pass zero columns in
 * dM_reg; their gradient entries are then zero.
 *
 * Args:
 *   map_sky, map_inv_var -- (nband, npix) maps.
 *   region -- (npix,) int64 region index per pixel, -1 to exclude a pixel.
 *   M_reg, dM_reg -- (nreg, nband, ncomp) mixing matrices and their parameter derivatives.
 *   nthreads -- number of threads.
 *   use_logdet -- marginal (true) or ridge (false) likelihood.
 *   want_grad -- whether to compute the gradient (skipped otherwise).
 *   want_extras -- also return the per-pixel conditional mean x = A^-1 b (ncomp, npix), the
 *                  diagonal of A^-1 (ncomp, npix) and the used-pixel mask (npix,), zeros / false
 *                  for pixels that did not enter the sum (region < 0, fewer than ncomp observed
 *                  bands, or A not positive definite). The trailing terms of a hierarchical
 *                  amplitude prior must be summed over exactly this pixel set.
 * Returns:
 *   (loglike (nreg,), grad (nreg, ncomp), nskipped): per-region log likelihood, its gradient, and
 *   the number of observed, well-determined pixels skipped because A was not positive definite;
 *   with want_extras, (loglike, grad, nskipped, x, diag_Ainv, used).
 */
py::tuple Py_marginal_loglike_regions(const CNpArr &map_sky_, const CNpArr &map_inv_var_,
                                      const CNpArr &region_, const CNpArr &M_reg_,
                                      const CNpArr &dM_reg_, size_t nthreads, bool use_logdet,
                                      bool want_grad, bool want_extras){
    auto map_sky = to_cmav<double,2>(map_sky_, "map_sky");
    auto map_inv_var = to_cmav<double,2>(map_inv_var_, "map_inv_var");
    auto region = to_cmav<int64_t,1>(region_, "region");
    auto M_reg = to_cmav<double,3>(M_reg_, "M_reg");
    auto dM_reg = to_cmav<double,3>(dM_reg_, "dM_reg");
    const size_t nband = map_sky.shape(0), npix = map_sky.shape(1);
    const size_t nreg = M_reg.shape(0), ncomp = M_reg.shape(2);
    MR_assert(map_inv_var.shape(0)==nband && map_inv_var.shape(1)==npix,
              "map_inv_var and map_sky shapes differ");
    MR_assert(region.shape(0)==npix, "region must have one entry per pixel");
    MR_assert(M_reg.shape(1)==nband, "M_reg must have shape (nreg, nband, ncomp)");
    MR_assert(dM_reg.shape(0)==nreg && dM_reg.shape(1)==nband && dM_reg.shape(2)==ncomp,
              "dM_reg must have the shape of M_reg");

    const size_t nthr = std::max<size_t>(1, nthreads);
    // Per-thread accumulators, so threads never write to the same region entry.
    std::vector<double> ll_thr(nthr * nreg, 0.0), grad_thr(nthr * nreg * ncomp, 0.0);
    std::vector<size_t> skipped(nthr, 0);
    // Extras (allocated at full size only when requested; otherwise empty placeholders).
    const size_t npix_x = want_extras ? npix : 0;
    auto x_out_ = make_Pyarr<double>({ncomp, npix_x});
    auto diag_out_ = make_Pyarr<double>({ncomp, npix_x});
    auto used_out_ = make_Pyarr<bool>({npix_x});
    auto x_out = to_vmav<double,2>(x_out_);
    auto diag_out = to_vmav<double,2>(diag_out_);
    auto used_out = to_vmav<bool,1>(used_out_);
    if (want_extras)
        for (size_t ipix = 0; ipix < npix; ipix++){
            used_out(ipix) = false;
            for (size_t c = 0; c < ncomp; c++){
                x_out(c, ipix) = 0.0;
                diag_out(c, ipix) = 0.0;
            }
        }
    execStatic(npix, nthr, 0, [&](Scheduler &sched){
        const size_t ithr = sched.thread_num();
        double *ll = ll_thr.data() + ithr * nreg;
        double *gr = grad_thr.data() + ithr * nreg * ncomp;
        vmav<double,2> A({ncomp, ncomp});
        vmav<double,1> rhs({ncomp}), x({ncomp}), u({ncomp}), y({ncomp});
        std::vector<double> resid(nband);
        size_t nskip = 0;
        while (auto rng = sched.getNext())
            for (size_t ipix = rng.lo; ipix < rng.hi; ipix++){
                const int64_t r = region(ipix);
                if (r < 0) continue;
                MR_assert(size_t(r) < nreg, "region index out of range");
                size_t nobs = 0;
                for (size_t b = 0; b < nband; b++)
                    if (map_inv_var(b, ipix) > 0.0) nobs++;
                if (nobs < ncomp) continue;
                for (size_t c1 = 0; c1 < ncomp; c1++){
                    rhs(c1) = 0.0;
                    for (size_t c2 = 0; c2 < ncomp; c2++)
                        A(c1, c2) = 0.0;
                }
                for (size_t b = 0; b < nband; b++){
                    const double inv_var = map_inv_var(b, ipix);
                    if (inv_var <= 0.0) continue;
                    const double weighted_d = inv_var * map_sky(b, ipix);
                    for (size_t c1 = 0; c1 < ncomp; c1++){
                        rhs(c1) += M_reg(r, b, c1) * weighted_d;
                        for (size_t c2 = 0; c2 <= c1; c2++)
                            A(c1, c2) += M_reg(r, b, c1) * inv_var * M_reg(r, b, c2);
                    }
                }
                for (size_t c1 = 0; c1 < ncomp; c1++)
                    for (size_t c2 = c1 + 1; c2 < ncomp; c2++)
                        A(c1, c2) = A(c2, c1);
                if (!solve_cholesky(A, rhs, x)){
                    nskip++;
                    continue;
                }
                double quad = 0.0, logdet = 0.0;
                for (size_t c = 0; c < ncomp; c++){
                    quad += rhs(c) * x(c);
                    logdet += 2.0 * std::log(A(c, c));
                }
                ll[r] += use_logdet ? 0.5 * quad - 0.5 * logdet : 0.5 * quad;
                if (want_extras){
                    used_out(ipix) = true;
                    for (size_t c = 0; c < ncomp; c++){
                        x_out(c, ipix) = x(c);
                        // diag(A^-1)_c = e_c^T A^-1 e_c through the Cholesky factor.
                        for (size_t c1 = 0; c1 < ncomp; c1++) u(c1) = (c1 == c) ? 1.0 : 0.0;
                        cholesky_solve(A, u, y);
                        diag_out(c, ipix) = y(c);
                    }
                }
                if (!want_grad) continue;
                for (size_t b = 0; b < nband; b++){
                    double model = 0.0;
                    for (size_t c = 0; c < ncomp; c++)
                        model += M_reg(r, b, c) * x(c);
                    resid[b] = map_sky(b, ipix) - model;
                }
                for (size_t c = 0; c < ncomp; c++){
                    double g_dot_res = 0.0;
                    bool any = false;
                    for (size_t b = 0; b < nband; b++){
                        const double inv_var = map_inv_var(b, ipix);
                        if (inv_var <= 0.0 || dM_reg(r, b, c) == 0.0) continue;
                        any = true;
                        g_dot_res += dM_reg(r, b, c) * inv_var * resid[b];
                    }
                    if (!any) continue;
                    double grad = x(c) * g_dot_res;
                    if (use_logdet){
                        for (size_t c1 = 0; c1 < ncomp; c1++){
                            u(c1) = 0.0;
                            for (size_t b = 0; b < nband; b++)
                                u(c1) += M_reg(r, b, c1) * map_inv_var(b, ipix) * dM_reg(r, b, c);
                        }
                        cholesky_solve(A, u, y);
                        grad -= y(c);
                    }
                    gr[r * ncomp + c] += grad;
                }
            }
        skipped[ithr] = nskip;
    });
    auto ll_ = make_Pyarr<double>({nreg});
    auto grad_ = make_Pyarr<double>({nreg, ncomp});
    auto ll_out = to_vmav<double,1>(ll_);
    auto grad_out = to_vmav<double,2>(grad_);
    size_t nskip = 0;
    for (size_t r = 0; r < nreg; r++){
        ll_out(r) = 0.0;
        for (size_t c = 0; c < ncomp; c++) grad_out(r, c) = 0.0;
    }
    for (size_t t = 0; t < nthr; t++){
        nskip += skipped[t];
        for (size_t r = 0; r < nreg; r++){
            ll_out(r) += ll_thr[t * nreg + r];
            for (size_t c = 0; c < ncomp; c++)
                grad_out(r, c) += grad_thr[(t * nreg + r) * ncomp + c];
        }
    }
    if (want_extras)
        return py::make_tuple(ll_, grad_, nskip, x_out_, diag_out_, used_out_);
    return py::make_tuple(ll_, grad_, nskip);
}


void add_compsep(py::module_ &msup)
  {
  using namespace py::literals;
  auto m = msup.def_submodule("compsep");
  m.doc() = "Compiled component-separation kernels.";

  m.def("marginal_loglike_regions", Py_marginal_loglike_regions,
        "Per-region marginal (or ridge) log likelihood and its gradient for a mixing matrix that "
        "is constant within regions: M_reg and dM_reg are (nreg, nband, ncomp), region is the "
        "(npix,) int64 region index (-1 excludes a pixel). Returns (loglike (nreg,), "
        "grad (nreg, ncomp), nskipped); with want_extras also the per-pixel conditional mean x "
        "(ncomp, npix), diag(A^-1) (ncomp, npix) and the used-pixel mask (npix,), zeros where a "
        "pixel did not enter the sum.",
        "map_sky"_a, "map_inv_var"_a, "region"_a, "M_reg"_a, "dM_reg"_a, "nthreads"_a=1,
        "use_logdet"_a=true, "want_grad"_a=true, "want_extras"_a=false);

  m.def("solve_perpix_regions", Py_solve_perpix_regions,
        "Sample component amplitudes (ncomp, npix) per pixel as solve_perpix, with the mixing "
        "matrix M_reg[region[p]] of the pixel's region (M_reg (nreg, nband, ncomp), region (npix,) "
        "int64, -1 gives zero amplitudes). Prior pseudo-bands appended to the stream make the "
        "draw N(A'^-1 b', A'^-1).",
        "map_sky"_a, "map_inv_var"_a, "region"_a, "M_reg"_a, "rand"_a, "nthreads"_a=1);

  m.def("marginal_loglike_perpix", Py_marginal_loglike_perpix,
        "Log marginal likelihood sum_pix [b^T A^-1 b / 2 - ln det A / 2] of the mixing matrix M "
        "(nband, ncomp) for the per-pixel model, with A = M^T N^-1 M and b = M^T N^-1 d built from "
        "map_sky and map_inv_var (nband, npix). use_logdet=False drops the ln det term (the ridge "
        "/ profile likelihood). Returns (loglike, nused, nskipped).",
        "map_sky"_a, "map_inv_var"_a, "M"_a, "nthreads"_a=1, "use_logdet"_a=true);

  m.def("solve_perpix", Py_solve_perpix,
        "Sample component amplitudes (ncomp, npix) by an independent least-squares solve per "
        "pixel, from map_sky and inverse noise variance map_inv_var (nband, npix), mixing matrix M "
        "(nband, ncomp) and standard normal random numbers rand (npix, nband). Unsolvable pixels "
        "get zero amplitudes.",
        "map_sky"_a, "map_inv_var"_a, "M"_a, "rand"_a, "nthreads"_a=1);
  }


}  // ends `namespace detail_pymodule_compsep`

using detail_pymodule_compsep::add_compsep;

}  // ends `namespace cmdr4`
