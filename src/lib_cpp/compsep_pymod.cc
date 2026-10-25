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
    execStatic(npix, nthreads, 0, [&](Scheduler &sched){
        // Per-thread workspace, reused for every pixel.
        vmav<double,2> A({ncomp, ncomp});
        vmav<double,1> rhs({ncomp}), x({ncomp});
        while (auto rng = sched.getNext())
            for (size_t ipix = rng.lo; ipix < rng.hi; ipix++){
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
                        rhs(c1) += M(b, c1) * weighted_d;
                        for (size_t c2 = 0; c2 < ncomp; c2++)
                            A(c1, c2) += M(b, c1) * inv_var * M(b, c2);
                    }
                }
                const bool success = solve_cholesky(A, rhs, x);
                for (size_t c = 0; c < ncomp; c++)
                    comp_maps(c, ipix) = success ? x(c) : 0.0;
            }
    });
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


void add_compsep(py::module_ &msup)
  {
  using namespace py::literals;
  auto m = msup.def_submodule("compsep");
  m.doc() = "Compiled component-separation kernels.";

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
