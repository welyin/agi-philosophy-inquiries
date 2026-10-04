"""765: the coupled covariant gauge complex, principal and adjoint checks.

No evolution, Hadamard-state construction or loop calculation is simulated.
Analytic arguments in note765 keep the full background-dependent Hessian.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as bg

HERE = Path(__file__).resolve().parent
TARGET = HERE / "joint_covariant_gauge_complex_results.json"
ETA = np.diag([-1., 1., 1., 1.])
NV, NW = 63, 16


def maximum(x):
    return float(np.max(np.abs(x)))


def basis4():
    tensors = []
    for i in range(4):
        t = np.zeros((4, 4)); t[i, i] = 1.
        tensors.append(t)
    for i in range(4):
        for j in range(i+1, 4):
            t = np.zeros((4, 4)); t[i, j] = t[j, i] = 1/np.sqrt(2)
            tensors.append(t)
    return np.asarray(tensors)


T = basis4()


def pack(h):
    return np.einsum("aij,ij->a", T, h)


def unpack(x):
    return np.einsum("a,aij->ij", x, T)


def trace_reverse(h):
    return h-.5*ETA*np.trace(ETA@h)


def scalar_metric(phi):
    F = 2.-phi@phi/6.
    return np.eye(5)/F+np.outer(phi, phi)/(6*F*F)


def pairings(phi, weights, kappa=1.):
    H = np.zeros((NV, NV)); B = np.zeros((NW, NW))
    H[:10, :10] = np.array([
        [np.sum(a*(ETA@trace_reverse(b)@ETA))/(4*kappa) for b in T]
        for a in T])
    for a, weight in enumerate(weights):
        H[10+4*a:14+4*a, 10+4*a:14+4*a] = weight*ETA
    H[58:, 58:] = scalar_metric(phi)
    B[:4, :4] = ETA/(2*kappa)
    B[4:, 4:] = np.diag(weights)
    return H, B


def principal(k, phi, weights, kappa=1.):
    """Real derivative symbol: replace partial_mu by k_mu, not i k_mu."""
    k = np.asarray(k, float); ku = ETA@k; k2 = float(k@ku)
    H, B = pairings(phi, weights, kappa)
    K = np.zeros((NV, NW))
    for a in range(4):
        xi_lower = ETA[:, a]
        K[:10, a] = pack(np.outer(k, xi_lower)+np.outer(xi_lower, k))
    for a in range(12):
        K[10+4*a:14+4*a, 4+a] = k
    # Minus sign is essential: formal adjoint of a first derivative.
    Kadj = -np.linalg.solve(B, K.T@H)
    P = np.zeros((NV, NV))
    for a, h in enumerate(T):
        F = trace_reverse(h)@ku
        P[:10, a] = pack(k2*h-np.outer(k, F)-np.outer(F, k))
    for a in range(12):
        P[10+4*a:14+4*a, 10+4*a:14+4*a] = k2*np.eye(4)-np.outer(k, ku)
    P[58:, 58:] = k2*np.eye(5)
    return H, B, K, Kadj, P, k2


def principal_check():
    phi = np.array([.13, .64, -.11, .08, .37])
    weights = np.repeat([.73, 1.17, .41], [8, 3, 1])
    rows = []
    for k in ((2., .1, .3, -.2), (1., 1., 0., 0.), (.2, .3, -.7, 1.)):
        H, B, K, Ka, P, k2 = principal(k, phi, weights, kappa=.81)
        D = P-K@Ka; D0 = -Ka@K
        errors = dict(self_adjoint=maximum(H@P-P.T@H),
                      noether=maximum(P@K),
                      field_wave_symbol=maximum(D-k2*np.eye(NV)),
                      gauge_wave_symbol=maximum(D0-k2*np.eye(NW)),
                      intertwining=maximum(D@K-K@D0))
        assert max(errors.values()) < 8e-15
        rows.append(dict(covector=list(k), lorentz_square=k2, errors=errors))
    eigen = np.linalg.eigvalsh(H)
    signature = dict(positive=int(sum(eigen > 1e-10)),
                     negative=int(sum(eigen < -1e-10)))
    assert signature == dict(positive=47, negative=16)
    return dict(rows=rows, field_pairing_signature=signature,
                numeric_coefficients="Declared nondegenerate calibration, not a fit to physical constants.",
                scope="Full second-order spacetime symbols in Einstein-frame variables; no lower-order Hessian simulation.")


def scalar_generators():
    mats = np.zeros((12, 5, 5))
    complex_generators = [1j*t for t in bg.geo.old.T] + [3j*np.eye(2)]
    for a, z in enumerate(complex_generators, start=8):
        mats[a, :4, :4] = np.block([[z.real, -z.imag], [z.imag, z.real]])
    return mats


def coupled_adjoint_check():
    q, psi, _, info, c = bg.completed(12, 1.)
    index = (0, 3, 1); phi = q["phi"][index]; conformal = float(psi[index])
    Ks = scalar_metric(phi)
    dphi = np.zeros((4, 5))
    dphi[0] = conformal**-6 * np.linalg.solve(Ks, q["p"][index])
    dphi[1:] = q["Dphi"][index]/conformal**2
    gen = scalar_generators()
    w = np.repeat(bg.geo.old.PAR["K"], [8, 3, 1])
    # Only the actual color magnetic subblock is needed for this witness.
    # Electric and electroweak curvature blocks are NOT claimed to be simulated.
    curvature = np.zeros((12, 4, 4))
    for a in range(8):
        for i in range(3):
            for j in range(3):
                curvature[a, 1+i, 1+j] = bg.inner(bg.T[a], c["F"][i, j])/conformal**4
    Kzero = np.zeros((NV, NW))
    for a in range(12):
        Kzero[10+4*a:14+4*a, :4] = curvature[a].T
    Kzero[58:, :4] = dphi.T
    Kzero[58:, 4:] = -np.einsum("aij,j->ia", gen, phi)
    H, B = pairings(phi, w)
    adj = np.linalg.solve(B, Kzero.T@H)
    rng = np.random.default_rng(765)
    u = rng.normal(size=NV); a = u[10:58].reshape(12, 4); chi = u[58:]
    direct = np.zeros(NW)
    source = np.einsum("a,arm,am->r", w, curvature, a@ETA)+dphi@Ks@chi
    direct[:4] = 2*ETA@source
    direct[4:] = -np.einsum("i,ij,ajk,k->a", chi, Ks, gen, phi)/w
    error = maximum(adj@u-direct)
    assert error < 1e-14
    v = np.zeros(NV); v[58:] = dphi[0]
    matter_to_diffeo = adj@v
    assert maximum(matter_to_diffeo[:4]) > .001
    assert maximum(adj[:4, 10:58]) > .01
    assert maximum(adj[4:, 58:]) > 1e-10
    # Check directly the derivative of K(phi) along t phi plus tensor terms.
    for t in gen:
        d = t@phi
        dmetric = (np.outer(d, phi)+np.outer(phi, d))/(6*q["F"][index]**2)
        assert maximum(dmetric+t.T@Ks+Ks@t) < 1e-14
    return dict(background_index=list(index), background_F=float(q["F"][index]),
                psi=conformal, scalar_metric_min_eigenvalue=float(np.linalg.eigvalsh(Ks)[0]),
                adjoint_error=error,
                color_magnetic_to_diffeomorphism_max=maximum(adj[:4, 10:58]),
                scalar_to_internal_gauge_max=maximum(adj[4:, 58:]),
                scalar_only_diffeomorphism_gauge_fix= matter_to_diffeo[:4].tolist(),
                blocks_simulated="All scalar gauge/translation zero-order blocks, and actual color magnetic curvature only.",
                all_curvature_blocks_or_PDE_simulated=False)


def gauge_cancellation_check():
    # A separate polynomial-symbol check detects wrong field/gauge adjoint weights.
    phi = np.array([.1, .7, .05, -.1, .4]); w = np.linspace(.6, 1.4, 12)
    H, B, K, Ka, P, k2 = principal([1.3, -.2, .6, .4], phi, w)
    wrong_adjoint = -K.T
    wrong = P-K@wrong_adjoint-k2*np.eye(NV)
    assert maximum(wrong) > .1
    # Parametric gauge coefficient: only the declared coefficient cancels nonminimal terms.
    right = P-K@Ka
    half = P-.5*K@Ka
    assert maximum(right-k2*np.eye(NV)) < 1e-14
    assert maximum(half-k2*np.eye(NV)) > .1
    return dict(wrong_Euclidean_adjoint_defect=maximum(wrong),
                half_gauge_coefficient_defect=maximum(half-k2*np.eye(NV)),
                normalized_gauge_coefficient_error=maximum(right-k2*np.eye(NV)),
                purpose="Independent failure witnesses for adjoint and gauge-normalization conventions; no Hadamard inference.")


def run():
    result = dict(round=765, tests_run=3, failures=0, errors=0,
                  covariant_principal_symbols=principal_check(),
                  original_background_gauge_mixing=coupled_adjoint_check(),
                  convention_failure_witnesses=gauge_cancellation_check())
    dependencies = ("research_note_573.md", "research_note_731.md",
                    "research_note_735.md", "research_note_752.md",
                    "research_note_753.md", "research_note_754.md",
                    "research_note_764.md", "joint_reference_constraint_strata.py",
                    "joint_gauss_einstein_initial_data.py")
    result["dependency_hashes"] = {
        n: hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies}
    result["scope"] = ("Coupled linear covariant gauge complex, causal physical algebra and "
                       "canonical matching proved analytically on the inherited smooth on-shell "
                       "Einstein-frame model. Numerics check second-order symbols and specified "
                       "zero-order gauge blocks. No physical Hadamard-state existence, quantum "
                       "Ward renormalization, all loops, original instrument matching or Q-to-E "
                       "continuum equivalence is claimed.")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--write-results", action="store_true")
    args = p.parse_args(); result = run()
    if args.write_results:
        with TARGET.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text("utf8"))
    print(json.dumps({k: result[k] for k in ("round", "tests_run", "failures", "errors")}))
