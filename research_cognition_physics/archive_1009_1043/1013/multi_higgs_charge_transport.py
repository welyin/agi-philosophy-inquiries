"""Round 1013: common unbroken generator and multi-doublet charge transport.

One finite calibration group. The general implications are proved in the note;
these finite matrices do not derive the gauge group or a vacuum dynamically.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "multi_higgs_charge_transport_results.json"
I2 = np.eye(2, dtype=complex)
SIGMA = [np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex),
         np.diag([1., -1.]).astype(complex)]
T = [s / 2 for s in SIGMA]
EPS = np.array([[0, 1], [-1, 0]], complex)


def op(a):
    return float(np.linalg.norm(a, ord=2))


def gram(vevs, charges, kinetic=None):
    """Quadratic form x^T G x = ||i(sum x_j T_j + x_3 Y)v||_K^2.

    Gauge couplings are absorbed in x. A positive gauge-field kinetic form can
    change nonzero mass eigenvalues but not this kernel. No 1/2 mass convention
    is claimed for the numerical G eigenvalues.
    """
    v = np.asarray(vevs, complex)
    h = np.asarray(charges, float)
    k = np.eye(len(v)) if kinetic is None else np.asarray(kinetic, complex)
    assert op(k - k.conj().T) < 1e-13
    assert op(k @ np.diag(h) - np.diag(h) @ k) < 1e-13, "kinetic mixes inequivalent charges"
    assert min(np.linalg.eigvalsh(k)) >= -1e-13
    cols = [np.concatenate([1j * t @ z for z in v]) for t in T]
    cols.append(np.concatenate([1j * q * z for q, z in zip(h, v)]))
    a = np.array(cols).T
    w = np.kron(k, I2)
    return np.real(a.conj().T @ w @ a), a, w


def nullity(a):
    return int(np.sum(np.linalg.svd(a, compute_uv=False) < 1e-10))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_old(number, filename):
    p = HERE.parent / str(number) / filename
    spec = importlib.util.spec_from_file_location("old_" + str(number), p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run():
    old1011 = load_old(1011, "flavor_charge_selection.py")
    rows = []
    k = np.array([[2., .2+.3j], [.2-.3j, 1.5]])
    cases = [
        ("aligned_equal_charge_nondiagonal_positive_kinetic", [[0, 1], [0, 1+1j]], [.5, .5], k, 1),
        ("opposite_charge_opposite_eigenlines", [[0, 1], [2, 0]], [.5, -.5], None, 1),
        ("same_charge_orthogonal_vevs_break_mixed_U1", [[0, 1], [1, 0]], [.5, .5], None, 0),
        ("same_charge_nonorthogonal_noncollinear_vevs", [[0, 1], [1, 1]], [.5, .5], None, 0),
        ("aligned_unequal_charge_magnitudes", [[0, 1], [0, 1]], [.5, 1.], None, 0),
        ("zero_vev_spectator_arbitrary_charge", [[0, 1], [0, 2], [0, 0]], [.5, .5, 7/3], None, 1),
        ("only_zero_charge_vevs_leave_pure_Y_exception", [[0, 1]], [0.], None, 1),
    ]
    for name, vevs, charges, kinetic, expected in cases:
        g, a, w = gram(vevs, charges, kinetic)
        assert nullity(g) == expected
        assert nullity(np.vstack([a.real, a.imag])) == expected
        # Independent realification of the complex kinetic quadratic form.
        ar = np.vstack([a.real, a.imag])
        wr = np.block([[w.real, -w.imag], [w.imag, w.real]])
        err = op(g - ar.T @ wr @ ar)
        assert err < 1e-12
        rows.append(dict(case=name,charges=charges,gram_eigenvalues=np.linalg.eigvalsh(g).tolist(),
                         common_kernel_dimension=expected,realification_error=err))

    # An invariant kinetic form cannot mix unequal charges in this complex basis.
    invalid_rejected = False
    try:
        gram([[0, 1], [0, 1]], [.5, 1.], [[1, .2], [.2, 1]])
    except AssertionError as e:
        invalid_rejected = "inequivalent" in str(e)
    assert invalid_rejected

    # Positive semidefinite is insufficient: the ignored doublet can transform.
    gsemi, asemi, _ = gram([[0, 1], [1, 0]], [.5, .5], [[1, 0], [0, 0]])
    q = np.array([0., 0., 1., 1.])
    assert abs(q @ gsemi @ q) < 1e-14
    assert np.linalg.norm(asemi @ q) > .99
    # The first three cases use genuine positive-definite kinetic matrices.
    assert min(np.linalg.eigvalsh(k)) > 1

    # Conjugating the opposite-charge field preserves the real mass quadratic.
    original = [np.array([0, 1], complex), np.array([2, 0], complex)]
    canonical = [original[0], EPS @ original[1].conjugate()]
    g0, _, _ = gram(original, [.5, -.5])
    g1, _, _ = gram(canonical, [.5, .5])
    conjugation_error = op(g0-g1)
    assert conjugation_error < 1e-14

    rotations = []
    for j in range(1, 17):
        axis = np.array([1., j/7, -.4])
        axis /= np.linalg.norm(axis)
        angle = j/9
        u = np.cos(angle/2)*I2 - 1j*np.sin(angle/2)*sum(x*s for x,s in zip(axis,SIGMA))
        nz = u @ T[2] @ u.conj().T
        n = np.array([float(np.real(2*np.trace(nz @ t))) for t in T])
        gr, ar, _ = gram([u@z for z in original], [.5, -.5])
        qr = np.r_[n, 1.]
        residual = float(np.linalg.norm(ar @ qr))
        spectral_error = float(np.max(np.abs(np.linalg.eigvalsh(gr)-np.linalg.eigvalsh(g0))))
        assert residual < 1e-12 and spectral_error < 1e-12
        rotations.append(dict(case=j,annihilation_residual=residual,spectral_error=spectral_error))

    determinants = []
    for rho, h in [([F(1),F(2)],[F(1,2),F(1,2)]),
                   ([F(1),F(2)],[F(1,2),F(-1,2)]),
                   ([F(2),F(3),F(5)],[F(1,2),F(1),F(3,2)]),
                   ([F(1),F(4),F(0)],[F(1,2),F(1,2),F(7,3)])]:
        s0, s1, s2 = sum(rho),sum(r*x for r,x in zip(rho,h)),sum(r*x*x for r,x in zip(rho,h))
        det = (s0*s2-s1*s1)/4
        pair = sum(rho[i]*rho[j]*(h[i]-h[j])**2 for i in range(len(h)) for j in range(i+1,len(h)))/4
        assert det == pair
        determinants.append(dict(weights=list(map(str,rho)),charges=list(map(str,h)),determinant=str(det)))

    # Three individually rank-one Yukawas generate full-rank physical masses.
    mu = np.diag([1.,2.,4.])
    md = np.array([[2.,1.,0.],[0.,2.,1.],[1.,0.,3.]])
    me = np.diag([.5,1.,3.])
    maps = []
    largest = 0.
    ql = np.eye(3)/6
    pl = -np.eye(3)/2
    for label,m,left,right,shift in [
        ("up",mu,ql,2*np.eye(3)/3,-.5),
        ("down",md,ql,-np.eye(3)/3,.5),
        ("charged_lepton",me,pl,-np.eye(3),.5)]:
        pieces = []
        for a in range(3):
            p = np.zeros((3,3)); p[a,a]=1
            y = p @ m
            assert np.linalg.matrix_rank(y) == 1
            err = op(left@y-y@(right+shift*np.eye(3)))
            largest = max(largest,err)
            pieces.append(y)
        assert np.linalg.matrix_rank(sum(pieces)) == 3
        assert op(sum(pieces)-m) == 0
        maps.append(dict(sector=label,individual_yukawa_ranks=[1,1,1],total_mass_rank=3))
    assert largest < 1e-12
    dim = old1011.commutant_dimension(mu@mu.T, md@md.T)
    assert dim == 1

    # Individual multi-Higgs kappa^{ab} need not be flavor symmetric.
    # Exchange of (flavor,Higgs) pairs and VEV contraction give symmetric Mnu.
    mnu = np.diag([0.,2.,5.]).astype(complex)
    a = np.array([[1j,2+1j,0],[.4,1,3j],[2,0,-1]], complex)
    kappas = [[mnu-a-a.T, a], [a.T, np.zeros((3,3),complex)]]
    effective = sum(kappas[i][j] for i in range(2) for j in range(2))
    assert op(a-a.T) > 1
    assert op(effective-effective.T) < 1e-13 and op(effective-mnu) < 1e-13
    assert np.max(np.abs(np.linalg.eigvalsh(effective.conj().T@effective)-[0,4,25])) < 1e-12

    # Nonstandard-charge branch transports too, with the old paired-mass obstruction.
    x = np.diag([1.,-1.,0.])
    b = np.array([[0,1+.2j,0],[.4j,0,0],[0,0,.7]],complex)
    assert op(x.T@b+b@x) == 0
    pair = b+b.T
    assert op(x.T@pair+pair@x) == 0
    vals = np.linalg.eigvalsh(pair.conj().T@pair)
    assert abs(vals[0]-vals[1]) < 1e-12

    sources = [
        "archive_935_955/research_note_949.md",
        "archive_1009_/research_note_1010.md",
        "archive_1009_/1010/hypercharge_neutrino_selection.py",
        "archive_1009_/research_note_1011.md",
        "archive_1009_/1011/flavor_charge_selection.py",
        "archive_956_989/981/drafts/common_parent_contract_v1.md",
    ]
    return dict(round=1013,date="2026-10-08",new_calibration_groups=1,
        cumulative_research_groups=3791,all_scientific_calibrations_passed=True,
        gram_cases=rows,kinetic_inequivalent_charge_mixing_rejected=invalid_rejected,
        semidefinite_counterexample=dict(gram_quadratic=0.,nonzero_vev_variation_norm=float(np.linalg.norm(asemi@q))),
        conjugate_doublet_gram_error=conjugation_error,rotated_vacua=rotations,
        exact_aligned_determinants=determinants,yukawa_transport=maps,
        largest_yukawa_intertwiner_residual=largest,total_quark_commutant_dimension=dim,
        multi_higgs_weinberg=dict(individual_kappa_nonsymmetry=op(a-a.T),
            effective_symmetry_error=op(effective-effective.T),
            effective_mass_squared=np.linalg.eigvalsh(effective.conj().T@effective).tolist(),
            nonstandard_charge_mass_squared=vals.tolist()),
        claim_boundaries=dict(SM_group_generated=False,generation_count_generated=False,
            Higgs_count_selected=False,vacuum_stability_proved=False,mass_mechanism_generated=False,
            new_cognitive_axiom=False,all_physics_generated=False,goal_completed=False,
            global_gauge_structure_certified=False),
        scope="Positive kinetic form and a stipulated common mixed unbroken generator constrain VEV directions/charges; physical total masses and a stipulated multi-Higgs Majorana operator transport 1010-1011.",
        historical_source_sha256={name:sha(BASE/name) for name in sources},code_sha256=sha(Path(__file__)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if k in ["round","new_calibration_groups",
         "cumulative_research_groups","all_scientific_calibrations_passed","total_quark_commutant_dimension"]},ensure_ascii=False))
