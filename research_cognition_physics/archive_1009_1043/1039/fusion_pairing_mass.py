"""1039: full fusion data, neutral-pair flip, and quadratic Weyl mass.

All finite group entries are Gaussian integers. Exact-equality checks below
involve only these integers and small dyadic averages; generic-basis and
perturbation checks use NumPy tolerances and do not prove universal claims.
"""
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "fusion_pairing_mass_results.json"
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
EPS = np.array([[0, 1], [-1, 0]], complex)
FLIP = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                 [0, 1, 0, 0], [0, 0, 0, 1]], complex)


def groups():
    ans = {}
    for name, r, s, form, sign in [
        ("D8", -EPS, Z, I2, 1),
        ("Q8", 1j * Z, 1j * X, EPS, -1),
    ]:
        labels = [(k, ell) for k in range(4) for ell in range(2)]
        mats = [np.linalg.matrix_power(r, k) @ np.linalg.matrix_power(s, ell)
                for k, ell in labels]
        chars = np.array([[(-1)**(a*k+b*ell) for k, ell in labels]
                          for a in range(2) for b in range(2)] +
                         [[np.trace(u) for u in mats]], complex)
        ans[name] = dict(mats=mats, generators=[r, s], labels=labels,
                         chars=chars, form=form, sign=sign)
    return ans


def serialize_matrix(a):
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def symmetric_basis(size):
    out = []
    for i in range(size):
        for j in range(i, size):
            e = np.zeros((size, size), complex)
            e[i, j] = e[j, i] = 1
            out.append(e)
    return out


def invariant_mass_nullity(gens, n):
    size = 2*n
    reps = [np.kron(g, np.eye(n)) for g in gens]
    basis = symmetric_basis(size)
    a = np.column_stack([np.concatenate([(u.T @ e @ u-e).ravel()
                                         for u in reps]) for e in basis])
    singular = np.linalg.svd(a, compute_uv=False)
    return len(basis)-int(np.sum(singular > 1e-9))


def mass_example(n, sign):
    if sign == 1:
        return np.diag(np.arange(1, n+1)).astype(complex)
    a = np.zeros((n, n), complex)
    for k in range(n//2):
        a[2*k, 2*k+1] = k+1
        a[2*k+1, 2*k] = -k-1
    return a


def run():
    gs = groups()
    exact = 0
    def eq(a, b):
        nonlocal exact
        assert np.array_equal(a, b), (a, b)
        exact += 1

    result = {}
    fusions = {}
    pplus, pminus = (np.eye(4)+FLIP)/2, (np.eye(4)-FLIP)/2
    pointer_u = np.kron(pplus, I2)+np.kron(pminus, X)
    pointer_h = np.pi/2*np.kron(pminus, I2-X)
    ev, vectors = np.linalg.eigh(pointer_h)
    exp_h = (vectors*np.exp(-1j*ev)) @ vectors.conj().T
    assert opnorm(exp_h-pointer_u) < 1e-14
    eq(pointer_u.conj().T @ pointer_u, np.eye(8))
    for name, data in gs.items():
        mats, chars, form, sign = [data[k] for k in
                                  ("mats", "chars", "form", "sign")]
        assert len({tuple(u.ravel()) for u in mats}) == 8
        eq(chars @ chars.conj().T/8, np.eye(5))
        eq(sum(int(chars[i, 0].real)**2 for i in range(5)), 8)
        table = []
        for u in mats:
            eq(u.conj().T @ u, I2)
            eq(u.T @ form @ u, form)
            table.append([])
            for v in mats:
                matches = [i for i, w in enumerate(mats)
                           if np.array_equal(u @ v, w)]
                assert len(matches) == 1
                table[-1].append(matches[0])
                exact += 1
        tensor_reps = [np.kron(u, u) for u in mats]
        p0 = sum(tensor_reps)/8
        vector = form.conj().ravel()/np.sqrt(2)
        # Avoid irrational normalization in exact projector checks.
        eq(p0, np.outer(form.conj().ravel(), form.ravel())/2)
        eq(p0 @ p0, p0)
        eq(np.trace(p0), 1)
        eq(FLIP @ p0, sign*p0)
        indicator = sum(np.trace(u @ u) for u in mats)/8
        eq(indicator, sign)
        fusion = np.empty((5, 5, 5), int)
        for i in range(5):
            for j in range(5):
                coeff = (chars[i]*chars[j]) @ chars.conj().T/8
                eq(coeff, np.round(coeff.real))
                fusion[i, j] = np.round(coeff.real).astype(int)
        fusions[name] = fusion
        for u in tensor_reps:
            eq(u @ p0, p0)
            eq(FLIP @ u, u @ FLIP)
            whole = np.kron(u, I2)
            eq(pointer_u @ whole, whole @ pointer_u)
        initial = np.kron(vector, np.array([1, 0], complex))
        final = pointer_u @ initial
        pointer_z = float(np.vdot(final, np.kron(np.eye(4), Z) @ final).real)
        assert abs(pointer_z-sign) < 1e-14
        masses = []
        for n in range(1, 9):
            a = mass_example(n, sign)
            m = np.kron(form, a)
            eq(m.T, m)
            for u in mats:
                rep = np.kron(u, np.eye(n))
                eq(rep.T @ m @ rep, m)
            nullity = invariant_mass_nullity(data["generators"], n)
            expected_dimension = n*(n+sign)//2
            assert nullity == expected_dimension
            svals = np.linalg.svd(m, compute_uv=False)
            null_mass_modes = int(np.sum(svals < 1e-10))
            assert null_mass_modes == (0 if sign == 1 else 2*(n % 2))
            masses.append(dict(n=n, mass_space_complex_dimension=nullity,
                               singular_values=svals.tolist(),
                               zero_quadratic_mass_components=null_mass_modes))
        result[name] = {
            "group_order": 8, "nonidentity_involutions":
                sum(np.array_equal(u @ u, I2) for u in mats)-1,
            "character_table_on_elements": chars.real.astype(int).tolist(),
            "multiplication_table": table, "fusion_coefficients": fusion.tolist(),
            "frobenius_schur_indicator": int(indicator.real),
            "trivial_pair_projector": serialize_matrix(p0),
            "pointer_contrast": pointer_z, "mass_cases": masses,
        }
    eq(fusions["D8"], fusions["Q8"])
    eq(gs["D8"]["chars"], gs["Q8"]["chars"])
    # Faithful V generates every simple: V^2 is the sum of all four characters.
    eq(fusions["D8"][4, 4], np.array([1, 1, 1, 1, 0]))
    wrong_p = np.outer(EPS.ravel(), EPS.ravel())/2
    d_p0 = sum(np.kron(u, u) for u in gs["D8"]["mats"])/8
    for u in gs["D8"]["mats"]:
        t = np.kron(u, u)
        eq(t @ wrong_p @ t.conj().T, wrong_p)
    eq(np.trace(d_p0 @ wrong_p), 0)
    eq(np.trace(FLIP @ wrong_p), -1)

    # Adverse leakage occupies the opposite flip eigenspace.
    errors = []
    for eta, epsilon in [(0, 0), (.1, .05), (.2, .1), (.49, .01), (.5, 0)]:
        margin = 1-2*eta-epsilon
        for name, data in gs.items():
            p = sum(np.kron(u, u) for u in data["mats"])/8
            q = wrong_p if name == "D8" else d_p0
            rho = (1-eta)*p+eta*q
            contrast = float(np.trace(FLIP @ rho).real)
            worst = contrast-data["sign"]*epsilon
            assert abs(data["sign"]*worst-margin) < 1e-14
        errors.append(dict(eta=eta, contrast_error=epsilon,
                           signed_lower_bound=margin,
                           strict_classification=bool(margin > 0)))

    # Generic complex basis and approximate Ward checks, fixed seed.
    rng = np.random.default_rng(1039)
    random_checks = []
    for n in [1, 3, 5, 7]:
        for trial in range(5):
            a = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
            a = (a-a.T)/2
            m0 = np.kron(EPS, a)
            v, _ = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))
            w = np.kron(v, np.eye(n))
            reps = [w.conj().T @ np.kron(g, np.eye(n)) @ w
                    for g in gs["Q8"]["mats"]]
            m0 = w.T @ m0 @ w
            noise = rng.normal(size=m0.shape)+1j*rng.normal(size=m0.shape)
            noise = (noise+noise.T)/2
            noise *= .03/opnorm(noise)
            m = m0+noise
            images = [u.T @ m @ u for u in reps]
            bar = sum(images)/8
            delta = max(opnorm(x-m) for x in images)
            distance = opnorm(m-bar)
            small = np.sort(np.linalg.svd(m, compute_uv=False))[:2]
            assert distance <= delta+1e-12
            assert max(small) <= delta+1e-12
            ward_bar = max(opnorm(u.T @ bar @ u-bar) for u in reps)
            assert ward_bar < 1e-12
            assert np.sort(np.linalg.svd(bar, compute_uv=False))[1] < 1e-12
            random_checks.append(dict(n=n, trial=trial, delta=delta,
                                      distance_to_invariant=distance,
                                      two_smallest_singular_values=small.tolist(),
                                      averaged_ward_residual=ward_bar))
    return {
        "round": 1039, "exact_equality_checks": exact, "groups": result,
        "same_full_fusion_ring": True,
        "density_invariance_negative_control": {
            "D8_antisymmetric_density_invariant": True,
            "trivial_sector_weight": 0, "flip_contrast": -1},
        "pointer_hamiltonian_min_eigenvalue": float(min(ev)),
        "pointer_hamiltonian_norm": opnorm(pointer_h),
        "pointer_unitary_exponential_residual": opnorm(exp_h-pointer_u),
        "preparation_and_readout_cases": errors,
        "approximate_mass_cases": random_checks,
        "universal_claim_basis": "analytic proof; numerical cases are calibration",
        "full_anomaly_or_QFT_completion": False,
        "independent_agent_review": "pending",
    }


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert abs(a-b) <= 1e-11*(1+abs(a))
    else:
        assert a == b, (a, b)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        compare(result, json.loads(OUT.read_text("utf8")))
    print(json.dumps({"passed": True, "round": 1039,
                      "exact_equality_checks": result["exact_equality_checks"],
                      "mass_cases": 16, "approximate_mass_cases": 20,
                      "fusion_ring_equal": True}, ensure_ascii=False))
