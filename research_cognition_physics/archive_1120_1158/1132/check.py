"""1132 normal moving support and the complete old clock as internal mass.

Default is read-only. --write-results writes only this round's results.json.
Gaussian quadrature and finite momentum-transfer matrices are independent
mechanism diagnostics, not simulations or proofs of the continuous QFT.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = {
    "archive_1120_/1129/proof.md": "c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa",
    "archive_1120_/1129/results.json": "69ed793bd5847d31003eb05b18205ecf72af34c26f9ad2919ebe548e7a25a02c",
    "archive_1106_1119/1106/proof.md": "3cf66b629a49e00a6503b68e516f67bdd5b9fecd258447e90156b15cf87fb8cb",
    "archive_1106_1119/1106/results.json": "1fa96fcf58ce22ccf529aebf0fd8d0eb8dd4edb1bdb449998e0ccaec2bbea7e2",
    "archive_1120_/1130/proof.md": "9e4c2c64f062db6ef61777f532753a2a9cf3b0ab13df4280dcc3c6c35519a6ee",
    "archive_1120_/1130/results.json": "e8cc423302b5fda96d509b5a3f04cbe572874cb22955877d68be5602e211ba62",
    "archive_935_955/research_note_942.md": "0f51b118cb170d547f9ec5edbdef1c78bca5c35ba8379b56e60904b98356374e",
    "archive_935_955/research_note_944.md": "8a74182ea5b0aaccc68bd850e06dd39069a310ad26d414ccf5d081b8ecf8e057",
    "archive_1120_/1125/proof.md": "434513dc989e62fe77bf6f4af3ea8a779c381bfea25a047273c96dcf5c02addf",
}


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e)) @ v.conj().T


def sources():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes == FROZEN
    old = json.loads((ROOT/"archive_1120_/1130/results.json").read_text(encoding="utf-8"))
    assert old["source_identity"]["frozen_old_clock_dimension"] == 288
    assert old["exact_budgets"]["fractions"]["time"] == "5/4"
    return dict(source_sha256=hashes, old_author_assets_modified=False,
                old_large_clock_reexecuted=False, complete_old_clock_dimension=288,
                same_normal_CAR_code=True, same_normal_scalar_source=True,
                calibration_W_applied=False,
                new_material_permission="M1132: normal independent quantum support, complete old D installed as its internal mass; original finite endpoint effects permitted.")


def budgets():
    eps, amp, lam, t = F(1, 10**5), F(10**34), F(11, 7), F(5, 4)
    sigma, mass, Bclock = F(1, 10**30), F(10**80), F(27)
    q, c = 8/eps**3, 8*(8/eps**3)*524
    grad_eta, ceta = F(2200, 7), 16*lam*8*q*F(2200, 7)
    dm, ds = 80*lam*eps, lam*c*512*3*t/(eps*amp)
    dy = ceta*(2*t*sigma+t*t/(2*mass*sigma))
    db = Bclock*t/(2*mass**2*sigma**2)
    dmov = dy+db
    r10 = lam**10/F(math.factorial(10))/(1-lam/11)
    position = 2*sigma+t/(mass*sigma)+ceta*t*t/(2*mass)
    vbound = 16*lam*c
    initial_kinetic = 3/(8*mass*sigma**2)
    free_energy = 268+1073*amp**2+27+initial_kinetic+2*vbound
    assert ceta == F(24780800000000000000000, 49)
    assert dy == F(38720000000000000000121, 30625000000000000000000000000000)
    assert db == F(16875, 10**103)
    assert dm+ds+dmov < F(1258155951, 10**12)
    assert 4*r10+2*ds+2*dmov < F(120109954, 10**12)
    assert position < F(2000001, 10**36)
    assert position**2 < F(4000005, 10**66)
    assert 5+t+eps < 7 and 2+t+eps < 4 and mass-vbound > 0
    vals = dict(epsilon=eps, amplitude=amp, time=t, mass=mass, sigma=sigma,
                q_upper=q, C_upper=c, eta_gradient_L1_upper=grad_eta,
                C_eta_upper=ceta, V_norm_upper=vbound,
                material_error=dm, source_error=ds, support_error=dy,
                clock_mass_correction_error=db, complete_moving_error=dmov,
                joint_handoff_error=dm+ds+dmov, finite_causal_error=4*r10+2*ds+2*dmov,
                old_clock_effect_error=dmov, initial_support_kinetic_mean=initial_kinetic,
                initial_support_kinetic_norm_upper=1/(2*mass*sigma**2),
                support_recoil_RMS_increment=ceta*t,
                support_output_momentum_RMS_upper=1/sigma+ceta*t,
                support_output_position_RMS_upper=position,
                support_outside_unit_ball_probability_upper=position**2,
                initial_source_energy_upper=1073*amp**2,
                free_energy_mean_upper=free_energy,
                physical_energy_mean_upper=mass+free_energy+vbound,
                physical_energy_lower_weak=mass-vbound)
    return dict(scope="Strict analytic rational bounds; pi<22/7, sqrt3<2, sqrt15<4. Tiny B is input-restricted, not an operator norm.",
                fractions={k: str(v) for k, v in vals.items()},
                decimal_views={k: float(v) for k, v in vals.items()},
                natural_old_clock_reference_probability=float(np.sin((np.pi/2+1.25)/2)**6),
                physical_rest_mass_included=True, support_tail_projected=False,
                old_endpoint_small_interaction_energy_inherited=False,
                spatial_translation_strongly_conserved=True,
                interacting_boost_constructed=False, complete_CO_SR_qualification=False)


def gaussian_diagnostics():
    # Independent Cartesian Hermite integration of standard normal momentum.
    # sigma=.5 implies each momentum coordinate has unit standard deviation.
    rows = []
    for n in [28, 48]:
        nodes, weights = np.polynomial.hermite.hermgauss(n)
        p = np.sqrt(2)*nodes
        w = weights/np.sqrt(np.pi)
        p2 = p[:, None, None]**2+p[None, :, None]**2+p[None, None, :]**2
        weights3 = w[:, None, None]*w[None, :, None]*w[None, None, :]
        avg = lambda x: float(np.sum(weights3*x))
        assert abs(avg(np.ones_like(p2))-1) < 1e-13
        assert abs(avg(p2)-3) < 1e-12 and abs(avg(p2*p2)-15) < 1e-11
        for mass in [3., 12., 100.]:
            energy = np.sqrt(mass*mass+p2)
            kinetic = p2/(energy+mass)  # stable, including the large mass limit
            v2 = p2/(mass*mass+p2)
            # Momentum derivative of e^{-it T(p)} psi(p) gives position norm.
            t, sig = 1.25, .5
            direct_position_sq = avg(4*sig**4*p2+t*t*v2)
            independently_heisenberg = 3*sig**2+t*t*avg(v2)
            assert abs(direct_position_sq-independently_heisenberg) < 1e-12
            assert avg(kinetic) <= 3/(2*mass)
            assert np.sqrt(avg(kinetic**2)) <= np.sqrt(15)/(2*mass)
            assert avg(v2) <= 3/mass**2 and np.max(v2) < 1
            # Direct square roots versus independent integral in internal mass.
            d = 2.3
            b_direct = np.sqrt((mass+d)**2+p2)-energy-d
            z, zw = np.polynomial.legendre.leggauss(36)
            b_integral = np.zeros_like(p2)
            for zz, ww in zip((z+1)*d/2, zw*d/2):
                b_integral += ww*((mass+zz)/np.sqrt((mass+zz)**2+p2)-1)
            err = float(np.max(np.abs(b_integral-b_direct)))
            assert err < 8e-14 and np.max(b_direct) <= 1e-13
            assert np.min(b_direct) >= -d-1e-13
            assert avg(b_direct*b_direct)**.5 <= d*np.sqrt(15)/(2*mass**2)
            rows.append(dict(order=n, diagnostic_mass=mass,
                kinetic_mean=avg(kinetic), kinetic_RMS=np.sqrt(avg(kinetic**2)),
                speed_squared_mean=avg(v2), position_squared=direct_position_sq,
                internal_mass_integral_residual=err,
                internal_mass_correction_RMS=np.sqrt(avg(b_direct*b_direct))))
    for a, b in zip(rows[:3], rows[3:]):
        assert abs(a["kinetic_mean"]-b["kinetic_mean"]) < 2e-9
        assert abs(a["speed_squared_mean"]-b["speed_squared_mean"]) < 2e-9
    # Actual enormous M and sigma require no subtraction of nearly equal energies.
    # Show the fourth-moment bound and normality analytically in exact_budgets.
    return dict(scope="Normal Gaussian continuum integrals at moderate diagnostic masses; exact actual-parameter certification is rational and analytic, no finite CCR truncation.",
                diagnostic_sigma=.5, momentum_second_moment=3., momentum_fourth_moment=15.,
                orders=[28, 48], rows=rows,
                free_position_checked_by_two_representations=True,
                all_actual_Gaussian_tails_retained=True)


def momentum_transfer_diagnostic():
    # Three physical factors: finite field momentum, finite support momentum,
    # a source two-level register with zero spatial momentum. This is NOT a
    # representation of the actual CAR or scalar CCR; only a transfer mechanism.
    kval, pval, nr = range(-2, 3), range(-3, 4), 2
    labels = [(k, p, r) for k in kval for p in pval for r in range(nr)]
    index = {label: i for i, label in enumerate(labels)}
    dim = len(labels)
    k = np.array([z[0] for z in labels], float)
    p = np.array([z[1] for z in labels], float)
    r = np.array([z[2] for z in labels], float)
    pf, py, pt = np.diag(k), np.diag(p), np.diag(k+p)
    md, gd, duration = 5., .65, .9
    h0 = np.diag(np.sqrt(1+k*k)+p*p/(np.sqrt(md*md+p*p)+md)+.7*r)
    v = np.zeros((dim, dim), complex)
    edges = []
    for i, (ki, pi, ri) in enumerate(labels):
        target = (ki+1, pi-1, 1-ri)
        if target in index:
            j = index[target]
            v[i, j] = v[j, i] = gd
            edges.append((i, j))
    assert all(labels[i][0]+labels[i][1] == labels[j][0]+labels[j][1] for i, j in edges)
    h = h0+v
    force = 1j*(h@py-py@h)
    field_force = 1j*(h@pf-pf@h)
    assert norm(force+field_force) < 1e-13
    assert norm(h@pt-pt@h) < 1e-13
    residuals = []
    for shift in [.13, .79, 2.21]:
        trans = np.diag(np.exp(-1j*shift*(k+p)))
        residuals.append(norm(trans@h@trans.conj().T-h))
    assert max(residuals) < 1e-13
    u = exp_h(h, duration)
    # Independent Gauss-Legendre time integral of the force vs endpoint change.
    xn, wn = np.polynomial.legendre.leggauss(40)
    integral = np.zeros_like(h)
    e, ev = np.linalg.eigh(h)
    for tt, ww in zip((xn+1)*duration/2, wn*duration/2):
        ut = (ev*np.exp(-1j*tt*e))@ev.conj().T
        integral += ww*(ut.conj().T@force@ut)
    recoil = u.conj().T@py@u-py
    residual = norm(recoil-integral)
    assert residual < 2e-13
    assert norm(recoil+u.conj().T@pf@u-pf) < 2e-13
    assert norm(u.conj().T@h@u-h) < 2e-13
    # Source is entangled with an untouched reference; every output retained.
    support = np.exp(-((np.array(list(pval))-.4)/1.3)**2/4)
    support = support/np.linalg.norm(support)
    psi = np.zeros((dim, 2), complex)
    for ip, pp in enumerate(pval):
        psi[index[(-1, pp, 0)], 0] = .6*support[ip]
        psi[index[(-1, pp, 1)], 1] = .8*support[ip]
    out = u@psi
    mean = lambda mat, state: float(np.real(np.einsum("ir,ij,jr->", state.conj(), mat, state)))
    dpy, dpf = mean(py, out)-mean(py, psi), mean(pf, out)-mean(pf, psi)
    assert abs(dpy+dpf) < 2e-13 and abs(dpy) > .01
    before = psi.reshape(5, 7, 2, 2)
    after = out.reshape(5, 7, 2, 2)
    support_rho = lambda a: np.einsum("kprz,kqrz->pq", a, a.conj())
    source_rho = lambda a: np.einsum("kprz,kpsz->rs", a, a.conj())
    ry0, ry1, ra0, ra1 = support_rho(before), support_rho(after), source_rho(before), source_rho(after)
    trace_distance = lambda a,b: float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)
    dy, da = trace_distance(ry0, ry1), trace_distance(ra0, ra1)
    purity_y = float(np.real(np.trace(ry1@ry1)))
    assert dy > .01 and da > .001 and purity_y < .99
    assert norm(out.conj().T@out-psi.conj().T@psi) < 1e-13
    return dict(scope="Finite open-boundary momentum-transfer diagnostic only. It has exact conserved momentum and a kept source/support/reference, not the continuous QFT or its normal vacuum.",
                dimension=dim, passive_reference_dimension=2, momentum_edges=len(edges),
                translation_conjugation_max_residual=max(residuals),
                recoil_integral_residual=residual,
                support_mean_momentum_change=dpy, field_mean_momentum_change=dpf,
                source_state_trace_distance=da, support_state_trace_distance=dy,
                output_support_purity=purity_y,
                total_energy_change=mean(h, out)-mean(h, psi),
                all_output_factors_retained=True,
                nonzero_actual_continuum_mean_recoil_claimed=False)


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual, expected, atol=2e-11, rtol=2e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = dict(round=1132, schema=1, status="PASS_MOVING_SUPPORT_INTERNAL_CLOCK_AND_SPATIAL_RECOIL",
                  new_adopted_axioms=0, scientific_calibration=0, stage="active",
                  numerical_check_groups=4, source_identity=sources(), exact_budgets=budgets(),
                  normal_Gaussian=gaussian_diagnostics(),
                  finite_momentum_transfer=momentum_transfer_diagnostic())
    target = HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print("PASS: four groups; only 1132/results.json registered explicitly.")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        print("PASS: four groups; nine frozen source hashes and recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
