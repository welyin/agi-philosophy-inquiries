"""1133 autonomous regional receiver: faithful full-CAR diagnostics and budgets.

Default is read-only. --write-results writes only this round's results.json.
No finite matrix below simulates the continuous source or replaces its proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = {
    "archive_1120_/1129/proof.md": "c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa",
    "archive_1120_/1129/results.json": "69ed793bd5847d31003eb05b18205ecf72af34c26f9ad2919ebe548e7a25a02c",
    "archive_1120_/1132/proof.md": "f05daf49972a0b7fe13f5d593a032998bf1c9008182c051f2e12bb6ea958df63",
    "archive_1120_/1132/results.json": "95dcbd80a2a9e5d3d7ba03e52f1d1aa7e931606c9735f0359623f7b4f969430c",
    "archive_1120_/1122/proof.md": "52cd6812b52dbbdfe16d4557d88918db239e80900e3ae70f8f6b05a35b52d8ce",
    "archive_1120_/1123/proof.md": "3d1786db4efc63e807fc8375d0c658b5b1e56f31cb301fb8ffe0ab52b2248401",
    "archive_935_955/research_note_939.md": "d7a6ce1aaf9231622999c2301732cd60d78a85ad33ea4e477ca64f8acf2d7eb3",
    "archive_935_955/research_note_943.md": "684bba213dc9b8521a5e91ec42ce2296a58880cc895868f949cb8dd9d41b8653",
}


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T


def sources():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes == FROZEN
    old = json.loads((ROOT/"archive_1120_/1132/results.json").read_text(encoding="utf-8"))
    assert old["source_identity"]["complete_old_clock_dimension"] == 288
    assert old["exact_budgets"]["fractions"]["mass"] == str(10**80)
    return dict(source_sha256=hashes, complete_old_clock_dimension=288,
                new_selector_dimension=2, new_pointer_dimension=2,
                common_positive_receiver_internal_mass=1,
                same_normal_scalar_source=True, same_original_normal_CAR_code=True,
                calibration_W_applied=False, previous_author_assets_modified=False,
                new_permission="M1133: explicit flavor-sensitive regional capture/record/continuation material, independent selector and pointer preparation, terminal physical pointer readout; no point-local claim.")


def budgets():
    eps, A, lam, kap, T = F(1, 10**5), F(10**34), F(11, 7), F(11, 14), F(4)
    sigma, M, dc, tc = F(1, 10**30), F(10**80), F(1, 10**6), F(5, 4)
    q = 8/eps**3
    C = 8*q*524
    Ceta = F(24780800000000000000000, 49)
    Crec = 68*42+8884+lam*(68*16+16*40)+kap*(68*16+16*80)
    assert Crec == 16316
    CY = Ceta+Crec
    source_derivative = lam*C*512+42*8/dc+(lam+kap)*512
    ds = 3*T/(eps*A)*source_derivative
    dy = CY*(2*T*sigma+T*T/(2*M*sigma))
    db = 28*T/(2*M*M*sigma*sigma)
    cap, dm = 42*dc, 80*lam*eps
    recv = cap+ds+dy+db
    total = dm+recv
    vrec = 42+16*lam+16*kap
    vall = 16*lam*C+vrec
    ypos = 2*sigma+T/(M*sigma)+CY*T*T/(2*M)
    rstar = 5+eps
    assert 1+T+eps < 7 and 7-rstar > F(3, 4)
    assert 6-eps > T and tc+dc < F(9, 4)
    assert vrec == F(558, 7) and total < F(1303, 10**6)
    assert 2*(dy+db) < F(81, 10**10)
    assert db == F(56, 10**100)
    vals = dict(epsilon=eps, amplitude=A, time=T, capture_time=tc,
                capture_width=dc, sigma=sigma, mass=M, new_receiver_internal_mass=F(1),
                old_plus_receiver_internal_mass_norm_upper=F(28),
                C_upper=C, C_eta_upper=Ceta, C_receiver_upper=Crec, C_Y_upper=CY,
                Q0_norm_upper=F(42), Q0_spatial_derivative_upper=F(8884),
                chi_gradient_L1_upper=F(68), grad_hg_upper=F(1801),
                source_derivative_budget=source_derivative,
                source_joint_error=ds, support_joint_error=dy,
                internal_mass_joint_error=db, capture_joint_error=cap,
                receiver_joint_error=recv, material_handoff_error=dm,
                complete_task_error=total, external_causal_error=2*(dy+db),
                V_receiver_norm_upper=vrec, V_all_norm_upper=vall,
                late_source_mean_lower=7-rstar,
                support_recoil_RMS_increment=CY*T,
                support_output_position_RMS_upper=ypos,
                support_outside_unit_ball_probability_upper=ypos*ypos,
                initial_support_kinetic_mean_upper=3/(8*M*sigma*sigma),
                initial_source_energy_upper=1073*A*A,
                output_free_energy_mean_upper=268+1073*A*A+28+3/(8*M*sigma*sigma)+2*vall,
                physical_energy_lower_weak=M-vall)
    return dict(scope="Exact rational upper bounds; dy/db quantify fixed independent Gaussian support, ds preserves the whole coherent source output.",
                fractions={k: str(v) for k, v in vals.items()},
                decimal_views={k: float(v) for k, v in vals.items()},
                record_window=[2.25, 2.75], continuation_window=[3.25, 3.75],
                source_error_needed_for_external_causal_certificate=False,
                Gaussian_tail_truncated=False, endpoint_small_interaction_energy_inherited=False,
                terminal_pointer_readout_permission_still_required=True,
                interacting_boost_or_complete_CO_SR_claimed=False)


def quad(fn, a, b, n=100):
    x, w = np.polynomial.legendre.leggauss(n)
    return float((b-a)/2*np.sum(w*fn((a+b)/2+(b-a)*x/2)))


def bump01(x):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    mask = (x > 0)&(x < 1)
    out[mask] = np.exp(-1/(x[mask]*(1-x[mask])))
    return out


IB = quad(bump01, 0, 1, 160)


def ramp(u):
    u = float(u)
    if u <= 0:
        return 0.
    if u >= 1:
        return 1.
    return quad(bump01, 0, u, 90)/IB


def gate(x):
    return 2*bump01(2*np.asarray(x)-.5)/IB


def geometry_and_kernels():
    q = lambda x: np.exp(-1/(1-x*x))
    qp = lambda x: -2*x/(1-x*x)**2*q(x)
    qpp = lambda x: (4*x*x/(1-x*x)**4-2/(1-x*x)**2-8*x*x/(1-x*x)**3)*q(x)
    n0, n1, n2 = quad(lambda x:q(x)**2, -1, 1, 180), quad(lambda x:qp(x)**2, -1, 1, 180), quad(lambda x:qpp(x)**2, -1, 1, 180)
    grad_f = np.sqrt(3*n1/n0)
    laplacian_sq = 3*n2/n0+6*(n1/n0)**2  # integration by parts in cross terms
    grad_hf = np.sqrt(laplacian_sq+grad_f**2)
    assert grad_f < 20 and np.sqrt(laplacian_sq) < 1800 and grad_hf < 1801
    rad = lambda r: np.exp(-1/(1-r*r))
    z = quad(lambda r:4*np.pi*r*r*rad(r), 0, 1, 180)
    gradchi = quad(lambda r:4*np.pi*r*r*(2*r/(1-r*r)**2)*rad(r), 0, 1, 180)/z
    assert z > .125 and gradchi < 68
    assert abs(quad(gate, .25, .75, 160)-1) < 1e-12
    assert abs(quad(lambda t:gate(t-2), 2.25, 2.75, 160)-1) < 1e-12
    assert abs(quad(lambda t:gate(t-3), 3.25, 3.75, 160)-1) < 1e-12
    ramp_area = quad(lambda u:np.array([1-ramp(z) for z in u]), 0, 1, 70)
    assert abs(ramp_area-.5) < 2e-12  # symmetry: sharper than the conservative proof
    # Independent smooth beta0 spherical averages, not just the cap formula.
    # beta0(r)=1-R(r-7); at radius r the massless Cauchy solution is t times
    # its spherical average. Integrate the exact smooth plateau numerically.
    rows = []
    cos_nodes, cos_weights = np.polynomial.legendre.leggauss(100)
    for radius in [0., 3., 5., 5.00001]:
        for time in [1.25, 1.99999, 2.5, 3., 4.]:
            radii = np.sqrt(radius**2+time**2+2*radius*time*cos_nodes)
            beta = np.array([1-ramp(rr-7) for rr in radii])
            mean = time*np.dot(cos_weights, beta)/2
            assert mean > .75
            if time <= 1.99999:
                assert abs(mean-time) < 1e-12
            else:
                rs = 5.00001
                lower = time/2+(49-rs*rs-time*time)/(4*rs)
                assert mean >= lower-1e-10
            rows.append(dict(radius=radius, time=time, smooth_source_mean=float(mean)))
    return dict(scope="Independent bump derivative integrals, actual normalized pulse areas, and smooth spherical source means; analytic inequalities remain the certificate.",
                f_gradient_norm=float(grad_f), f_laplacian_norm=float(np.sqrt(laplacian_sq)),
                gradient_hf_norm=float(grad_hf), rho_normalization=z,
                chi_gradient_L1=gradchi, normalized_gate_area=quad(gate, .25, .75, 160),
                capture_ramp_missing_area=ramp_area,
                old_source_same_normal_preparation=True, smooth_source_rows=rows)


def annihilator(mode, modes=4):
    out = np.zeros((2**modes, 2**modes), complex)
    for state in range(2**modes):
        if (state >> mode)&1:
            out[state^(1 << mode), state] = (-1)**((state&((1 << mode)-1)).bit_count())
    return out


def car_instruments():
    # Faithful full algebra of two captured and two leakage modes. The theorem
    # is all-input, so this check never replaces the actual sea by an empty mode.
    c = [annihilator(i) for i in range(4)]
    ident = np.eye(16)
    lift = lambda a: sum((a[i,j]*c[i].conj().T@c[j] for i in range(4) for j in range(4)), np.zeros((16,16), complex))
    a = .8
    h1 = np.array([[0,0,a,0],[0,0,0,a],[a,0,.7,0],[0,a,0,1.2]], complex)
    coupling = h1.copy()
    coupling[2,2] = coupling[3,3] = 0
    offset = -np.sum(np.minimum(np.linalg.eigvalsh(h1), 0))
    hf, q0 = lift(h1)+offset*ident, lift(coupling)
    hold = hf-q0
    assert np.linalg.eigvalsh(hf).min() > -1e-13
    holding_residual = max(norm(hold@ci-ci@hold) for ci in c[:2])
    assert holding_residual < 1e-13 and abs(norm(q0)-2*a) < 1e-13
    nz = c[0].conj().T@c[0]
    cx = (c[0]+c[1])/np.sqrt(2)
    nx = cx.conj().T@cx
    yg = 1j*(c[1].conj().T@c[0]-c[0].conj().T@c[1])
    number = c[0].conj().T@c[0]+c[1].conj().T@c[1]
    assert abs(norm(nz@nx-nx@nz)-.5) < 1e-13
    uc = exp_h(yg, np.pi/4)
    uh = exp_h(hold, 4-1.25)
    xp = np.array([[0,1],[1,0]], complex)
    ready = np.array([[1],[0]], complex)
    inject = np.kron(ident, ready)
    rows = []
    for label, n in [("Z",nz),("X",nx)]:
        record_u = exp_h(np.kron(n,xp), np.pi/2)
        actual_iso = np.kron(uh@uc, np.eye(2))@record_u@inject
        exact_iso = np.kron(uh@uc@(ident-n), ready)-1j*np.kron(uh@uc@n, np.array([[0],[1]]))
        err = norm(actual_iso-exact_iso)
        assert err < 1e-13
        kraus = [actual_iso.reshape(16,2,16)[:,r,:] for r in range(2)]
        assert norm(sum(k.conj().T@k for k in kraus)-ident) < 1e-13
        assert norm(kraus[1].conj().T@kraus[1]-n) < 1e-13
        assert max(norm(k@number-number@k) for k in kraus) < 1e-13
        # Every full-CAR matrix unit and its passive-reference counterpart.
        matrix_unit_residual = 0.
        for i in range(16):
            for j in range(16):
                for r in range(2):
                    expected_k = uh@uc@([ident-n,-1j*n][r])
                    got = np.outer(kraus[r][:,i],kraus[r][:,j].conj())
                    expected = np.outer(expected_k[:,i],expected_k[:,j].conj())
                    matrix_unit_residual = max(matrix_unit_residual,norm(got-expected))
        assert matrix_unit_residual < 1e-13
        # Maximally entangled full 16D input/reference tests complete output,
        # retaining coherent pointer phase; dropping -i is observably different.
        bell_out = (actual_iso/4).reshape(-1)
        bell_target = (exact_iso/4).reshape(-1)
        bell_error = float(np.linalg.norm(bell_out-bell_target))
        wrong_phase = np.kron(uh@uc@(ident-n), ready)+np.kron(uh@uc@n,np.array([[0],[1]]))
        assert norm(actual_iso-wrong_phase) > 1
        if label == "Z":
            input_state = c[0].conj().T[:,0]
        else:
            input_state = cx.conj().T[:,0]
        conditional = kraus[1]@input_state
        same_effect_after = float(np.real(conditional.conj()@n@conditional))
        assert abs(same_effect_after-.5) < 1e-13
        rows.append(dict(setting=label, all_input_isometry_residual=err,
            full_matrix_unit_instrument_residual=matrix_unit_residual,
            entangled_reference_dimension=16, coherent_reference_residual=bell_error,
            eigenbranch_same_effect_after_continuation=same_effect_after,
            wrong_coherent_phase_operator_error=norm(actual_iso-wrong_phase)))
    # Finite-duration capture diagnostic using time ordering, independent of
    # simply exponentiating the claimed target. A moderate width resolves error.
    width = .03
    ramp_propagators = []
    for steps in [80,160]:
        out = ident.astype(complex)
        dt = width/steps
        for j in range(steps):
            out = exp_h(hf-ramp((j+.5)/steps)*q0, dt)@out
        ramp_propagators.append(out)
    hard = exp_h(hold,width)
    capture_err = norm(ramp_propagators[-1]-hard)
    assert capture_err < norm(q0)*width
    convergence = norm(ramp_propagators[1]-ramp_propagators[0])
    assert convergence < 2e-7
    return dict(scope="Faithful four-mode/full-16D CAR algebra, all occupancy sectors and all matrix units; actual continuum vacuum/source proof remains analytic.",
                CAR_modes=4, full_CAR_dimension=16, logical_only_test=False,
                holding_commutator_residual=holding_residual,
                two_effect_commutator_norm=norm(nz@nx-nx@nz), rows=rows,
                positive_width_capture=dict(diagnostic_width=width, steps=[80,160],
                    full_propagator_error=capture_err, Duhamel_upper=norm(q0)*width,
                    step_doubling_residual=convergence),
                records_and_all_field_poststates_retained=True,
                continuum_source_or_finite_CCR_simulation_claimed=False)


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i,(a,b) in enumerate(zip(actual,expected)):
            compare(a,b,path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual,expected,atol=2e-10,rtol=2e-10), path
    else:
        assert actual == expected, path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write-results",action="store_true")
    args = p.parse_args()
    result = dict(round=1133,schema=1,status="PASS_AUTONOMOUS_REGIONAL_RECEIVER_FULL_CAR_INSTRUMENTS",
                  new_adopted_axioms=0,scientific_calibration=0,stage="active",
                  numerical_check_groups=4,source_identity=sources(),exact_budgets=budgets(),
                  kernels_and_long_source_window=geometry_and_kernels(),full_CAR_instruments=car_instruments())
    target = HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("PASS: four groups; only 1133/results.json registered explicitly.")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
        print("PASS: four groups; eight frozen source hashes and results match; no files written.")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
