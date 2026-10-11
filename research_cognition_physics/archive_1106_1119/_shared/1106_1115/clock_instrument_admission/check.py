"""Non-numbered admission: finite relational clock and neutral internal records.
Default recomputes and checks saved results without writing.
This is not an autonomous incoming/outgoing collision realization.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def op(a):
    return float(np.linalg.norm(a, 2))


def unitary(h, t=1):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T


def spectral(h, f):
    e, v = np.linalg.eigh(h)
    return (v*f(e))@v.conj().T


def run():
    i2, i4 = np.eye(2), np.eye(4)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    n = np.diag([0, 1])
    omega, dr, dm = .8, 2., 3.
    hc, ha = omega*np.kron(n, i2), omega*np.kron(i2, n)
    hca = hc+ha
    b = np.zeros((4, 4), complex)
    b[1, 2], b[2, 1] = 1j, -1j
    fplus, fminus = (i4+b)/2, (i4-b)/2
    aplus = spectral(fplus, lambda z: np.sqrt(np.clip(z, 0, None)))
    aminus = spectral(fminus, lambda z: np.sqrt(np.clip(z, 0, None)))
    theta = spectral(b, lambda z: np.arccos(np.clip(z, -1, 1))/2)
    kr = np.kron(theta, y)/dr
    ur = unitary(kr, dr)
    ur_formula = np.kron(aplus, i2)-1j*np.kron(aminus, y)
    km = np.pi/(2*dm)*np.kron(n, x)
    um = unitary(km, dm)
    w = np.kron(i4, um)@np.kron(ur, i2)
    hbdry = np.kron(hca, i4)
    ready = np.array([1, 0, 0, 0], complex)
    src = np.kron(i4, ready[:, None])
    target = np.column_stack([np.kron(aplus[:, j], ready)-1j*np.kron(aminus[:, j], np.array([0, 0, 0, 1])) for j in range(4)])
    g1 = dict(read_generator_norm=op(kr), copy_generator_norm=op(km),
        total_phase_readout_dimension=16, relational_effect_energy_commutator=op(b@hca-hca@b),
        readout_formula_error=op(ur-ur_formula), unitary_error=op(w.conj().T@w-np.eye(16)),
        complete_instrument_identity_error=op(w@src-target),
        complete_boundary_energy_intertwining_error=op(hbdry@w-w@hbdry))
    assert max(g1[k] for k in ("relational_effect_energy_commutator", "readout_formula_error", "unitary_error", "complete_instrument_identity_error", "complete_boundary_energy_intertwining_error")) < 1e-12

    rows = []
    for phi in (-np.pi/3, -np.pi/6, 0., np.pi/6, np.pi/3):
        psi = np.array([0, 1, np.exp(-1j*phi), 0], complex)/np.sqrt(2)
        state = w@np.kron(psi, ready)
        amps = state.reshape(4, 2, 2)
        pplus = float(np.sum(abs(amps[:, 0, 0])**2))
        mismatch = float(np.sum(abs(amps[:, 0, 1])**2)+np.sum(abs(amps[:, 1, 0])**2))
        expected = (1+np.sin(phi))/2
        er_var = float(np.vdot(psi, theta@theta@psi).real/(dr*dr))
        ec_var = (np.pi/(2*dm))**2*(1-expected)
        assert abs(pplus-expected) < 1e-12 and mismatch < 1e-12
        rho_deph = np.diag(abs(psi)**2)
        assert abs(np.trace(fplus@rho_deph).real-.5) < 1e-12
        rows.append(dict(phase=phi, reported_plus=pplus, analytic_plus=float(expected),
                         record_disagreement_probability=mismatch,
                         dephased_reference_plus=.5,
                         first_contact_energy_variance=er_var,
                         copying_contact_energy_variance=float(ec_var)))
    g2 = dict(readout_duration=dr, internal_record_duration=dm,
              phase_branch=[-np.pi/3, np.pi/3], omega=omega,
              minimum_probability_slope=omega/4,
              probability_error_to_time_error_multiplier=4/omega,
              two_candidate_phases=[-np.pi/6, np.pi/6],
              one_record_equal_prior_discrimination_error=.25,
              initial_total_clock_energy=omega, initial_total_clock_energy_variance=0.,
              each_clock_energy_mean=omega/2, each_clock_energy_variance=omega*omega/4,
              clock_reference_energy_covariance=-omega*omega/4, cases=rows)

    # Unknown Q is a genuine arbitrary finite spectator, not a chosen test input.
    hq = np.diag([0., .3, .7])
    full_w = np.kron(w, np.eye(3))
    full_h = np.kron(hbdry, np.eye(3))+np.kron(np.eye(16), hq)
    kraus = (np.kron(aplus, np.eye(3)), -1j*np.kron(aminus, np.eye(3)))
    complete = sum(a.conj().T@a for a in kraus)
    assert op(complete-np.eye(12)) < 1e-12
    assert op(full_h@full_w-full_w@full_h) < 1e-12
    # Fixed boundary rest masses are declared independently; all h, including Q, is retained.
    baseline_mass = 10.+11.+12.+13.+17.
    mass = baseline_mass*np.eye(48)+full_h
    velocities = (np.array([1., 0., 0., 0.]), np.array([1.25, .75, 0., 0.]))
    four_errors = [op((u_component*mass)@full_w-full_w@(u_component*mass))
                   for u in velocities for u_component in u]
    # Two different reference histories; free common waiting at the read vertex does not change B.
    free_wait = unitary(hca, dr+dm)
    free_read_error = op(free_wait.conj().T@b@free_wait-b)
    g3 = dict(payload_dimension=3, all_input_kraus_completeness_error=op(complete-np.eye(12)),
              complete_payload_energy_intertwining_error=op(full_h@full_w-full_w@full_h),
              independently_declared_baseline_mass=baseline_mass,
              minimum_boundary_mass=float(np.min(np.linalg.eigvalsh(mass))),
              common_velocity_four_momentum_intertwining_max_error=max(four_errors),
              common_wait_relational_effect_error=free_read_error,
              arbitrary_passive_reference_follows_from_factor_identity=True,
              actual_fresh_payload_loading_implemented=False,
              energy_of_unimplemented_trigger_or_support_not_signed_off=True)
    assert max(four_errors+[free_read_error]) < 1e-11
    paths = ["archive_935_955/research_note_939.md", "archive_956_989/research_note_976.md",
             "archive_370_428/research_note_420.md", "archive_1106_/1108/proof.md"]
    return dict(candidate="finite_relational_clock_record_instrument", formal_round=None,
                scientific_calibration_increment=0, numerical_groups=3,
                status="PASS_OPERATOR_AND_BOUNDARY_LEDGER_ONLY",
                group1_instrument=g1, group2_clock_readings_and_resources=g2,
                group3_payload_and_mass_boundary=g3,
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                scope=dict(effective_positive_duration_collision_is_an_extra_material_permission=True,
                           autonomous_incoming_exit_trigger_not_constructed=True,
                           complete_local_SR_experiment_not_completed=True,
                           preparation_phase_reference_and_record_blanks_are_resources=True,
                           no_postselection=True, no_battery_optimization=True))


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys()==b.keys(), path
        for k in a: compare(a[k], b[k], path+"."+k)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a)==len(b), path
        for i, (aa, bb) in enumerate(zip(a,b)): compare(aa,bb,path+str(i))
    elif isinstance(a, float):
        assert math.isfinite(a) and math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-11),path
    else: assert a==b,path


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args();data=run();target=HERE/"results.json"
    if args.write_results:
        assert not target.exists(),"Refusing to overwrite saved candidate results"
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    compare(data,json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps({"candidate":data["candidate"],"status":data["status"],
                      "groups":data["numerical_groups"],"saved_results_match":True,
                      "formal_round":None,"scientific_calibration_increment":0},ensure_ascii=False))
