"""1233: actual regional energy, boundary flux, and autonomous controller source.

Finite CAR matrices diagnose algebra only. They are NOT a cutoff simulation of
the infinite Dirac vacuum or of the original Gaussian clocks. The latter domain
and full-history assertions are proved in proof.md. Default and --verify read
only; --write-results creates a new receipt without overwriting existing files.
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
TARGET = HERE / "results.json"
SOURCES = [
    "archive_1120_1158/1128/proof.md",
    "archive_1120_1158/1138/proof.md",
    "archive_1120_1158/1150/energy_balance.md",
    "archive_1219_/1228/proof.md",
    "archive_1219_/1229/proof.md",
    "archive_1219_/1230/proof.md",
    "archive_1219_/1230/check.py",
    "archive_1219_/1230/results.json",
    "archive_1219_/_admission/stage_scope_audit_after1230.md",
    "archive_1219_/1233/proof.md",
    "archive_1219_/1233/check.py",
]


def comm(a, b):
    return a @ b - b @ a


def maxabs(a):
    return float(np.max(np.abs(a)))


def tensor(*args):
    out = np.array([[1.0 + 0j]])
    for a in args:
        out = np.kron(out, a)
    return out


def annihilators(n):
    out = []
    for j in range(n):
        a = np.zeros((2**n, 2**n), complex)
        for b in range(2**n):
            if b & (1 << j):
                a[b ^ (1 << j), b] = (-1)**((b & ((1 << j)-1)).bit_count())
        out.append(a)
    return out


def car_and_one_particle():
    n = 4
    cs = annihilators(n)
    ident = np.eye(2**n)
    car_error = max(maxabs(cs[i] @ cs[j].conj().T + cs[j].conj().T @ cs[i]
                           - (ident if i == j else 0))
                    for i in range(n) for j in range(n))
    assert car_error < 1e-13
    # A finite matrix diagnostic, not the spatial Dirac discretization.
    h = np.diag([-2., -1., 1., 2.])
    for i in range(3):
        h[i, i+1] = h[i+1, i] = .7
    ev, vec = np.linalg.eigh(h)
    S = vec[:, ev < 0] @ vec[:, ev < 0].conj().T

    def raw(a):
        return sum((a[i, j]*cs[i].conj().T @ cs[j]
                    for i in range(n) for j in range(n)), np.zeros_like(ident, complex))

    def normal(a):
        return raw(a) - np.trace(S @ a)*ident

    def create(u):
        return sum((u[i]*cs[i].conj().T for i in range(n)), np.zeros_like(ident, complex))

    vacuum = np.zeros(2**n, complex)
    vacuum[0] = 1
    for u in vec[:, ev < 0].T:
        vacuum = create(u) @ vacuum
    assert abs(np.linalg.norm(vacuum)-1) < 1e-13
    hf = normal(h)
    assert np.linalg.eigvalsh(hf)[0] > -1e-13
    assert np.linalg.norm(hf @ vacuum) < 1e-13
    chi = np.diag([1., 1., 0., 0.])
    a = (chi @ h + h @ chi)/2
    j = 1j*comm(h, a)
    source_residual = 0.
    source_bound_ratios = []
    for u in (np.eye(n)[:, 0], np.eye(n)[:, 3]):
        nu = create(u) @ create(u).conj().T
        rhs = create(a @ u) @ create(u).conj().T - create(u) @ create(a @ u).conj().T
        source_residual = max(source_residual, maxabs(comm(normal(a), nu)-rhs))
        denom = 2*np.linalg.norm(a @ u)
        if denom:
            source_bound_ratios.append(float(np.linalg.norm(rhs, 2)/denom))
        else:
            assert np.linalg.norm(rhs) < 1e-13
    j_residual = maxabs(j-.5j*comm(h @ h, chi))
    normal_free_residual = maxabs(1j*comm(hf, normal(a))-normal(j))
    assert max(source_residual, j_residual, normal_free_residual) < 1e-12
    assert max(source_bound_ratios) <= 1+1e-13
    data = dict(status="PASS", diagnostic_modes=n, fock_dimension=2**n,
                filled_negative_modes=int(np.sum(ev < 0)), car_residual=car_error,
                CAR_source_residual=source_residual,
                one_particle_boundary_identity_residual=j_residual,
                normal_ordered_free_identity_residual=normal_free_residual,
                bounded_source_norm_to_bound_ratios=source_bound_ratios,
                HS_kernel_not_assumed_bounded_Fock_generator=True)
    return data, dict(h=h, chi=chi, A=a, J=j, HF=hf, Q=normal(a), FJ=normal(j),
                      n1=cs[0].conj().T @ cs[0], n2=cs[3].conj().T @ cs[3], vacuum=vacuum)


def autonomous_history(d):
    # Fock(16) x controller1(2) x record1(2) x controller2(2) x record2(2).
    # These finite controls ONLY test exact Ward algebra; not the real clocks.
    i2, i16 = np.eye(2), np.eye(16)
    hc1 = np.array([[1., .25], [.25, 1.6]])
    hc2 = np.array([[1.2, -.2], [-.2, 1.8]])
    gc1, gc2 = np.diag([.2, .9]), np.diag([.3, .8])
    k = np.pi/2*np.array([[1., -1.], [-1., 1.]])
    hf = tensor(d["HF"], i2, i2, i2, i2)
    c1 = tensor(i16, hc1, i2, i2, i2)
    c2 = tensor(i16, i2, i2, hc2, i2)
    v1 = tensor(d["n1"], gc1, k, i2, i2)
    v2 = tensor(d["n2"], i2, i2, gc2, k)
    q = tensor(d["Q"], i2, i2, i2, i2)
    j = tensor(d["FJ"], i2, i2, i2, i2)
    h = hf+c1+c2+v1+v2+4*np.eye(256)
    source = 1j*comm(v1+v2, q)
    regional = q+c1+v1
    total_field_source = 1j*comm(v1+v2, hf)
    matrix_errors = dict(
        regional_field=maxabs(1j*comm(h, q)-j-source),
        regional_with_controller=maxabs(1j*comm(h, regional)-j),
        sources_commute=maxabs(comm(v1, v2)),
        global_field=maxabs(1j*comm(h, hf)-total_field_source),
        global_controller=maxabs(1j*comm(h, c1+c2+v1+v2)+total_field_source),
        complementary_boundary=maxabs(1j*comm(h, hf-q+c2+v2)+j),
    )
    assert max(matrix_errors.values()) < 1e-11
    plus = np.array([1., 1.])/np.sqrt(2)
    phase = np.array([1., 1j])/np.sqrt(2)
    zero = np.array([1., 0.])
    initial = tensor(d["vacuum"], plus, zero, phase, zero).ravel()
    assert abs(np.linalg.norm(initial)-1) < 1e-13
    eig, unitary = np.linalg.eigh(h)
    coeff = unitary.conj().T @ initial
    observables = dict(q=q, j=j, source=source, regional=regional,
                       hf=hf, c=c1+c2, v=v1+v2, h=h)
    transformed = {key:unitary.conj().T @ value @ unitary for key,value in observables.items()}

    def expectation(key, time):
        state = np.exp(-1j*time*eig)*coeff
        return float(np.vdot(state, transformed[key] @ state).real)

    t = .7
    nodes, weights = np.polynomial.legendre.leggauss(24)
    integrals = {key:sum(float(w)*t/2*expectation(key, t*(float(x)+1)/2)
                         for x,w in zip(nodes,weights)) for key in ("j", "source")}
    delta = {key:expectation(key,t)-expectation(key,0) for key in observables}
    balance_errors = dict(field_integrated=abs(delta["q"]-integrals["j"]-integrals["source"]),
                          controller_region_integrated=abs(delta["regional"]-integrals["j"]),
                          complete_energy=abs(delta["hf"]+delta["c"]+delta["v"]))
    assert max(balance_errors.values()) < 1e-10
    ts = .37
    at = {key:expectation(key,ts) for key in ("j", "source")}
    deletion = dict(status="PASS", observation_time=ts,
                    true_field_rate=at["j"]+at["source"],
                    incorrectly_free_field_rate=at["j"],
                    deleted_actual_source_residual=abs(at["source"]),
                    time_integrated_deleted_source=abs(integrals["source"]),
                    interaction_energy_change=delta["v"],
                    omitted_interaction_energy_residual=abs(delta["hf"]+delta["c"]),
                    finite_example_only=True)
    assert deletion["deleted_actual_source_residual"] > 1e-4
    assert deletion["time_integrated_deleted_source"] > 1e-4
    assert deletion["omitted_interaction_energy_residual"] > 1e-4
    data = dict(status="PASS", full_diagnostic_dimension=256,
                all_finite_controllers_and_records_retained=True,
                original_Gaussian_clock_simulation=False,
                original_infinite_Dirac_simulation=False,
                operator_identity_residuals=matrix_errors,
                actual_history_duration=t, quadrature_points=24,
                integrated_balance_residuals=balance_errors,
                actual_energy_changes={key:delta[key] for key in ("q","regional","hf","c","v")},
                actual_boundary_inflow_integral=integrals["j"],
                actual_control_source_integral=integrals["source"])
    return data, deletion


def continuum_inputs_and_bounds():
    gamma, delta, T = F(1,10**8), F(1,10**24), F(2)
    pi_upper = F(22,7)
    # Fixed, independently chosen spatial plateaus for both old geometry menus.
    assert F(1,5)+gamma < F(3,10)
    assert F(1,100)+2*gamma < F(1,10)
    assert F(1)-gamma-F(1,100)-gamma > F(1,5)
    # Old mode support propagated by |s-t| <= 2 remains in the larger plateau.
    assert F(1,5)+gamma+T < 3
    assert F(1,2)+F(1,100)+2*gamma+T < 3
    c0, c1, c2 = F(1),F(24),F(900)
    theta_prime = F(15,8)
    assert theta_prime**2*(25+25+100) < c1**2
    assert 6*(25+25+100) == c2
    bA, bJ = c0*400+c1/2, c1*400+c2/2
    q_bound = 4*pi_upper*T/delta*bA
    flux_bound = 4*pi_upper*T/delta*bJ
    source_bound = 4*pi_upper/delta*bA
    old = json.loads((ROOT/"archive_1219_/1230/results.json").read_text(encoding="utf8"))
    assert old["status"] == "PASS" and old["round"] == 1230
    old_information = old["groups"]["actual_finite_information_response"]
    assert old_information["actual_large_minus_small_MI_lower"] > 4e-22
    return dict(status="PASS", original_mode_H_regularities=8,
                original_mode_H_norm_upper=400, original_geometry_tolerance=str(gamma),
                boundary_defined_before_energy_and_entropy=True,
                fixed_small_plateau_centers_z=[2,3],
                small_plateau_halfwidths=[.3,.3,.1],
                small_transition_widths=[.2,.2,.1],
                larger_plateau_center=[0,0,2.5], larger_plateau_halfwidths=[3,3,3],
                all_actual_geometry_supports_and_time_halo_contained=True,
                chi_norm_bounds=[float(c0),float(c1),float(c2)],
                A_on_original_transported_mode_bound=str(bA),
                J_on_original_transported_mode_bound=str(bJ),
                uniform_actual_region_energy_absolute_bound=float(q_bound),
                uniform_actual_boundary_flux_absolute_bound=float(flux_bound),
                uniform_actual_control_source_absolute_bound=float(source_bound),
                bounds_are_finite_sufficient_not_optimized=True,
                regulator_only_used_in_analytic_definition=True,
                finite_R_CAR_quadratic_trace_class=True,
                Q_chi_on_vacuum_not_assumed=True,
                no_C2_quantitative_cutoff_rate_claim=True,
                existing_independent_geometry_actual_MI_gap=old_information["actual_large_minus_small_MI_lower"],
                reused_geometry_record_certificate_not_recomputed=True)


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert isinstance(b,dict) and a.keys()==b.keys(), path
        for key in a:
            compare(a[key],b[key],path+"."+key)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b), path
        for i,(aa,bb) in enumerate(zip(a,b)):
            compare(aa,bb,path+f"[{i}]")
    elif isinstance(a,float):
        assert isinstance(b,(float,int)) and abs(a-b) <= 1e-11+1e-10*abs(a), path
    else:
        assert a==b, path


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-results",action="store_true")
    mode.add_argument("--verify",action="store_true")
    args=parser.parse_args()
    car, data = car_and_one_particle()
    history, deletion = autonomous_history(data)
    continuum = continuum_inputs_and_bounds()
    result=dict(round=1233,status="PASS",check_groups=4,
                groups=dict(CAR_source_boundary_identity=car,
                            actual_autonomous_history_and_regional_balance=history,
                            omitted_source_and_interaction_energy_counterexample=deletion,
                            same_continuum_material_geometry_and_finite_bounds=continuum),
                scope=dict(actual_full_driven_history=True,
                    old_1230_Dirac_material_reference_controls_and_records_reused=True,
                    continuum_domain_proved_analytically_not_numerically=True,
                    local_energy_expectation_difference_not_new_sharp_operator=True,
                    effective_regional_gate_permission_input=True,
                    spatial_apparatus_worldtubes=False,
                    all_apparatus_four_momentum=False,
                    independent_geometry_information_response_reused=True,
                    area_response_dynamics=False,Einstein_dynamics=False,
                    new_energy_measurement_device=False,Unruh_device=False,
                    empirical_experiment=False,stage_complete=False),
                source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES})
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(result,json.loads(TARGET.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1233,status="PASS",check_groups=4,
        maximal_operator_residual=max(history["operator_identity_residuals"].values()),
        maximal_integrated_residual=max(history["integrated_balance_residuals"].values()),
        deleted_source_residual=deletion["deleted_actual_source_residual"],
        omitted_interaction_energy_residual=deletion["omitted_interaction_energy_residual"],
        read_only=not args.write_results),ensure_ascii=False))


if __name__=="__main__":
    main()
