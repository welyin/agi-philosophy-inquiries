"""953: common neutral-body worldline sources and exact internal propagation.

Smooth prescribed fields / leading monopole source identities only.
Does not certify point-source quantum field dynamics or inherit 947 errors.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / "worldline_material_bridge_results.json"

def read(p): return json.loads(p.read_text("utf-8-sig"))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a): return float(np.linalg.norm(a, 2))
def unitary(h):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j*w)) @ v.conj().T

def run():
    old = read(STAGE / "947/protocol_field_transport_results.json")
    portal = read(STAGE / "946/portal_common_process_results.json")
    body = read(STAGE / "951/body_sector_bridge_results.json")
    assert all(r["all_scientific_checks_passed"] for r in (old, portal, body))
    spec = importlib.util.spec_from_file_location("protocol929", STAGE / "929/joint_protocol_selection.py")
    protocol = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(protocol)
    core = np.eye(protocol.DIM, 4, dtype=complex)
    for gate in protocol.circuit(np.pi/2, (2, 3))[:2]:
        core = protocol.act(core, gate)
    assert np.linalg.norm(core[4:]) == 0
    gcore = core[:4]
    C = gcore.conj().T @ np.diag([1., 1., -1., -1.]) @ gcore
    scale = old["same_physical_protocol"]["internal_clock_scale"]
    off = .5*np.sqrt(np.arange(1,31)*np.arange(30,0,-1))
    hc = scale*(15*np.eye(31)+np.diag(off,1)+np.diag(off,-1))
    h = np.kron(hc, np.eye(4))
    B = np.kron(np.eye(31), C)
    identity = np.eye(124)
    checks = dict(actual_B_square_error=norm(B@B-identity),
                  actual_h_B_commutator_error=norm(h@B-B@h))
    assert max(checks.values()) < 1e-12
    full_dim = 31*2**17
    assert full_dim == old["same_physical_protocol"]["original_work_and_fault_and_clock_dimension"]
    T, charge, diagnostic_mass = 30., .5, 20.
    # Nonuniform parameterization of one path. Fields below are diagnostic
    # prescribed histories; this is not a self-consistent Einstein solution.
    nodes, weights = np.polynomial.legendre.leggauss(120)
    lam, weights = .5*T*(nodes+1), .5*T*weights
    tau = lam + .12*np.sin(2*np.pi*lam/T)
    ein = 1 + .12*2*np.pi/T*np.cos(2*np.pi*lam/T)
    phi = .3 + .2*np.sin(2*np.pi*tau/T)
    lapse_variation = .15 + .1*np.sin(2*np.pi*lam/T)
    field_variation = .4 + .1*np.cos(2*np.pi*tau/T)
    proper_duration = float(weights @ ein)
    integral_phi = float(weights @ (ein*phi))
    action = (diagnostic_mass*identity+h)*proper_duration + charge*B*integral_phi
    U = unitary(action)
    expected = np.exp(-1j*diagnostic_mass*T)*np.kron(unitary(T*hc), unitary(charge*.3*T*C))
    checks["reparameterized_proper_duration_error"] = abs(proper_duration-T)
    checks["reparameterized_field_integral_error"] = abs(integral_phi-.3*T)
    checks["full_internal_unitary_error"] = norm(U-expected)
    checks["unitarity_error"] = norm(U.conj().T@U-identity)
    start = np.zeros((124,4),complex)
    start[:4] = np.eye(4)
    finish = np.zeros((124,4),complex)
    finish[-4:] = np.exp(-1j*diagnostic_mass*T)*unitary(charge*.3*T*C)
    checks["all_input_clock_endpoint_isometry_error"] = norm(U@start-finish)
    # Vary geometric proper-time density and scalar source together.
    def action_at(s):
        e = ein + s*lapse_variation
        p = phi + s*field_variation
        return (diagnostic_mass*identity+h)*float(weights@e)+charge*B*float(weights@(e*p))
    analytic = ((diagnostic_mass*identity+h)*float(weights@lapse_variation)
        +charge*B*float(weights@(lapse_variation*phi+ein*field_variation)))
    step = 1e-4
    fd = (action_at(-2*step)-8*action_at(-step)+8*action_at(step)-action_at(2*step))/(12*step)
    checks["independent_action_source_variation_error"] = norm(fd-analytic)
    dud = (unitary(action_at(-2*step))-8*unitary(action_at(-step))
           +8*unitary(action_at(step))-unitary(action_at(2*step)))/(12*step)
    checks["unitary_source_insertion_relative_error"] = norm(1j*U.conj().T@dud-analytic)/norm(analytic)
    omitted = charge*B*float(weights@(lapse_variation*phi))
    wrong_source_gap = norm(omitted)
    # Compare covariant four-force to independent coordinate Hamilton equations
    # for every energy/charge branch of the original internal operators.
    eta = np.diag([-1.,1.,1.,1.])
    v = np.array([.25,-.1,.05])
    gamma = 1/math.sqrt(1-v@v)
    u = gamma*np.r_[1.,v]
    grad = np.array([.17,.11,-.2,.07])  # covector (dt p, dx p, ...)
    phi_here = .35
    force_error = 0.
    transverse_error = 0.
    bad_residual = 0.
    for energy in scale*np.arange(31):
        for sign in (-1.,1.):
            q = charge*sign
            mass = diagnostic_mass+energy+q*phi_here
            dmass = q*float(u@grad)
            acceleration = -q/mass*(eta+np.outer(u,u))@grad
            covariant_dp = mass*acceleration+u*dmass
            momentum = mass*u[1:]
            E = math.sqrt(float(momentum@momentum)+mass**2)
            coordinate_dp = gamma*np.r_[q*mass/E*grad[0], -q*mass/E*grad[1:]]
            force_error = max(force_error, float(np.linalg.norm(covariant_dp-coordinate_dp)))
            force_error = max(force_error, float(np.linalg.norm(covariant_dp+q*eta@grad)))
            transverse_error = max(transverse_error, abs(float(u@eta@acceleration)))
            # Using only m+h in T while keeping the actual trajectory/source
            # leaves this local coefficient of the distributional Ward defect.
            missing = q*phi_here*acceleration+u*dmass
            bad_residual = max(bad_residual, float(np.linalg.norm(missing)))
    checks["covariant_vs_Hamilton_four_force_error"] = force_error
    checks["acceleration_orthogonality_error"] = transverse_error
    # The original two scalar mass matrix belongs to the same invariant portal.
    pars = portal["parameters"]
    l, vH, a, mu = (pars[k] for k in ("lambda_H","v","portal_alpha","mu"))
    K = np.array([[mu**2+2*l*a*a,2*l*a*vH],[2*l*a*vH,2*l*vH*vH]])
    checks["same_gauge_invariant_portal_Hessian_error"] = norm(K-np.array(pars["mass_matrix"]))
    rng = np.random.default_rng(953)
    doublet = rng.normal(size=2)+1j*rng.normal(size=2)
    mat = rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    Q, _ = np.linalg.qr(mat)
    def potential(z,p):
        return float(l*(np.vdot(z,z).real-vH*vH/2+a*p)**2+.5*mu*mu*p*p)
    checks["portal_U2_rotation_error"] = abs(potential(Q@doublet,.23)-potential(doublet,.23))
    # Finite-size matching: a point vertex is not the original Gaussian vertex.
    shape_rows = []
    for k in (.1,1.,4.):
        x = .5*k*k # sigma=1, inherited from 944
        single = -math.expm1(-x)
        double = -math.expm1(-2*x)
        assert single <= x and double <= 2*x
        shape_rows.append(dict(k_sigma=k,single_vertex_relative_difference=single,
            two_vertex_kernel_relative_difference=double,
            single_vertex_upper=x,two_vertex_upper=2*x))
    assert checks["full_internal_unitary_error"] < 2e-11
    assert checks["all_input_clock_endpoint_isometry_error"] < 2e-11
    assert checks["independent_action_source_variation_error"] < 1e-7
    assert checks["unitary_source_insertion_relative_error"] < 1e-7
    assert force_error < 1e-14 and transverse_error < 1e-14
    assert checks["same_gauge_invariant_portal_Hessian_error"] < 1e-14
    assert checks["portal_U2_rotation_error"] < 1e-12
    assert wrong_source_gap > .5 and bad_residual > .05
    files = [Path(__file__),STAGE/"929/joint_protocol_selection.py",
        STAGE/"947/protocol_field_transport_results.json",STAGE/"946/portal_common_process_results.json",
        STAGE/"951/body_sector_bridge_results.json",STAGE/"952/portal_lapse_source_results.json"]
    return dict(round=953,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(actual_internal_dimension=full_dim,computed_core_dimension=124,
            spectator_multiplicity=2**15,diagnostic_base_mass=diagnostic_mass,
            inherited_physical_base_mass=1e14,duration=T,g_source=charge,
            proper_field_integral=integral_phi,same_portal_K=K.tolist()),
        checks=checks,negative_controls=dict(omitted_interaction_geometric_source_norm=wrong_source_gap,
            omitted_interaction_mass_Ward_coefficient_max=bad_residual),
        finite_size_comparison=shape_rows,
        scope=dict(same_actual_h_B_worldline_internal_unitary=True,
            full_unknown_internal_state_and_passive_reference_transport=True,
            scalar_and_metric_source_exchange_identity_at_monopole_order=True,
            neutral_body_preserves_bulk_gauge_invariance=True,
            no_fundamental_Dirac_species_per_internal_state_required=True,
            point_body_not_identified_with_original_smeared_body=True,
            original_947_probability_bound_inherited=False,
            coupled_quantum_body_field_gravity_process_certified=False,
            self_force_or_point_quantum_UV_completed=False,
            microscopic_SM_binding_constructed=False,
            dimensions_or_Einstein_or_SM_derived_from_cognition=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args=parser.parse_args()
    if args.write:
        assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2)
            f.write("\n")
    else:
        old=read(TARGET)
        assert old["source_hashes"]==out["source_hashes"]
        assert old["scope"]==out["scope"]
        for k,val in out["checks"].items():
            assert abs(val-old["checks"][k])<1e-8,k
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))
