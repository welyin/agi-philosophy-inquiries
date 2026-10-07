"""950: finite-scale cost of the literal 942/947 Dirac-flavor realization.

Reuses the subtracted free-fermion spectrum from round 630. One closed
fermion loop, flat vacuum, spacelike two-point response only. No all-order
control or parent-to-protocol matching is inferred from small residuals.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'material_vacuum_matching_results.json'

def read(p): return json.loads(p.read_text('utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a): return float(np.linalg.norm(a, 2))

def remainder(masses, multiplicity, coupling, q, order):
    # v = sqrt(1 - 4 M^2/s). Integrate the positive residual and its
    # deficit from the leading moment separately, avoiding cancellation.
    x, w = np.polynomial.legendre.leggauss(order)
    v, w = (x + 1)/2, w/2
    t = q*q/(4*masses[:, None]**2)
    denominator = 1 + t*(1-v*v)
    integral = (v**4/denominator) @ w
    deficit_integral = (v**4*(1-v*v)/denominator) @ w
    residual = multiplicity*coupling**2*q**4/(16*np.pi**2)*np.sum(integral/masses**2)
    deficit = multiplicity*coupling**2*q**6/(64*np.pi**2)*np.sum(deficit_integral/masses**4)
    return float(residual), float(deficit)

def run():
    old = read(STAGE/'947/protocol_field_transport_results.json')
    bridge = read(STAGE/'948/neutral_matter_bridge_results.json')
    weak = read(STAGE/'949/joint_weak_domain_results.json')
    assert all(r['all_scientific_checks_passed'] for r in (old, bridge, weak))
    spec = importlib.util.spec_from_file_location('protocol929', STAGE/'929/joint_protocol_selection.py')
    protocol = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(protocol)
    # These gates only act on qubits 0 and 1. Reuse their actual action;
    # the other eleven work and four fault qubits are exact degeneracies.
    core = np.eye(protocol.DIM, 4, dtype=complex)
    for gate in protocol.circuit(np.pi/2, (2, 3))[:2]:
        core = protocol.act(core, gate)
    assert np.linalg.norm(core[4:]) == 0
    unitary = core[:4]
    z = np.diag([1., 1., -1., -1.])
    c = unitary.conj().T @ z @ unitary
    b_eigen_error = float(np.max(abs(np.linalg.eigvalsh(c)-[-1,-1,1,1])))
    assert b_eigen_error < 1e-12 and abs(np.trace(c)) < 1e-12
    n = 30
    off = .5*np.sqrt(np.arange(1, n+1)*np.arange(n, 0, -1))
    scale = old['same_physical_protocol']['internal_clock_scale']
    hc = scale*(15*np.eye(31)+np.diag(off,1)+np.diag(off,-1))
    eig, vec = np.linalg.eigh(hc)
    exact = scale*np.arange(31)
    spectrum_error = float(np.max(abs(eig-exact)))
    assert spectrum_error < 1e-12
    degeneracy = 2**17
    flavors = degeneracy*len(exact)
    assert flavors == old['same_physical_protocol']['original_work_and_fault_and_clock_dimension']
    binomial = np.array([math.comb(30,e)/2**30 for e in range(31)])
    weight_error = float(np.max(abs(abs(vec[0])**2-binomial)))
    assert weight_error < 1e-12
    mass = 1e14
    masses = mass+exact
    coupling = bridge['parameters']['g_source']
    K = np.array(bridge['parameters']['inherited_K'])
    cw = bridge['parameters']['W_mass_squared_derivative']
    cf = bridge['parameters']['fermion_mass_derivative']
    C = np.array([[coupling,0.,0.],[0.,cw,cf]])
    L4 = float(degeneracy*coupling**2/(80*np.pi**2)*np.sum(1/masses**2))
    L6 = float(degeneracy*coupling**2/(1120*np.pi**2)*np.sum(1/masses**4))
    wrong_state_weight_L4 = float(coupling**2/(80*np.pi**2)*np.sum(binomial/masses**2))
    assert abs(L4/wrong_state_weight_L4/flavors-1) < 1e-12
    # Independent Dirac projector trace fixes the normalization in 630,
    # with scalar Yukawa g replacing that round's mass-source m.
    pauli = [np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    beta = np.diag([1.,1.,-1.,-1.])
    alpha = [np.block([[np.zeros((2,2)),s],[s,np.zeros((2,2))]]) for s in pauli]
    spin_error = 0.
    for k in (np.array([.2,.3,.4]), np.array([1.,-.2,.7])):
        mtest = 1.7
        energy = math.sqrt(mtest*mtest+k@k)
        hdirac = beta*mtest+sum(a*x for a,x in zip(alpha,k))
        plus, minus = (np.eye(4)+hdirac/energy)/2, (np.eye(4)-hdirac/energy)/2
        measured = np.trace(minus@(coupling*beta)@plus@(coupling*beta)).real
        spin_error = max(spin_error, abs(measured-2*coupling**2*(k@k)/energy**2))
    assert spin_error < 1e-13
    quadrature_error = 0.
    rows = []
    for q in (.25,1.,4.):
        a,b = remainder(masses,degeneracy,coupling,q,32)
        aa,bb = remainder(masses,degeneracy,coupling,q,64)
        quadrature_error = max(quadrature_error,abs(a-aa)/aa,abs(b-bb)/bb)
        upper = L4*q**4
        deficit_upper = L6*q**6
        assert aa <= upper*(1+1e-12) and bb <= deficit_upper*(1+1e-12)
        assert abs(aa+bb-upper)/upper < 1e-12
        D = np.linalg.inv(q*q*np.eye(2)+K)
        eta = upper*D[0,0]
        rows.append(dict(spacelike_Q=q,one_loop_subtracted_residual=aa,
                         leading_positive_upper=upper,leading_deficit=bb,
                         deficit_upper=deficit_upper,leading_relative_deficit_upper=deficit_upper/upper,
                         shared_response_relative_change_upper=eta/(1-eta)))
    # A resolvable diagnostic checks the integrals beyond the heavy-mass
    # regime; these masses are NOT used as the physical material.
    test_mass = np.array([2.,3.])
    test_q = 3.
    a,b = remainder(test_mass,3,.4,test_q,32)
    aa,bb = remainder(test_mass,3,.4,test_q,64)
    diagnostic_lead = 3*.4**2*test_q**4/(80*np.pi**2)*np.sum(1/test_mass**2)
    assert abs(aa+bb-diagnostic_lead) < 1e-14
    assert bb > 1e-4 and abs(a-aa)+abs(b-bb) < 1e-14
    assert quadrature_error < 1e-12
    # Full Q in [0,4] bound, using D_pp <= 1/lambda_min(K).
    window_kernel_bound = L4*4**4
    window_eta = window_kernel_bound/np.linalg.eigvalsh(K)[0]
    window_response_bound = window_eta/(1-window_eta)
    # Different local matching produces a finite difference even while
    # the SAME spectral residual is tiny. Once measured K is held fixed,
    # this is no longer a renormalization-scheme freedom.
    A = np.eye(2)+K
    delta = .1
    shifted = A + np.diag([delta,0.])
    D0,D1 = np.linalg.inv(A),np.linalg.inv(shifted)
    R0,R1 = -C.T@D0@C,-C.T@D1@C
    local_change = abs(R1[0,1]/R0[0,1]-1)
    local_formula = delta*D0[0,0]/(1+delta*D0[0,0])
    assert abs(local_change-local_formula) < 1e-14 and local_change > .04
    assert norm(R1-R1.T) < 1e-14
    assert abs(R1[0,1]/cw-R1[0,2]/cf) < 1e-14
    # The rank-one identity is checked at a resolvable kernel shift;
    # the physical 1e-25 change cannot be obtained by floating subtraction.
    sherman = D0-delta*np.outer(D0[:,0],D0[0,:])/(1+delta*D0[0,0])
    sherman_error = norm(sherman-D1)
    assert sherman_error < 1e-14
    parameters = []
    for epsilon in (1.,.001):
        collective = flavors*(epsilon*coupling)**2/(16*np.pi**2)
        parameters.append(dict(epsilon=epsilon,
            collective_loop_factor_diagnostic=collective,
            one_loop_Q1_residual_upper=epsilon**2*L4,
            dimension_two_local_scale_diagnostic=collective*mass**2))
    assert parameters[0]['collective_loop_factor_diagnostic'] > 6000
    assert parameters[1]['collective_loop_factor_diagnostic'] < .007
    files = [Path(__file__),STAGE/'929/joint_protocol_selection.py',
             STAGE/'947/protocol_field_transport_results.json',
             STAGE/'948/neutral_matter_bridge_results.json',
             STAGE/'949/joint_weak_domain_results.json',
             STAGE.parent/'archive_629_652/research_note_630.md']
    return dict(round=950,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(literal_Dirac_flavors=flavors,clock_eigenvalues=exact.tolist(),
            each_energy_multiplicity=degeneracy,each_energy_each_B_sign_multiplicity=degeneracy//2,
            base_mass=mass,yukawa_coupling=coupling,K=K.tolist(),spacelike_Q_window=[0.,4.]),
        checks=dict(actual_B_core_eigenvalue_error=b_eigen_error,actual_clock_spectrum_error=spectrum_error,
            state_binomial_weight_error=weight_error,Dirac_projector_trace_error=spin_error,
            positive_quadrature_relative_agreement=quadrature_error,Sherman_Morrison_error=sherman_error,
            finite_local_matching_response_identity_error=abs(local_change-local_formula)),
        one_loop=dict(L4=L4,L6=L6,wrong_state_probability_weight_L4=wrong_state_weight_L4,
            true_to_wrong_species_count_ratio=L4/wrong_state_weight_L4,rows=rows,
            entire_Q_window_kernel_upper=window_kernel_bound,
            entire_Q_window_common_response_relative_upper=float(window_response_bound),
            odd_p_vertices_cancel_in_pure_fermion_determinant=True),
        local_matching_counterexample=dict(delta_K_pp=delta,Q=1.,
            record_W_before=float(R0[0,1]),record_W_after=float(R1[0,1]),
            record_fermion_before=float(R0[0,2]),record_fermion_after=float(R1[0,2]),
            common_relative_response_change=float(local_change)),
        weak_family_diagnostics=parameters,
        scope=dict(one_closed_fermion_loop_only=True,subtraction_constants_are_matching_inputs=True,
            local_mass_scale_is_not_a_physical_prediction=True,all_orders_certified=False,
            full_metric_source_or_parent_matching_certified=False,
            same_state_weights_as_vacuum_multiplicities=False,
            composite_internal_levels_identified_as_fundamental_species=False,
            new_finite_protocol_probability_bound=False,physical_cutoff_or_UV_completion_inferred=False,
            literal_flavor_candidate_selected_as_final=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        # Small physical effects must not pass because of an absolute 1e-10 tolerance.
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=0.),(a,b)
    else: assert a==b,(a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else: compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','parameters')},ensure_ascii=False,indent=2))
