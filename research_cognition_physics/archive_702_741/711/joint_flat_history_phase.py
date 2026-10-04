"""711: flat scalar history phases on the original finite Gauss process.

Uses the original quotient group, tree dictionary, 32 CAR representation and
non-diagonal electric coefficients. Tests exact finite-model identities only;
classification, operator domains and heat traces are proved in the note.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_tree_matter_source as tree
import joint_charge_changing_vertex as vertex

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_flat_history_phase_results.json'


def character(g):
    return np.linalg.det(g[2]**3*g[1])


def primitive_loop(t):
    C=np.diag(np.exp(2j*np.pi*t*np.array([1,1,-2])/3))
    W=np.diag(np.exp(1j*np.pi*t*np.array([1,-1])))
    return C,W,np.exp(1j*np.pi*t/3)


def integer_cycles(edges,vertices,tree_count):
    D=np.zeros((len(edges),vertices),int)
    for e,(s,t) in enumerate(edges):D[e,s]=1;D[e,t]=-1
    non_tree=len(edges)-tree_count
    C=np.zeros((len(edges),non_tree),int);C[tree_count:]=np.eye(non_tree,dtype=int)
    tree_block=D[:tree_count,1:].T
    assert tree_block.shape==(vertices-1,vertices-1)
    assert abs(round(np.linalg.det(tree_block)))==1
    C[:tree_count]=np.rint(np.linalg.solve(tree_block,-D[tree_count:,1:].T)).astype(int)
    assert np.array_equal(D.T@C,np.zeros((vertices,non_tree),int))
    return D,C


def topology_normalization():
    ts=np.linspace(0,1,65)
    chis=np.array([character(primitive_loop(t)) for t in ts])
    winding=float((np.unwrap(np.angle(chis))[-1]-np.unwrap(np.angle(chis))[0])/(2*np.pi))
    R=tree.matter.representation(*primitive_loop(1))
    endpoint=float(np.linalg.norm(R-np.eye(32)))
    faithful_dictionary=0.
    for t in ts:
        g=primitive_loop(t);r=tree.matter.representation(*g)
        faithful_dictionary=max(faithful_dictionary,float(abs(character(g)-1/r[28,28])))
    cover=np.exp(12j*np.pi*ts)
    cover_winding=float(np.unwrap(np.angle(cover))[-1]/(2*np.pi))
    weak=np.array([character((np.eye(3),np.diag(np.exp(2j*np.pi*t*np.array([1,-1]))),1)) for t in ts])
    weak_winding=float(np.unwrap(np.angle(weak))[-1]/(2*np.pi))
    assert abs(winding-1)<1e-14 and endpoint<1e-13 and faithful_dictionary<1e-13
    assert abs(cover_winding-6)<1e-14 and abs(weak_winding)<1e-14
    rows=[]
    for edges in (tree.ENDS,tree.ENDS+((0,3),)):
        D,C=integer_cycles(edges,4,3)
        node_windings=np.array([2,-1,3,0],int)
        period_pairing=C.T@D@node_windings
        assert np.array_equal(period_pairing,np.zeros(len(edges)-3,int))
        rows.append(dict(edges=len(edges),nodes=4,configuration_cycle_rank=len(edges)-3,
            incidence=D.tolist(),integral_cycle_basis=C.tolist(),gauge_orbit_winding_pairing=period_pairing.tolist()))
    return dict(primitive_quotient_loop_winding=winding,original32_endpoint_identity_error=endpoint,
        primitive_character_dictionary_error=faithful_dictionary,cover_U1_full_turn_winding=cover_winding,
        weak_SU2_loop_winding=weak_winding,graph_rows=rows,
        no_claim_graph_cycle_rank_equals_spatial_manifold_betti_number=True)


def loop_character(edges):
    _,_,loop=tree.rooted(edges,np.zeros((4,5)))
    return character(loop)


def gauss_and_actual_gradients():
    rng=np.random.default_rng(7112);D,C=integer_cycles(tree.ENDS,4,3);c=C[:,0]
    gauge_error=lift_error=gradient_error=0.
    step=1e-6
    for _ in range(5):
        edges=[tree.old.sample(rng) for e in tree.ENDS]
        gs=[tree.old.sample(rng) for _ in range(4)]
        transformed=[tree.mul(gs[s],u,tree.old.inverse(gs[t])) for u,(s,t) in zip(edges,tree.ENDS)]
        base=loop_character(edges)
        gauge_error=max(gauge_error,float(abs(loop_character(transformed)-base)))
        direct=np.prod([character(u)**int(k) for u,k in zip(edges,c)])
        assert abs(direct-base)<1e-13
        for e in range(4):
            for sector in (0,1,2):
                plus=list(edges);minus=list(edges)
                plus[e]=tree.mul(tree.one_parameter(sector,0,step),edges[e])
                minus[e]=tree.mul(tree.one_parameter(sector,0,-step),edges[e])
                got=-1j*(loop_character(plus)-loop_character(minus))/(2*step)
                target=6*c[e]*base if sector==2 else 0.
                gradient_error=max(gradient_error,float(abs(got-target)))
        zeta=primitive_loop(1)
        for e in range(4):
            lifted=list(edges);lifted[e]=tree.mul(zeta,edges[e])
            lift_error=max(lift_error,float(abs(loop_character(lifted)-base)))
    assert max(gauge_error,lift_error)<1e-13 and gradient_error<2e-9
    theta=.71
    return dict(orientation_c=c.tolist(),Gauss_invariance_error=gauge_error,
        independent_central_lift_error=lift_error,original_electric_generator_FD_error=gradient_error,
        primitive_chord_holonomy=[float(np.cos(theta)),float(np.sin(theta))],
        weak_loop_holonomy=[1.,0.],cycle_form_is_horizontal=np.array_equal(D.T@c,np.zeros(4,int)))


FLOW=np.array([[-1,1,0],[0,-1,0],[1,0,0]],float)


def electric_coefficient(t=.8,sigma=.12):
    K=tree.coefficients(t,sigma)[2]
    return float(np.einsum('vi,ij,vj->',FLOW,K,FLOW))


def physical_sources():
    zeta=.37;m=1;h=2e-5
    k=electric_coefficient()
    # Difference of complete H expectations in two legal states. All other
    # terms cancel analytically; neither state is claimed to be an eigenstate.
    E=lambda z,s=.12:36*(m-z)**2*electric_coefficient(sigma=s)
    gap=E(zeta)-36*zeta**2*k
    source=-72*(m-zeta)*k;second=72*k
    fd=(E(zeta+h)-E(zeta-h))/(2*h)
    fdd=(E(zeta+h)-2*E(zeta)+E(zeta-h))/h**2
    geom=(E(zeta,.12+h)-E(zeta,.12-h))/(2*h)
    shear_expected=36*(m-zeta)**2*float(np.einsum('vi,ij,vj->',FLOW,tree.coefficients()[2]@tree.SHAPE,FLOW))
    shear_fd=36*(m-zeta)**2*(electric_coefficient(.8+h)-electric_coefficient(.8-h))/(2*h)
    rel=lambda a,b:float(abs(a-b)/max(1.,abs(b)))
    errors=dict(phase=rel(fd,source),phase_second=rel(fdd,second),conformal=rel(geom,-2*E(zeta)),shear=rel(shear_fd,shear_expected))
    assert max(errors.values())<1e-6 and k>0
    periodic=abs(36*(m+1-(zeta+1))**2*k-E(zeta));assert periodic<1e-13
    # Actual original coefficient test: adding scalar bosonic derivatives
    # cannot modify any original CAR quark-charge commutator.
    q=np.diag(np.r_[np.ones(24),np.zeros(8)])
    rng=np.random.default_rng(7113);charge_error=0.
    for _ in range(5):
        hm,dm=tree.matter.mass_matrices(.24*rng.normal(size=5))
        charge_error=max(charge_error,float(np.linalg.norm(q@hm-hm@q)),float(np.linalg.norm(q@dm+dm@q.T)))
    assert charge_error<1e-13
    # Pure QQQL charge still permits exactly the710 phase conjugation; use
    # original sparse CAR, not a two-level surrogate Hamiltonian.
    state=vertex.old.normalized(vertex.old.add({0:1.},vertex.operation({0:1.},True),1j))
    alpha=.23;g=.19*np.exp(.41j)
    got=vertex.phase(vertex.vertex(vertex.phase(state,alpha),g),-alpha)
    expected=vertex.vertex(state,g*np.exp(3j*alpha))
    covariance=vertex.distance(got,expected);assert covariance<1e-13
    return dict(original_non_diagonal_K_U1=tree.coefficients()[2].tolist(),hbar=1.,
        normalized_phase=zeta,integer_character=m,electric_cycle_coefficient=k,
        actual_electric_expectation=E(zeta),full_H_gap_to_m_zero=gap,
        original_conformal_source_gap=-2*gap,phase_source=source,phase_source_contact=second,
        shear_source=shear_expected,finite_difference_relative_errors=errors,
        joint_integer_phase_state_shift_error=periodic,original_quark_charge_error=charge_error,
        QQQL_single_phase_covariance_error=covariance,
        flat_phase_not_second_baryon_changing_amplitude=True,
        finite_energy_physical_states_not_full_H_eigenstates=True,
        thermal_existence_by_full_operator_proof_not_truncated_spectrum=True)


def run():
    names=('joint_tree_matter_source.py','joint_tree_matter_source_results.json',
        'joint_charge_changing_vertex.py','joint_charge_changing_vertex_results.json',
        'research_note_623.md','research_note_629.md','round711_drafts/topology_sampling_entry_results.json')
    return dict(round=711,tests_run=3,failures=0,errors=0,topology=topology_normalization(),
        Gauss_and_generators=gauss_and_actual_gradients(),full_process_sources=physical_sources(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope='All smooth scalar flat phase modifications on original unreduced/based finite configuration, retaining original Gauss lift. Abelian graph-cycle classes and actual full-H source matching; no weak c2 angle or anomaly generation, no general continuum or nonflat no-go.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
