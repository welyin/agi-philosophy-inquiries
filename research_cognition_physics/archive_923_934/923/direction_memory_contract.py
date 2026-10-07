"""923: one direction-report contract tested with opposite propagation sectors.
Finite operational model, not a derivation of Dirac dynamics, GR or field content.
"""
from pathlib import Path
import json,argparse,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
TARGET=HERE/'direction_memory_contract_results.json'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
P=(I,X,Y,Z);ALPHA=[np.kron(Z,a) for a in P[1:]];BETA=np.kron(X,I)
def mx(a):return float(np.max(abs(a)))
def channel(rho,nu):return (np.trace(rho)*I+nu*sum(np.trace(rho@a)*s for a,s in zip(ALPHA,P[1:])))/2

def choi(nu):
    blocks=[]
    for i in range(4):
        row=[]
        for j in range(4):
            e=np.zeros((4,4),complex);e[i,j]=1;row.append(channel(e,nu))
        blocks.append(row)
    return np.block(blocks)

def kraus(nu):
    out=[]
    for branch,lam in ((0,nu),(1,-nu)):
        E=np.zeros((2,4),complex);E[:,2*branch:2*branch+2]=I
        weights=[(1+3*lam)/4]+[(1-lam)/4]*3
        for w,s in zip(weights,P):
            assert w>=-1e-15
            if w>1e-14:out.append(np.sqrt(w)*s@E)
    return out

def evolve(H,t):
    e,V=np.linalg.eigh(H);return (V*np.exp(-1j*t*e))@V.conj().T

def algebra(gens):
    basis=[]
    def add(a):
        v=a.ravel().astype(complex).copy()
        for q in basis:v-=np.vdot(q,v)*q
        n=np.linalg.norm(v)
        if n>1e-10:basis.append(v/n);return True
        return False
    add(np.eye(4))
    for a in gens:add(a)
    changed=True
    while changed:
        changed=False
        for q in list(basis):
            for a in gens:changed=add(q.reshape(4,4)@a) or changed
    return len(basis)

def run():
    rows=[]
    for nu in (0.,1/3,.5,1.):
        actual=np.linalg.eigvalsh(choi(nu))
        predicted=np.sort([(1+3*nu)/2]+[(1-nu)/2]*3+[(1-3*nu)/2]+[(1+nu)/2]*3)
        err=mx(actual-predicted);assert err<1e-14
        rows.append(dict(visibility=nu,choi_eigenvalues=actual.tolist(),analytic_spectrum_error=err))
    assert np.linalg.eigvalsh(choi(1.))[0]<-.99
    nu=1/3;Ks=kraus(nu);complete=sum(k.conj().T@k for k in Ks);assert mx(complete-np.eye(4))<1e-14
    worst=0.;rng=np.random.default_rng(923)
    for _ in range(12):
        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=a@a.conj().T;rho/=np.trace(rho)
        worst=max(worst,mx(channel(rho,nu)-sum(k@rho@k.conj().T for k in Ks)))
    assert worst<1e-14
    # The two input states have identical report AND identical chirality-dephased state.
    spin=(I+Z)/2;states=[np.kron((I+s*Y)/2,spin) for s in (1,-1)]
    mass=1.;H=mass*BETA;t=np.pi/(4*mass);U=evolve(H,t)
    states_t=[U@q@U.conj().T for q in states]
    velocities0=[[float(np.trace(q@a).real) for a in ALPHA] for q in states]
    velocitiest=[[float(np.trace(q@a).real) for a in ALPHA] for q in states_t]
    probability=[float(np.trace(q@(np.eye(4)+ALPHA[2])/2).real) for q in states_t]
    dephase=lambda q:sum(E@q@E for E in (np.diag([1,1,0,0]),np.diag([0,0,1,1])))
    assert mx(dephase(states[0])-dephase(states[1]))==0
    assert mx(channel(states[0],nu)-channel(states[1],nu))==0
    assert abs(probability[0]-probability[1]-1)<1e-14
    future_report_distance=float(np.linalg.norm(channel(states_t[0],nu)-channel(states_t[1],nu),ord='nuc')/2)
    assert abs(future_report_distance-nu)<1e-14
    dims=dict(instantaneous_closed_direction_algebra=algebra(ALPHA),mass_mixed_future_algebra=algebra(ALPHA+[U.conj().T@a@U for a in ALPHA]))
    assert list(dims.values())==[8,16]
    # Pure dilation: report2 x internal environment8; seven nonzero Kraus slots.
    V=np.zeros((2,8,4),complex)
    for j,k in enumerate(Ks):V[:,j,:]=k
    V=V.reshape(16,4);Henc=V@H@V.conj().T;Uenc=evolve(Henc,t)
    iso=mx(V.conj().T@V-np.eye(4));intertwine=mx(Uenc@V-V@U)
    future=[];partial_errors=[]
    for rho in states:
        full=V@rho@V.conj().T
        report=np.einsum('aebe->ab',full.reshape(2,8,2,8))
        partial_errors.append(mx(report-channel(rho,nu)))
        full_t=Uenc@full@Uenc.conj().T
        effect=V@(np.eye(4)+ALPHA[2])@V.conj().T/2
        future.append(float(np.trace(full_t@effect).real))
    assert max(iso,intertwine,*partial_errors)<1e-14
    assert mx(np.array(future)-probability)<1e-14
    # The declared paired principal symbol has a relativistic square; this is an input model check.
    p=np.array([.2,-.3,.4]);cp=sum(a*k for a,k in zip(ALPHA,p))+H
    square=mx(cp@cp-(np.dot(p,p)+mass*mass)*np.eye(4));assert square<1e-14
    return dict(round=923,date='2026-10-06',all_scientific_checks_passed=True,choi_rows=rows,
      cp_range_for_nonnegative_visibility=[0.,1/3],sharp_report_cp=False,sharp_normalized_reference_witness_minimum=-.25,
      legal_boundary_report_kraus_rank=len(Ks),kraus_completeness_error=mx(complete-np.eye(4)),kraus_channel_identity_error=worst,
      initial_velocities=velocities0,final_velocities=velocitiest,final_original_direction_probabilities=probability,
      initial_report_equal=True,initial_classical_chirality_summary_equal=True,future_report_trace_distance=future_report_distance,
      original_future_probability_prediction_minimax_error_lower_bound=.5,
      finite_task_algebra_dimensions=dims,
      dilation_input_dimension=4,dilation_report_dimension=2,dilation_environment_dimension=8,nonzero_environment_dimension=len(Ks),
      isometry_error=iso,same_Hamiltonian_intertwining_error=intertwine,report_partial_trace_errors=partial_errors,
      retained_global_memory_future_probabilities=future,paired_symbol_square_identity_error=square,
      positive_pair_model_and_registers_are_explicit_inputs=True,Dirac_dynamics_derived_from_cognition=False,
      apparatus_autonomous_implementation_proved=False,energy_cost_of_preparing_dilation_computed=False,
      whole_stage_completed=False,full_goal_completed=False,
      source_hashes={str(Path(__file__).relative_to(STAGE)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run()
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(d,ensure_ascii=False,indent=2))
