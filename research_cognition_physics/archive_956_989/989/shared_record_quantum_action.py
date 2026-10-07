"""989: classical records can select, but cannot replace, a quantum contact.

Reuses the 958 material, time, window and all-input bound. No parameter search.
LOCC and Schmidt-overlap facts are inherited, not claimed as new theorems.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
OUT=HERE/'shared_record_quantum_action_results.json'

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def pack(x):return dict(exact=str(x),value=float(x))

def run():
    saved=json.loads((STAGE/'958/capacitive_material_write_results.json').read_text('utf-8'))
    interface=saved['exact_interface']
    material=load('material968',STAGE/'968/internal_relay.py')
    m,h,q,_,W,material_checks=material.material()
    phase=load('phase965',STAGE/'965/material_field_window.py').phases
    J=m['J'];T=interface['controlled_phase_time'];window=interface['read_window_half_width']
    chi=interface['conditional_energy_chi'];I=np.eye(4)
    Hfree=np.kron(h,I)+np.kron(I,h)
    contact=.2*np.kron(q,q)
    H=Hfree+contact
    # The same V0=.3 common density energy in 958; fixed in this task.
    common_energy=1.2
    physical_H=H+common_energy*np.eye(16)
    code=np.kron(W,W)
    plus=np.ones(2)/math.sqrt(2)
    initial=code@np.kron(plus,plus)
    rot=np.diag([phase(np.array([-J]),T)[0],1.])
    local=W@rot
    target=np.kron(local,local)@np.array([-1.,1.,1.,1.])/2
    projector=np.outer(target,target.conj())
    P=local@local.conj().T
    X=local@np.array([[0,1],[1,0]])@local.conj().T
    Y=local@np.array([[0,-1j],[1j,0]])@local.conj().T
    Z=local@np.diag([1.,-1.])@local.conj().T
    effects=[np.kron(P,P),np.kron(X,Z),np.kron(Z,X),np.kron(Y,Y)]
    signs=[1,-1,-1,1]
    witness=sum(s*A for s,A in zip(signs,effects))/4
    decomposition_error=norm(witness-projector)
    assert decomposition_error<1e-14
    schmidt=np.linalg.svd(target.reshape(4,4),compute_uv=False)**2
    assert np.max(abs(schmidt-np.array([.5,.5,0,0])))<1e-14
    # Independent rational proof that eta=(sqrt(29)-5)/(sqrt(29)*(sqrt(29)+3))<.00853.
    assert 29*97441**2<524737**2
    eta_upper=F(853,100000)
    state_error=2*eta_upper
    drift_upper=F(1,50) # (2J+chi)*.01/(J+chi)<.02.
    lower=1-state_error-drift_upper
    classical_limit=F(1,2)
    mixture_lower=(lower+F(1,4))/2
    assert float(state_error)>interface['all_time_unknown_input_error_upper']
    assert (2*J+chi)*window<float(drift_upper)
    assert lower==F(48147,50000) and mixture_lower==F(60647,100000)

    vals,V=np.linalg.eigh(H)
    vals0,V0=np.linalg.eigh(Hfree)
    certificates=[]
    rows=[]
    for t in (T-window,T,T+window):
        U=(V*phase(vals,t))@V.conj().T
        U0=(V0*phase(vals0,t))@V0.conj().T
        psi=U@initial;free=U0@initial
        expectation=lambda vec,A:float(np.vdot(vec,A@vec).real)
        correlators=[expectation(psi,A) for A in effects]
        fidelity=abs(np.vdot(target,psi))**2
        free_fidelity=abs(np.vdot(target,free))**2
        free_exact=(1+math.sin(J*(t-T))**2)/4
        # A fully classical, conserved 0/1 label, with equal weights.
        rho=.5*np.outer(psi,psi.conj())+.5*np.outer(free,free.conj())
        mixed_fidelity=float(np.trace(projector@rho).real)
        eig_res=norm(H@V-V*vals)
        gram=norm(V.conj().T@V-np.eye(16))
        numerical_budget=gram+abs(t)*(eig_res+1e-13)*norm(V)+1e-8
        assert numerical_budget<1e-6
        assert fidelity>float(lower) and mixed_fidelity>float(mixture_lower)
        assert abs(free_fidelity-free_exact)<1e-9
        assert abs(sum(s*x for s,x in zip(signs,correlators))/4-fidelity)<1e-13
        drift=expectation(psi,physical_H)-expectation(initial,physical_H)
        assert abs(drift)<1e-12
        rows.append(dict(time=float(t),four_local_statistics=correlators,
            native_fidelity=float(fidelity),free_fidelity=float(free_fidelity),
            free_exact_formula=float(free_exact),unread_classical_flag_fidelity=mixed_fidelity,
            energy_drift=drift,contact_source=expectation(psi,np.kron(q,q))))
        certificates.append(dict(eigen_residual=eig_res,gram_defect=gram,
                                 operator_arithmetic_budget=numerical_budget))

    # Explicit block-controlled H: the classical flag is not the quantum mediator.
    Hflag=np.zeros((32,32));Hflag[:16,:16]=Hfree+common_energy*np.eye(16)
    Hflag[16:,16:]=physical_H
    flag=np.diag([0.]*16+[1.]*16)
    flag_commutator=norm(Hflag@flag-flag@Hflag)
    assert flag_commutator==0 and np.linalg.eigvalsh(Hflag)[0]>1.
    # A separable state reaches the 1/2 ceiling; this is not an LOCC optimization.
    a=local[:,0];b=target.reshape(4,4).T@a.conj();b=b/np.linalg.norm(b)
    product=np.kron(a,b)
    separable_fidelity=float(abs(np.vdot(target,product))**2)
    assert abs(separable_fidelity-.5)<1e-14
    # Fully dephase local logical inputs before the native interaction.
    # It is not asserted that every conceivable classical channel is dephasing.
    center=T
    U=(V*phase(vals,center))@V.conj().T
    diagonal=.25*(U@code)@(U@code).conj().T
    dephased_score=float(np.trace(projector@diagonal).real)
    assert dephased_score<.26

    sources=[STAGE/'research_note_958.md',STAGE/'958/capacitive_material_write.py',
        STAGE/'958/capacitive_material_write_results.json',STAGE/'968/internal_relay.py',
        STAGE/'research_note_929.md',STAGE/'research_note_987.md',STAGE/'research_note_988.md',
        STAGE/'957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE.parent/'archive_217_222/research_note_220.md',
        STAGE.parent/'archive_259_300/research_note_269.md',
        HERE/'drafts/STATUS.md',HERE/'drafts/mechanism_adoption_decision.md']
    return dict(round=989,all_scientific_checks_passed=True,
        inputs=dict(U=1,v=.1,kappa=.2,common_density_energy=common_energy,
            time_from_958=T,half_window_from_958=window,chi_from_958=chi,
            initial_local_codes=['plus','plus'],shared_flag_one_probability=.5),
        material_restriction_checks=material_checks,
        analytic=dict(eta_upper=pack(eta_upper),native_state_error_upper=pack(state_error),
            fixed_effect_window_drift_upper=pack(drift_upper),
            native_witness_lower=pack(lower),all_LOCC_separable_ceiling=pack(classical_limit),
            finite_prediction_gap_lower=pack(lower-classical_limit),
            unread_classical_flag_witness_lower=pack(mixture_lower),
            overlap_proof_integer_margin=524737**2-29*97441**2),
        target_schmidt_squares=schmidt.tolist(),local_decomposition_error=decomposition_error,
        rows=rows,numerical_diagnostics=certificates,
        flag_commutator_norm=flag_commutator,
        explicit_product_state_witness=separable_fidelity,
        fully_locally_dephased_code_score=dephased_score,
        decision=dict(evidence_policy_is_optional_internal_capability=True,
            universal_prediction_optimization_adopted=False,
            replace_all_quantum_contacts_by_classical_messages_rejected=True,
            every_relation_label_must_be_coherent=False,
            classical_label_selecting_quantum_contact_is_allowed=True,
            additional_shared_entanglement_excluded_only_in_classical_comparator=True,
            full_goal_completed=False,gravitational_quantization_proved=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sources},
        references=['https://arxiv.org/abs/quant-ph/9911117','https://arxiv.org/abs/1408.4740'])

def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'/{i}')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-8),(path,a,b)
    else:assert a==b,(path,a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(OUT.read_text('utf-8')))
    print(json.dumps({k:result[k] for k in ('round','all_scientific_checks_passed','analytic','rows')},
        ensure_ascii=False,indent=2))
