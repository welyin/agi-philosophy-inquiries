"""One bounded test of program-source versus actual-rate-source identity.

Uses the established gate-dressed history clock.  No gravity is simulated.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'autonomous_source_identity_results.json'
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
K=np.diag([0.,.5]).astype(complex)
HC=(I+X)/2
T=np.pi
W=np.vstack([I,np.zeros((2,2),complex)])
OUT=np.vstack([np.zeros((2,2),complex),I])
PLUS=np.array([1,1],complex)/np.sqrt(2)
RHO=np.outer(PLUS,PLUS.conj())
E=(I+X)/2

def norm(a):return float(np.linalg.norm(a))
def exph(h,t):
    v,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*v))@u.conj().T
def gate(u):return exph(K,T*(1+u))
def dressed(u):
    g=gate(u)
    return .5*np.block([[I,g.conj().T],[g,I]])
def emb(source,u):return exph(source(u),T)@W
def prob(v):return float(np.trace(np.kron(I,E)@v@RHO@v.conj().T).real)
def gram(source,um,up):return emb(source,um).conj().T@emb(source,up)
def clock_cf(du):return np.exp(-1j*T*du/2)*np.cos(T*du/2)

def run():
    h0=dressed(0)
    hd=np.kron(HC,I)+np.kron(I,K)
    fam={'program':dressed,'actual_rate':lambda u:(1+u)*h0,'direct_rate':lambda u:(1+u)*hd}
    endpoint={name:norm(emb(f,0)+OUT@gate(0)) for name,f in fam.items()}
    assert max(endpoint.values())<2e-14
    moments=[]
    for n in range(1,5):
        a=W.conj().T@np.linalg.matrix_power(h0,n)@W
        b=W.conj().T@np.linalg.matrix_power(hd,n)@W
        expected_direct=(np.linalg.matrix_power(K,n)+np.linalg.matrix_power(K+I,n))/2
        assert norm(a-.5*I)<1e-13 and norm(b-expected_direct)<1e-13
        moments.append(dict(order=n,compiled_diagonal=np.diag(a).real.tolist(),direct_diagonal=np.diag(b).real.tolist()))
    rows=[]
    for u in (-.1,-.02,0.,.02,.1):
        target=(1+np.cos(T*(1+u)/2))/2
        actual=.5+.5*np.cos(T*(1+u)/2)**2
        ps={name:prob(emb(f,u)) for name,f in fam.items()}
        assert abs(ps['program']-target)<1e-13 and abs(ps['direct_rate']-target)<1e-13
        assert abs(ps['actual_rate']-actual)<1e-13
        rows.append(dict(source=u,**ps,target_probability=float(target)))
    grams=[]
    for um,up in ((0.,.07),(-.03,.04),(.02,-.08)):
        du=up-um;logical=exph(K,T*du);cf=clock_cf(du)
        a=gram(fam['program'],um,up);b=gram(fam['actual_rate'],um,up);c=gram(fam['direct_rate'],um,up)
        errs=[norm(a-logical),norm(b-cf*I),norm(c-cf*logical)]
        assert max(errs)<2e-13
        grams.append(dict(left_source=um,right_source=up,identity_errors=errs,
                          program_actual_rate_difference=norm(a-b)))
    deriv=[]
    for eps in (1e-3,5e-4,2.5e-4):
        ds={name:(prob(emb(f,eps))-prob(emb(f,-eps)))/(2*eps) for name,f in fam.items()}
        deriv.append(dict(step=eps,**ds,program_error=abs(ds['program']+np.pi/4)))
        assert abs(ds['actual_rate'])<1e-10
        assert abs(ds['program']+np.pi/4)<1e-6 and abs(ds['direct_rate']+np.pi/4)<1e-6
    # Analytic first parameter derivative of the dressed off-diagonal gate.
    g=gate(0);dg=-1j*T*K@g
    dp=.5*np.block([[np.zeros((2,2)),dg.conj().T],[dg,np.zeros((2,2))]])
    program_initial_force=W.conj().T@dp@W
    actual_initial_force=W.conj().T@h0@W
    assert norm(program_initial_force)<1e-14 and norm(actual_initial_force-.5*I)<1e-14
    energy=[]
    for name,h in [('compiled',h0),('direct',hd)]:
        a=W.conj().T@h@W;v=exph(h,T)@W
        err=norm(v.conj().T@h@v-a)
        assert err<1e-13
        energy.append(dict(model=name,initial_energy_operator_diagonal=np.diag(a).real.tolist(),
                           conservation_error=err,spectrum=np.linalg.eigvalsh(h).tolist()))
    paths=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/unified_operation_hypotheses_v0_1.md',
           HERE/'drafts/joint_model_identity_review.md',HERE.parent/'research_note_854.md',
           HERE.parent/'research_note_898.md',HERE.parent/'research_note_929.md',
           ROOT/'research_cognition_physics/archive_342_369/research_note_346.md']
    return dict(round=935,date='2026-10-07',all_scientific_checks_passed=True,
        one_bounded_common_source_identity_test=True,endpoint_isometry_errors=endpoint,
        full_input_energy_moments=moments,source_rows=rows,two_branch_matrix_rows=grams,
        derivative_checks=deriv,analytic_derivatives=dict(program=-float(np.pi/4),actual_rate=0.,direct_rate=-float(np.pi/4)),
        source_generator_difference=norm(dp-h0),energy_rows=energy,
        datum_energy_K_is_nontrivial=True,clock_baseline_accounted=True,
        exact_uniform_universe_time_rescaling_not_used_as_observable=True,
        pure_gate_dressed_independent_clock_class_only=True,
        original_Qeff_rejected=False,all_autonomous_implementations_rejected=False,
        gravity_or_standard_model_restored=False,full_goal_completed=False,
        clock_optimization_stopped=True,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--write',action='store_true');args=a.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
