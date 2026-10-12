"""1147: exact budgets plus independent finite algebra diagnostics.
Default readonly. Continuous Fock/Dirac statements are analytic, not simulated here.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent

def pack(q): return {'fraction':str(q),'float':float(q)}
def budget():
    # Use only the three analytic inequalities in the fixed lambda choice.
    q=Q(1,40000); p=q/(1-q); full=4*p
    relative_eta=Q(8,300000); clock=Q(4,300000)
    pi_hi=Q(22,7); dt=Q(1,10000)
    phase=pi_hi*(dt+(1+dt)*relative_eta)
    rec=full+phase; instrument=rec+Q(2,1000)
    early=full+pi_hi/2*clock
    # A conservative all-parameter T0 lower bound; C2<3, R<=1/60.
    lambda_hi=Q(1,600000); time_lo=Q(3)/(3*lambda_hi**2)
    assert full<Q(100003,10**9)
    assert rec<Q(499,10**6) and instrument<Q(25,10000)
    assert 1-2*instrument>Q(995,1000)
    assert early<Q(121,10**6) and Q(1,2)-early>Q(499879,10**6)
    assert clock<1 and time_lo>1
    return {'pass':True,'kind':'exact_rational_uniform_budget_no_spectral_truncation',
      'lambda_times_D_upper':pack(q),'Riesz_projection_upper':pack(p),
      'all_time_complete_isometry_upper':pack(full),'relative_conditional_energy_upper':pack(relative_eta),
      'relative_actual_clock_shift_upper':pack(clock),'timing_relative_tolerance':pack(dt),
      'record_isometry_upper':pack(rec),'with_phase_and_instrument_upper':pack(instrument),
      'same_effect_record_probability_gap_lower':pack(1-2*instrument),
      'actual_clock_isometry_upper':pack(early),'clock_readout_change_lower':pack(Q(1,2)-early),
      'gate_time_lower_not_feasibility_estimate':pack(time_lo),
      'exact_spectral_constants_not_numerically_certified':['Delta_D','C_epsilon','gamma'],
      'no_truncated_spectrum_used_in_error_certificate':True}

def second_order():
    # Only one virtual boson is required for this second-order coefficient.
    # Retain one excited material state in each invariant branch. This tests
    # the cross-term identity, not the infinite material/field approximation.
    omega=2.; ga=.75; gb=2/3; j=[0.,.4]; eye=np.eye(2)
    X=np.array([[0.,1.],[1.,0.]])/math.sqrt(2*omega)
    gamma=ga*gb/omega**2; want=-gamma*(j[1]-j[0])**2
    rows=[]
    for pair_strength in (.2,.7):
        coeff={}; singles={}
        for a in range(2):
            for b in range(2):
                da=3.+a; db=4.+2*b
                HA=np.diag([0.,da]); HB=np.diag([0.,db]); HF=np.diag([0.,omega])
                JA=np.array([[j[a],pair_strength],[pair_strength,.17]])
                JB=np.array([[j[b],pair_strength+.1],[pair_strength+.1,-.11]])
                H0=np.kron(np.kron(HA,eye),eye)+np.kron(np.kron(eye,HB),eye)+np.kron(np.kron(eye,eye),HF)
                V=ga*np.kron(np.kron(JA,eye),X)+gb*np.kron(np.kron(eye,JB),X)
                resolvent=np.diag([0.]+[1/x for x in np.diag(H0)[1:]])
                e2=-float((V@resolvent@V)[0,0]); coeff[str(a)+str(b)]=e2
                LA=-ga**2/(2*omega)*(j[a]**2/omega+pair_strength**2/(da+omega))
                LB=-gb**2/(2*omega)*(j[b]**2/omega+(pair_strength+.1)**2/(db+omega))
                predicted=LA+LB-gamma*j[a]*j[b]
                assert abs(e2-predicted)<1e-14
                singles[str(a)+str(b)]=[LA,LB]
        eta=coeff['11']-coeff['10']-coeff['01']+coeff['00']
        assert abs(eta-want)<1e-14 and coeff['00']<0
        rows.append({'pair_strength':pair_strength,'branch_second_order':coeff,
                     'local_self_terms':singles,'four_branch_difference':eta,
                     'cross_identity_residual':abs(eta-want)})
    assert abs(rows[0]['branch_second_order']['00']-rows[1]['branch_second_order']['00'])>.001
    # Independent same bare readout and complete unknown-input record algebra.
    T=1.37; EA=.31; EB=.47; E00=-.23
    phase=np.exp(-1j*T*np.array([E00,E00+EB,E00+EA,E00+EA+EB]))
    U=np.diag(phase*np.array([1,1,1,-1]))
    plus=np.array([1.,1.])/math.sqrt(2)
    W=np.kron(np.eye(2),plus[:,None])
    v=np.array([1,np.exp(-1j*EB*T)])/math.sqrt(2)
    F=np.kron(np.eye(2),np.outer(v,v.conj()))
    decoded=(U@W).conj().T@F@(U@W)
    assert np.linalg.norm(decoded-np.diag([1.,0.]))<1e-13
    # A small phase mis-setting changes the full two-outcome Lüders isometry.
    eps=.001; w=np.array([1,np.exp(-1j*(EB*T+eps))])/math.sqrt(2)
    P=np.outer(v,v.conj()); Qp=np.outer(w,w.conj())
    delta=np.vstack([P-Qp,-P+Qp]); instrument_norm=float(np.linalg.norm(delta,2))
    assert instrument_norm<eps
    return {'pass':True,'kind':'second_order_four_branch_and_same_bare_record_algebra',
      'branch_diagnostic_dimension':8,'gamma':gamma,'predicted_mixed_coefficient':want,
      'cases':rows,'E00_response_preserved':True,
      'bare_record_effect_residual':float(np.linalg.norm(decoded-np.diag([1.,0.]))),
      'phase_misset_Lueders_isometry_norm':instrument_norm,
      'not_a_continuous_field_simulation':True,'no_long_time_truncation_claim':True}

def load_material():
    spec=importlib.util.spec_from_file_location('round1147_material',HERE/'material_check.py')
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    actual=mod.run(); saved=json.loads((HERE/'material_results.json').read_text(encoding='utf-8'))
    mod.compare(saved,actual)
    return actual

def run():
    sources=json.loads((HERE/'source_hashes.json').read_text(encoding='utf-8'))
    for name,sha in sources.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,name
    material=load_material()
    return {'round':1147,'all_passed':True,'check_groups':4,'scientific_check_increment':0,
            'frozen_source_count':len(sources),'groups':[budget(),second_order()],
            'material_results':material}

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+'/'+str(i))
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
    else:assert a==b,(path,a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    actual=run(); target=HERE/'results.json'
    if args.write:
        assert not target.exists(),'saved results already exist'
        target.write_text(json.dumps(actual,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: compare(json.loads(target.read_text(encoding='utf-8')),actual)
    print(json.dumps({'round':1147,'status':'PASS','check_groups':4,'saved_results_match':not args.write,'frozen_sources':actual['frozen_source_count']}))
