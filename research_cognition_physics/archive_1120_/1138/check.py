"""Round 1138: explicit conditional finite-resolution energy certificate.
Default execution is read-only. Numerical diagnostics do not prove field domains.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, subprocess, sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FROZEN = {'1137/proof.md': 'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb', '1137/results.json': '1fb7f63266bc240d80ab880af707414d4e41a92c2075e14b69459e53b819ceb4', '1137/acceptance.json': '8101ad968bca9dcc40ccb19c35280a69d64c23cf84603e11788cea7b1b69adad', '_shared/after1137_direction_decision.md': '42c0ec818799ac24bf81310ef8a628929a2d4e3b0c558625c2e3af534eb2ee83', '_shared/after1137_local_balance_preflight.md': '95915f871317e8638d23c3ee29cec779a0ff52aad5f76873a0ff896f6a2fcbe0', '_shared/local1138_dynamics_preflight.md': '56a559135a1288336d38907cfdf14e7fa5d23b4f3d1ec7268f83b97936d7f1b9', '_shared/local1138_dynamics_preflight.py': '53ca56547c0a017d35df822b7a6746a52ff4f5e5af612cb3fef884bbd1959a77', '_shared/local1138_dynamics_preflight.json': '3482a6c485289c414e4bae28243ac052902ae067c4eb987609c14aaf3f1a13db', '_shared/local1138_kernel_completion.md': '963b64fe27a133ec2c3eac580252e669ec95d4454d4bc85de91f9601323a24fd', '_shared/local1138_kernel_completion.py': 'b5b70ea0ef7e167cbc62806335638503eb7662f23b3c41a90755a484286343ae', '_shared/local1138_kernel_completion.json': 'ab363ecc68499a6f9539c376d973efaf8e378914e9c49281c3081acfcc328ca1', '_shared/local1138_scope_preflight.md': 'aac764bba2a971a175a7dcd38e52cc7c06d4a4c884f77186d7f3b22e9cad0bff', '_shared/local1138_joint_preflight.md': 'b281cc1c63ef64c95f0f22bfb591048ab1ba7ea5f71614e4c494a0355a6064dc', '_shared/after1137_local_balance_admission.md': 'abcd142ce2aec9a6c30b808e176af8f017c4855c753e688fc7367ca3c27647db'}

def pack(d):
    return {'fractions':{k:str(v) for k,v in d.items()},
            'decimal':{k:float(v) for k,v in d.items()}}

def exact():
    A=F(10**34); M=F(10**80); sig=F(1,10**30); eps=F(1,10**5)
    lam=F(11,7); kap=F(11,14); cs=F(33536000000000000000)
    cy=F(24780800000000000000000,49)+16316
    v=16*lam*cs+F(558,7)
    dv=lam*cs*512+42*F(8,1)*10**6+(lam+kap)*512
    ds=F(12)/(eps*A)*dv
    dy=cy*(8*sig+8/(M*sig))
    dr=F(3242,10**9)
    assert ds+dy<dr and dv<28*10**21
    material=F(268)+56+2*v+4*dv+F(1)
    assert material<2*10**23
    B=F(1074*10**68); de=F(10**60); ce=F(650); hl=F(3201,10**11)
    source_graph=1073*A*A+94*A
    assert source_graph+material+1<B
    NE=ce*B+de; NJ=4*NE
    y_two=12876*sig*A*A
    material_two=2*(261*material+de)
    # de here overestimates all finite internal/pair/V additive terms.
    other=F(10**66)
    assert y_two+material_two<other
    baseline=F(27,8)*F(314,100)*A*A
    I=baseline-4*dr*NE-other
    J_lower=2*I-2*hl*NJ
    balance=4*hl*NE+hl*NJ
    S=F(25*10**67); N=F(28*10**73)
    assert J_lower>S and NJ<N and balance<18*10**66
    tau=S/(2*N*N); bias=3*S*S/(16*N*N); disturbance=S/(4*N)
    calibration=F(1,10**14)
    old=F(1446465,10**9); full=old+disturbance
    assert bias>calibration and full<F(1447,10**6)
    return pack(dict(delta_source=ds,delta_Y=dy,delta_reference_safe=dr,
      derivative_material_bound=dv,reference_material_graph_bound=material,
      reference_source_graph_bound=source_graph,common_B_graph=B,
      normal_Y_two_endpoints_bound=y_two,other_subtraction=other,
      cE=ce,dE=de,N_E=NE,N_J=NJ,source_baseline_lower=baseline,
      reference_integrated_flux_lower=I,actual_J_lower=J_lower,
      three_point_balance_bound=balance,S=S,N=N,tau=tau,
      ideal_readout_bias_lower=bias,probability_calibration_budget=calibration,
      bias_after_calibration=bias-calibration,full_readout_disturbance_bound=disturbance,
      complete_old_task_with_readout_bound=full))

def evolve(H,t,psi):
    w,u=np.linalg.eigh(H)
    return u @ (np.exp(-1j*t*w)[:,None]*(u.conj().T@psi))

def average(H,O,t,h):
    w,u=np.linalg.eigh(H); x=w[:,None]-w[None,:]
    return u @ ((u.conj().T@O@u)*np.exp(1j*t*x)*np.sinc(x*h/(2*np.pi))) @ u.conj().T

def three_point():
    rng=np.random.default_rng(1138)
    n=4; r=3; a=25/8; b=29/8; t=27/8; h=.5
    L=np.diag([.2,1.,2.4,4.1])
    x=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    R=x.conj().T@x; R=R/np.linalg.norm(R,2)*.07
    H=L+R
    z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    E=(z+z.conj().T)/5
    J0=-1j*(L@E-E@L)
    J=average(L,J0,0,h)
    psi=rng.normal(size=(n,r))+1j*rng.normal(size=(n,r))
    psi=psi/np.linalg.norm(psi)
    expect=lambda O,v: float(np.real(np.vdot(v,O@v)))
    lb=expect(E,evolve(L,a,psi))-expect(E,evolve(L,b,psi))-h*expect(J,evolve(L,t,psi))
    hb=expect(E,evolve(H,a,psi))-expect(E,evolve(H,b,psi))-h*expect(J,evolve(H,t,psi))
    ds={s:float(np.linalg.norm(evolve(H,s,psi)-evolve(L,s,psi))) for s in [a,b,t]}
    bound=2*np.linalg.norm(E,2)*(ds[a]+ds[b])+2*h*np.linalg.norm(J,2)*ds[t]
    actual_history=average(H,J0,t,h)
    # Compare distinct operators, not only an accidentally cancelling state.
    w,u=np.linalg.eigh(H)
    Ut=u@np.diag(np.exp(1j*t*w))@u.conj().T
    history_mismatch=np.linalg.norm(actual_history-Ut@J@Ut.conj().T,2)
    assert abs(lb)<1e-12 and abs(hb)<=bound+1e-12 and history_mismatch>1e-5
    return dict(reference_balance_residual=float(abs(lb)),
      actual_three_point_residual=float(abs(hb)),comparison_bound=float(bound),
      H_history_vs_L_calibrated_center_operator_difference=float(history_mismatch),
      entangled_reference_dimension=r,
      field_domain_proof_supplied_by_this_diagnostic=False)

def reader():
    rng=np.random.default_rng(11380)
    x=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    J=x.conj().T@x/5+np.eye(4)*.1
    w,u=np.linalg.eigh(J)
    psi=rng.normal(size=(4,3))+1j*rng.normal(size=(4,3)); psi/=np.linalg.norm(psi)
    S=float(np.real(np.vdot(psi,J@psi))); N=float(np.linalg.norm(J@psi))
    tau=S/(2*N*N)
    c=u@np.diag(np.cos(tau*w/2))@u.conj().T
    s=u@np.diag(np.sin(tau*w/2))@u.conj().T
    ks=[(c+s)/np.sqrt(2),(c-s)/np.sqrt(2)]
    complete=sum(k.conj().T@k for k in ks)
    p=float(np.linalg.norm(ks[0]@psi)**2)
    out=np.stack([k@psi for k in ks]); ideal=np.stack([psi,psi])/np.sqrt(2)
    disturbance=float(np.linalg.norm(out-ideal))
    lower=3*S*S/(16*N*N)
    residual=float(np.linalg.norm(complete-np.eye(4),2))
    assert residual<1e-12 and p-.5>=lower-1e-12 and disturbance<=S/(4*N)+1e-12
    return dict(completeness_residual=residual,actual_bias=p-.5,
      certified_bias=lower,complete_reference_preserving_disturbance=disturbance,
      certified_disturbance=S/(4*N),arbitrary_reference_dimension=3,
      autonomous_detector_claimed=False)

def frozen():
    out={}
    for rel,expected in FROZEN.items():
        actual=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
        assert actual==expected, rel
        out[rel]=actual
    return out

def constituent_replays():
    rows={}
    for name in ['local1138_dynamics_preflight.py','local1138_kernel_completion.py']:
        p=ROOT/'_shared'/name
        run=subprocess.run([sys.executable,'-B','-X','utf8',str(p)],capture_output=True,text=True,encoding='utf-8')
        assert run.returncode==0,(name,run.stdout,run.stderr)
        rows[name]={'exit_code':0,'saved_results_verified':True}
    return rows

def compute():
    return dict(schema='round1138_finite_three_point_certificate_v1',
      round=1138,new_axioms=0,scientific_calibration_increment=0,
      groups=dict(frozen_sources=frozen(),constituent_replays=constituent_replays(),
        exact_finite_budget=exact(),noncommuting_three_point_diagnostic=three_point(),
        complete_weak_reader_diagnostic=reader()),
      scope=dict(actual_H_on_same_bare_L_calibrated_observables=True,
        actual_H_continuous_history_flux_claimed=False,
        endpoint_E_readouts_quantitatively_certified=False,
        J_reader_is_explicit_new_effective_permission=True,
        full_poststate_conclusions_at_given_Kraus_level=True,
        autonomous_detector_work_closed=False,
        local_stress_and_full_field_stage_complete=False))

def check_close(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: check_close(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):check_close(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=3e-12,abs_tol=1e-13),(a,b)
    else: assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    data=compute(); result=HERE/'results.json'
    if args.write_results:
        assert not result.exists(),'Refuse to overwrite frozen result'
        result.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: check_close(data,json.loads(result.read_text(encoding='utf-8')))
    print('Round1138 PASS: 5 groups; frozen sources, constituent replays, exact budget, noncommuting three-point and full-reference weak reader. Conditional effective scope only.')
