"""Nonformal original-input multileg-tail qualification. Default replay is readonly.
Exact budgets check the written analysis; they do not prove continuous operator domains.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
SOURCES={
    '1129/proof.md':'c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa',
    '1133/proof.md':'476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
    '1137/proof.md':'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb',
    '1141/proof.md':'d71543b186660711ac1602a15495c5700ba2430577366d3057f70249fb650c96',
    '1142/proof.md':'9b25627319b3c0fb2a976aa8c27e8ea9c723bbe174539c7a4b29d79a7753012f',
}
def encode(x):
    if isinstance(x,F):return {'fraction':str(x),'float':float(x)}
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,list):return [encode(v) for v in x]
    return x

def calculate():
    found={k:hashlib.sha256((BASE/k).read_bytes()).hexdigest() for k in SOURCES}
    assert found==SOURCES
    eps=F(1,10**5); A=F(10**34); M=F(10**80); T=F(5,4)
    e=F(15000); y=F(3,10**30); mu=5/eps
    # 2 sqrt(mu e) <174000, no floating-point square root in proof budget.
    ann=F(174000); assert ann**2>4*mu*e
    fluct=ann+3/eps; mean_tail=T*A*y
    rms=F(600000); assert fluct+mean_tail<rms
    r=rms/A; LS=F(8*10**6); LF=F(512)
    lam=F(11,7); kap=F(11,14); L2S=F(128*10**12); L2F=F(32768)
    vrec_state=r*(42*LS+(lam+kap)*LF)
    force_raw=3*r*(42*L2S+(lam+kap)*L2F)
    force_column=force_raw/A
    power_res=T*24*A*force_column
    storage_source=T*4*A*force_column
    resolution_residual=power_res+storage_source
    assert vrec_state<F(21,10**21)
    assert force_column<F(1,10**46)
    assert power_res<F(3,10**11)
    assert storage_source<F(5,10**12)
    assert resolution_residual<F(4,10**11)
    # Physical free-source pi mean is A on the near-Y support; global |pi_cl|<=11A.
    # Volume of B_(1+2eps) < 5. sqrt(5)<3, relative column <4e12.
    assert (1+2*eps)**3*F(88,21)<5
    assert 3*A*(1+12*y)+4*10**12<4*A
    dirac=r*(7204*LS+(lam*42+kap*84)*LF)
    # V1129's one-particle generator in each commuting source fiber has norm<=16 lam.
    v1=16*lam
    assert 8*v1*21==4224
    old_handoff=r*(4224*LS+(lam*(2*v1)+kap*(4*v1))*LF)
    receiver_self=2*F(558,7)*vrec_state
    angle_column=2*r*(29032*LS+(lam*132+kap*264)*LF)
    Bstar=F(1074,1000)*10**71; inertia=M; Lbound=F(5*10**75)
    assert Lbound**2>2*inertia*Bstar
    rotor=Lbound*angle_column/inertia
    P2=F(1201,1000)*10**72
    SY_residual=2*P2/M
    support=2*SY_residual*vrec_state
    Rbound=F(3,10**8)
    assert F(7,4)*P2/M<Rbound
    bt=2*Rbound*vrec_state
    assert T*dirac<F(5,10**18)
    assert T*old_handoff<F(3,10**18)
    assert T*receiver_self<F(5,10**18)
    assert T*rotor<F(2,10**21)
    assert T*support<F(2,10**27)
    assert T*bt<F(2,10**27)
    Eold=F(1447,10**6)
    reference_S_deficit=LS*3/(eps*A)+y
    actual_S4=1-2*Eold-reference_S_deficit
    assert actual_S4>F(997,1000)
    # Algebraic source-resolution identity is checked symbolically as a word decomposition.
    # In the analysis: Q= <w pi,C'>/A and P=-<K^2(w pi),C'>/A.
    # Its residual is ((I-K^2)w)pi = w(I-K^2)pi - [K^2,w]pi.
    # This finite exact matrix diagnostic only checks this algebra, not field physics.
    import numpy as np
    K=np.array([[1.,.25,0.],[.25,.75,.125],[0.,.125,.5]])
    W=np.diag([1.,.5,0.]); pi=np.array([1.,-2.,3.]); c=np.array([2.,1.,-1.])
    K2=K@K
    lhs=np.dot(W@pi,c)-np.dot(K2@W@pi,c)
    rhs=np.dot(W@(np.eye(3)-K2)@pi-(K2@W-W@K2)@pi,c)
    assert abs(lhs-rhs)<1e-12
    return encode({
      'status':'PASS','formal_round_reserved':False,'new_axioms':0,'science_increment':0,
      'scope':'original normal inputs, w, b=0, 0<=t<=5/4; receiver tails only',
      'frozen_sources':found,
      'actual_source_rms_numerator_upper':rms,'dimensionless_source_rms_upper':r,
      'actual_receiver_storage_vector_upper':vrec_state,
      'actual_receiver_storage_endpoint_change_upper':2*vrec_state,
      'receiver_source_force_column_upper':force_column,
      'resolved_receiver_power_L1_upper':power_res,
      'full_source_to_receiver_storage_L1_upper':storage_source,
      'source_resolution_cut_residual_L1_upper':resolution_residual,
      'other_storage_exchange_L1_upper':{
         'Dirac_free':T*dirac,'old_handoff':T*old_handoff,
         'receiver_self':T*receiver_self,'rotor':T*rotor,
         'support_kinetic':T*support,'BT_H_minus_Ho':T*bt,'Z':F(0)},
      'actual_holding_coefficient_at_T4_lower':actual_S4,
      'source_resolution_identity_diagnostic_residual':float(abs(lhs-rhs)),
      'continuous_domains_proved_by_code':False,
      'active_receiver_localized':False,'complete_local_stress_proved':False,
      'normal_sea_unknown_Q_and_reference_retained':True,
      'original_parameters_changed':False,
    })
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=calculate()
    out=HERE/'after1142_multileg_results.json'
    if args.write_results:
        assert not out.exists(),'new qualification must not overwrite an existing result'
        out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(out.read_text(encoding='utf-8'))==result,'saved result mismatch'
    print(json.dumps({'status':'PASS','frozen_sources':len(SOURCES),'exact_budget_groups':3,'saved_results_match':not args.write_results,'full_local_stress_proved':False},ensure_ascii=False))
