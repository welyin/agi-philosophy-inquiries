"""Default read-only replay; --write-results explicitly writes this module's JSON."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
EXPECTED={
 'local1138_dynamics_preflight.md':'56a559135a1288336d38907cfdf14e7fa5d23b4f3d1ec7268f83b97936d7f1b9',
 'local1138_scope_preflight.md':'aac764bba2a971a175a7dcd38e52cc7c06d4a4c884f77186d7f3b22e9cad0bff',
 '../1137/proof.md':'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb',
}
def arctan_enclosure(x:F,n:int):
    p=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    q=p+(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(p,q),max(p,q)
def rational_bounds():
    a,b=arctan_enclosure(F(1,5),16);c,d=arctan_enclosure(F(1,239),5)
    pl,pu=16*a-4*d,16*b-4*c
    assert F(3141592653589,10**12)<pl<pu<F(3141592653590,10**12)
    width=F(49,50)
    assert F(35,16)/width<3
    assert F(84,5)**2/5/width**4<8**2
    assert F(105,2)/width**3<84
    M0,M1,M2=F(1437),F(3117),F(4797)
    assert 4*pu*7**3/3<M0
    assert 9*4*pu*(7**3-6**3)/3<M2
    assert M0*M2<M1**2
    vol=4*pu*127/3
    D2=81*vol;D3=F(521,6)**2*vol
    f0sq=(M0+2*M2+D2)/(8*pl)
    f1sq=3*(M0+3*M2+3*D2+D3)/(32*pl)
    assert f0sq<50**2 and f1sq<420**2
    ip={n:4*16**(3+n)*(F(1,3+n)+F(1,3-n)) for n in (0,1,2)}
    den=8*pl**2
    k2=ip[0]*M0/den
    eta2=(ip[1]*M0+ip[0]*M1/2)/den
    beta2=(ip[2]*M0+ip[1]*M1+ip[0]*M2/4)/den
    bf2=4*ip[0]*M2/den
    assert k2<447**2 and eta2<4*10**6 and beta2<10**8 and bf2<3*10**6
    assert F(10**8,2)<7072**2 and 20000+7072<30000
    assert 547**2+261**2+1<607**2<650**2
    v=16*F(11,7)*33536000000000000000+F(558,7)
    assert v+28+30000+3000001<10**60
    values={'pi_lower':pl,'pi_upper':pu,'F0_squared_upper':f0sq,'F1_squared_upper':f1sq,
      'kappa_squared_upper':k2,'eta_squared_upper':eta2,'beta_squared_upper':beta2,
      'dirac_beta_squared_upper':bf2,'bounded_V_upper':v}
    return {'exact':{k:str(v) for k,v in values.items()},'display':{k:float(v) for k,v in values.items()},
      'certificate':{'scalar_c':547,'scalar_d':30000,'dirac_c':261,'dirac_d':3000001,
      'energy_c':650,'energy_d':str(10**60),'averaged_flux_c':2600,'combined_c_squared':367331}}
def number_graph_diagnostic():
    h=np.diag([0.4,1.7,3.1]);a=np.array([[0.2,0.4j,0.1],[-0.4j,-0.3,0.2],[0.1,0.2,0.8]],complex)
    hi=np.diag(1/np.sqrt(np.diag(h)))
    c0=np.max(np.abs(np.linalg.eigvalsh(hi@a@hi)))
    c1=np.linalg.norm(a@np.linalg.inv(h),2);cap=max(c0,c1)
    at=np.kron(a,np.eye(3))+np.kron(np.eye(3),a)
    ht=np.kron(h,np.eye(3))+np.kron(np.eye(3),h)
    slack=np.linalg.eigvalsh(cap**2*(ht@ht)-at@at).min()
    assert slack>-1e-10 and np.linalg.norm(a@h-h@a)>0.1
    return {'c0':float(c0),'c1':float(c1),'quadratic_slack_min':float(slack)}
def ccr_diagnostic():
    n=7;ann=np.diag(np.sqrt(np.arange(1,n)),1);eye=np.eye(n)
    aa=[np.kron(ann,eye),np.kron(eye,ann)]
    om=np.array([0.7,1.3]);B=np.array([[0.4,0.2j],[0.2j,-0.3]],complex)
    hf=sum(om[j]*(aa[j].conj().T@aa[j]) for j in range(2))
    Q=sum(B[i,j]*(aa[i]@aa[j])/2 for i in range(2) for j in range(2))
    bb=B.conj().T@B
    expected=sum(bb[i,j]*(aa[i].conj().T@aa[j]) for i in range(2) for j in range(2))+np.linalg.norm(B)**2*np.eye(n*n)/2
    inds=[i*n+j for i in range(n) for j in range(n) if i+j<=3]
    select=np.eye(n*n)[:,inds]
    defect=np.linalg.norm((Q@Q.conj().T-Q.conj().T@Q-expected)@select)
    assert defect<1e-12
    kap=np.linalg.norm(B/np.sqrt(om[:,None]*om[None,:]))
    eta=np.linalg.norm(B/np.sqrt(om[:,None]));beta=np.linalg.norm(B)
    rng=np.random.default_rng(1138);ratio=0.0
    for _ in range(40):
        z=rng.normal(size=len(inds))+1j*rng.normal(size=len(inds));z/=np.linalg.norm(z);psi=select@z
        num=np.linalg.norm((Q+Q.conj().T)@psi)
        den=kap*np.linalg.norm(hf@psi)+eta*math.sqrt(np.vdot(psi,hf@psi).real)+beta/math.sqrt(2)
        ratio=max(ratio,num/den)
    assert ratio<1
    return {'input_max_total_particles':3,'output_supported_total_particles':5,'commutator_defect':float(defect),'maximum_bound_ratio':float(ratio)}
def car_diagnostic():
    sm=np.array([[0,1],[0,0]],complex);zz=np.diag([1,-1]);ee=np.eye(2);ops=[]
    for j in range(4):
        x=np.array([[1]],complex)
        for k in range(4):x=np.kron(x,zz if k<j else sm if k==j else ee)
        ops.append(x)
    s=[0.7,0.2];Q=sum(s[j]*(ops[2*j+1]@ops[2*j]) for j in range(2))
    nn=[a.conj().T@a for a in ops];ident=np.eye(16)
    exact=sum(s[j]**2*(ident-nn[2*j]-nn[2*j+1]) for j in range(2))
    defect=np.linalg.norm(Q@Q.conj().T-Q.conj().T@Q-exact)
    assert defect<1e-12
    return {'all_CAR_dimension':16,'pair_commutator_defect':float(defect)}
def generate():
    observed={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in EXPECTED}
    assert observed==EXPECTED,{'frozen_source_mismatch':observed}
    return {'status':'PASS_QUALIFICATION_ONLY_NOT_ROUND_ACCEPTANCE','frozen_sources':observed,
      'proof_sha256':hashlib.sha256((HERE/'local1138_kernel_completion.md').read_bytes()).hexdigest(),
      'rational_bounds':rational_bounds(),'diagnostics':{'noncommuting_number_graph':number_graph_diagnostic(),
      'CCR_full_output_on_input_subspace':ccr_diagnostic(),'CAR_complete_finite_algebra':car_diagnostic()},
      'scope':'Analytic all-Fock/domain proof is in markdown; numerical diagnostics are finite checks only.'}
def main():
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    fresh=generate();dest=Path(__file__).with_suffix('.json')
    if args.write_results:dest.write_text(json.dumps(fresh,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(dest.read_text(encoding='utf-8'))
        assert saved==fresh,'saved certificate differs from read-only recomputation'
    print(json.dumps({'status':fresh['status'],'mode':'write-results' if args.write_results else 'read-only-replay',
      'frozen_source_count':len(EXPECTED),'groups':['rational_kernel_bounds','noncommuting_dGamma','CCR','CAR'],
      'cE':650,'dE':'1e60'},ensure_ascii=False))
if __name__=='__main__':main()
