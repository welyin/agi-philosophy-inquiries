"""Finite resolution qualification: exact budgets and finite CAR algebra, not field simulation."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
import numpy as np
BASE=Path(__file__).resolve().parent
SOURCES={
 'after1142_joint_contract_review.md':'dcb3316bcedcb86a4e0a7416652445b43d63e382d5e7c488a24929d950891b59',
 '../1133/proof.md':'476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
 '../1137/proof.md':'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb',
}

def run():
    hashes={n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in SOURCES}
    assert hashes==SOURCES
    r=F(17,4); lam=F(11,7); kap=F(11,14)
    coeff=[2*r*84,2*lam*16*r,2*kap*16*r*2]
    assert coeff==[F(714),F(1496,7),F(1496,7)] and sum(coeff)==F(7990,7)
    response=coeff[1]*F(1,1000)/25
    support_response=response*(r+F(1,10**5))/r
    assert response==F(187,21875) and support_response==F(4675011,546875000)
    assert support_response<F(854860,10**8)
    weak=F(1,10**4); relative_task=16*lam*weak
    original_task=2*F(1446465,10**9)
    assert relative_task+original_task<F(5408,10**6)<F(1,100)
    # One finite CAR representation only. Not the original continuum normal sea.
    n=4; dim=2**n; annih=[]
    for j in range(n):
        c=np.zeros((dim,dim),complex)
        for ket in range(dim):
            if (ket>>j)&1:
                sign=(-1)**((ket & ((1<<j)-1)).bit_count())
                c[ket^(1<<j),ket]=sign
        annih.append(c)
    eye=np.eye(dim)
    for i in range(n):
        for j in range(n):
            assert np.max(np.abs(annih[i]@annih[j].conj().T+annih[j].conj().T@annih[i]-(eye if i==j else 0)))<1e-14
    def q(a,S=None):
        out=sum((a[i,j]*annih[i].conj().T@annih[j] for i in range(n) for j in range(n)),np.zeros((dim,dim),complex))
        if S is not None: out-=np.trace(S@a)*eye
        return out
    def sn(a): return float(np.max(np.abs(np.linalg.eigvalsh(a))))
    def expu(a,t):
        e,u=np.linalg.eigh(a)
        return (u*np.exp(-1j*t*e))@u.conj().T
    rng=np.random.default_rng(114203)
    positions=np.array([-3.25,-1.25,1.25,3.25]); Rg=3.25; Rchi=1.
    x=np.array([-.75,0.,.75]); weights=np.array([.25,.5,.25]); L=.04
    wm=np.diag(.5+L*positions); ws=.5+L*x
    # Commuting source spectral fibers; phi functional calculus stays explicit.
    f=np.array([[.2,.8,.5],[1.,.3,.7],[.6,.4,.9]])
    C=np.diag(weights@f); Cw=np.diag((weights*ws)@f)
    srcdim=3; max_comp=0.; max_ratio=0.; max_duhamel=0.
    for trial in range(8):
        z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); A=(z+z.conj().T)/2
        S=np.diag([1.,0.,1.,0.]) if trial%2 else None
        As=(wm@A+A@wm)/2
        E_s=np.kron(Cw,q(A,S)); E_m=np.kron(C,q(As,S))
        V=np.kron(C,q(A,S))
        E_s_comp=np.kron(C-Cw,q(A,S))
        E_m_comp=np.kron(C,q(A-As,S))
        err=max(sn(E_s+E_s_comp-V),sn(E_m+E_m_comp-V));max_comp=max(max_comp,err)
        trace=float(np.sum(np.abs(np.linalg.eigvalsh(A))))
        bound=2*L*(Rchi+Rg)*trace
        gap=sn(E_s-E_m); assert gap<=bound+1e-12
        max_ratio=max(max_ratio,gap/bound)
        H=V+np.kron(np.diag([.1,.4,.9]),eye)+np.kron(np.eye(srcdim),q(np.diag([-.4,-.1,.2,.5]),S))
        strength=.1; time=.5
        U_s=expu(H+strength*E_s,time);U_m=expu(H+strength*E_m,time)
        actual=float(np.linalg.norm(U_s-U_m,2)); bd=strength*time*gap
        assert actual<=bd+2e-13
        max_duhamel=max(max_duhamel,actual/bd if bd else 0.)
    assert max_comp<1e-12
    return dict(status='PASS',sources=hashes,groups={
      'exact_contract_budgets':dict(coefficients=[str(v) for v in coeff],large_response=str(response),source_support_response=str(support_response),weak_response=str(support_response/10),task_with_response=str(relative_task+original_task)),
      'finite_CAR_algebra':dict(trials=8,car_modes=n,dimension=dim*srcdim,complement_residual=float(max_comp),allocation_bound_ratio=float(max_ratio),duhamel_bound_ratio=float(max_duhamel))},
      normal_continuum_field_simulated=False,new_response_materialized=False,same_four_momentum_response_proved=False,
      full_regional_balance_proved=False,formal_round_reserved=False,new_axioms=0,science_increment=0)

def main():
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    got=run();out=BASE/'after1142_joint_contract_results.json'
    if a.write_results:
        assert not out.exists()
        out.write_text(json.dumps(got,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: assert json.loads(out.read_text(encoding='utf-8'))==got
    print(json.dumps(dict(status=got['status'],groups=len(got['groups']),saved_results_match=not a.write_results,formal_round_reserved=False)))
if __name__=='__main__':main()
