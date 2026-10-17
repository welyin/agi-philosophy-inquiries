"""Hard energy-support diagnostic for the same fixed reference trace."""
from pathlib import Path
import importlib.util
import json
import math
import sys
import numpy as np

BASE=Path(__file__).with_name('finite_observer_entropy_check_after1270.py')
spec=importlib.util.spec_from_file_location('base_entropy',BASE)
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
RESULT=Path(__file__).with_suffix('.json')


def cdf(k,b):
    return 1-math.exp(-b)*sum(b**n/math.factorial(n) for n in range(k))


def quad(cuts,order):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    return (np.concatenate([(a+b)/2+(b-a)*nodes/2 for a,b in zip(cuts[:-1],cuts[1:])]),
            np.concatenate([(b-a)*weights/2 for a,b in zip(cuts[:-1],cuts[1:])]))


def calculate(order):
    e=np.array([0.,.5,1.5]); beta=math.log(2); z=float(np.exp(-beta*e).sum())
    psi=np.array([1,1j,2])/math.sqrt(6)
    pure=np.outer(psi,psi.conj()); diag=np.diag([.2,.3,.5]).astype(complex)
    states={'commuting':diag,'pure_coherent':pure,'mixed_coherent':.65*pure+.35*diag}
    rows=[];k=4
    for b in [12.,20.]:
        prob=cdf(k,b)
        xq,wq=quad(list(np.arange(0,b+1,1)),order)
        pq=base.gamma_pdf(xq,k)/prob
        h_unit=float(-np.dot(wq,pq*np.log(pq)))
        for width in [4.,16.]:
            shifts=e/width
            cuts=sorted(set([0.,*np.arange(1.,b+1.,1.),*shifts,*(b+shifts)]))
            xs,ws=quad(cuts,order)
            gp=np.stack([np.where(xs-shift<=b,base.gamma_pdf(xs-shift,k)/prob,0)
                         for shift in shifts],axis=1)
            for name,rho in states.items():
                pi=rho.diagonal().real; root=base.sqrtm(rho);srho=base.entropy(rho)
                py=gp@pi; mean=width*float(np.dot(ws,xs*py))
                s_fiber=avg_r=iyj=kl=0.
                for x,w,g,pr in zip(xs,ws,gp,py):
                    if pr<1e-280:continue
                    amp=np.sqrt(g);d=amp[:,None]*rho*amp[None,:]
                    s_fiber+=w*base.entropy(d)
                    ru=root.T@np.diag(g)@root.conj()
                    avg_r+=w*pr*base.entropy(ru/pr)
                    mask=g>0
                    iyj+=w*float(np.dot(pi[mask]*g[mask],np.log(g[mask]/pr)))
                    g0=math.exp((k-1)*math.log(x)-x-math.lgamma(k))
                    kl+=w*float(np.dot(pi[mask]*g[mask],np.log(g[mask]/g0)))
                exact=s_fiber+math.log(width)+math.log(beta)-math.log(z)-beta*mean
                mean_q=width*k*cdf(k+1,b)/prob
                leading=srho-beta*float(pi@e)-math.log(z)-beta*mean_q+h_unit+math.log(width)+math.log(beta)
                gap=iyj-srho+avg_r
                bound=-math.log(prob)+float(pi@e**2)/(2*(k-2)*width**2*prob)
                assert abs(float(ws@py)-1)<2e-11
                assert abs(mean-mean_q-float(pi@e))<2e-10
                assert abs(exact-leading-gap)<2e-10
                assert -2e-10<=gap<=iyj+2e-10 and iyj<=kl+2e-10 and kl<=bound+2e-10
                if name=='commuting':assert abs(gap)<2e-10
                if name=='pure_coherent':assert abs(iyj-gap)<2e-10
                rows.append(dict(state=name,k=k,width=width,cutoff_multiple=b,
                    energy_hard_cap=b*width,observer_mean=mean_q,retained_probability=prob,
                    entropy_exact=exact,entropy_leading=leading,gap=gap,
                    I_Y_energy=iyj,weighted_shift_KL=kl,bound=bound,
                    identity_error=abs(exact-leading-gap)))
    return rows


def main():
    a,b=calculate(32),calculate(64)
    keys=['entropy_exact','entropy_leading','gap','I_Y_energy','weighted_shift_KL','bound']
    err=max(abs(x[k]-y[k]) for x,y in zip(a,b) for k in keys)
    assert err<2e-9
    result=dict(status='PASS_FINITE_HARD_OBSERVER_SUPPORT',scope='State energy support only; original corner and trace kept; no finite conjugate-momentum claim',
        cases=len(b),quadrature_comparison_error=err,max_identity_error=max(x['identity_error'] for x in b),rows=b)
    if sys.argv[1:]==['--write-once']:
        with RESULT.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        action='WRITE_ONCE'
    elif not sys.argv[1:]:
        old=json.loads(RESULT.read_text('utf-8'))
        assert old['status']==result['status'] and old['cases']==result['cases']
        for x,y in zip(old['rows'],b):
            for k in keys:assert abs(x[k]-y[k])<2e-9
        action='READ_ONLY_MATCH'
    else:raise SystemExit('Expected no arguments or --write-once')
    print(json.dumps({k:result[k] for k in ['status','cases','quadrature_comparison_error','max_identity_error']}|{'action':action}))


if __name__=='__main__':main()
