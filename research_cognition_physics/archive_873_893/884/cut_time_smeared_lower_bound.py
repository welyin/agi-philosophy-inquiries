"""884 working: positive diagonal lower bound survives any fixed real
C-infinity compact time test with integral1. It is a lower bound, not the full
smeared variance. Original64 matrices and paired opposite momenta are kept.
"""
from pathlib import Path
import argparse,json
import numpy as np
import cut_source_variance_probe as probe
old=probe.old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'cut_time_smeared_lower_bound_results.json'

def diagonal_lower(H,G):
    ev,v=np.linalg.eigh(H)
    # Stable thermal occupation; the product n(1-n) is even in energy.
    r=np.exp(-old.BETA*abs(ev))
    weight=r/(1+r)**2
    g=v.conj().T@G@v
    return float(np.sum(weight*abs(np.diag(g))**2))

def run():
    m=old.load();q,_=old.charges(m)
    x,weights=np.polynomial.legendre.leggauss(32)
    ws=.75*x+1.25;weights=.75*weights
    rows=[]
    for N in (17,33,65,129,257,1025,4097):
        value=0.;continuum=0.
        for w,weight in zip(ws,weights):
            u,wp=probe.solve_u(N,w);alpha=-(.5-u)/6
            v=0.;vc=0.
            for k in (N//2,-(N//2)):
                hn=old.H(m,q,N,[k,0,0],alpha)
                gn=m.GAMMA[0]@np.diag(q*probe.source(N,k+alpha*q))
                hc=old.H(m,q,N,[k,0,0],alpha,kind='continuum')
                gc=m.GAMMA[0]@np.diag(q)
                v+=diagonal_lower(hn,gn)/2
                vc+=diagonal_lower(hc,gc)/2
            value+=weight*v/(6*wp)
            continuum+=weight*vc/(6*wp)
        lo,_=probe.solve_u(N,.5);hi,_=probe.solve_u(N,2.)
        rows.append(dict(N=N,strip_width=(hi-lo)/6,
            time_smeared_variance_lower_bound=value,lower_bound_over_logN_squared=value/np.log(N)**2,
            continuum_diagonal_lower_same_strip=continuum))
    expected=12/old.ETA*float(weights@(ws*np.exp(-old.BETA*ws)/(1+np.exp(-old.BETA*ws))**2))
    assert all(b['time_smeared_variance_lower_bound']>a['time_smeared_variance_lower_bound']
               for a,b in zip(rows,rows[1:]))
    return dict(working_round=884,formal_round_still=883,cumulative_groups_still=3668,
         scope='Uniform classical flat hypercharge average; fixed smooth compact time test integral1. All original mass matrices. No full dynamical Gauss Gibbs assertion.',
         bound='sum n(E)(1-n(E)) abs(G_ii)^2 <= full time-smeared variance since fhat(0)=1',
         expected_largeN_coefficient_pending_proof=expected,rows=rows,
         full_smeared_variance_numerically_integrated=False,
         asymptotic_proof_completed=False,full_dynamical_gauge_model_refuted=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
