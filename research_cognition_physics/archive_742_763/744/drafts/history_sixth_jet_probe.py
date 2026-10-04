"""Follow actual native parity history to the next allowed time order."""
from pathlib import Path
from math import comb
import json
import numpy as np
import history_jet_probe as old

def more_jets(state):
    result=old.jets(state)
    result.append(old.H(result[-1]));result.append(old.H(result[-1]))
    return result

def measure(jet,order):
    u,wu=np.polynomial.laguerre.laggauss(order);v,wv=np.polynomial.hermite.hermgauss(order)
    r,s=np.meshgrid(np.sqrt(u/old.ALPHA),v/np.sqrt(old.ALPHA),indexing='ij')
    weight=wu[:,None]*wv[None,:]*u[:,None]/np.sqrt(1+(r*r+s*s)/6);weight/=weight.sum()
    A=np.sin(np.sqrt(2/(1+(r*r+s*s)/6))*s)**2/4
    values=[old.evaluate(p,r,s) for p in jet]
    integrand=-sum((-1)**k*comb(6,k)*old.inner(values[6-k],values[k]) for k in range(7))
    return dict(sixth=float(np.sum(weight*A*integrand).real),identity=float(np.sum(weight*integrand).real))

if __name__=='__main__':
    a=more_jets(0);b=more_jets((1<<30)|(1<<31))
    result=dict(terms_empty=[len(x) for x in a],terms_pair=[len(x) for x in b],rows=[])
    for n in (48,64,80):
        av=measure(a,n);bv=measure(b,n)
        result['rows'].append(dict(order=n,empty=av,pair=bv,difference_sixth=bv['sixth']-av['sixth']))
    with Path(__file__).with_name('history_sixth_jet_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))
