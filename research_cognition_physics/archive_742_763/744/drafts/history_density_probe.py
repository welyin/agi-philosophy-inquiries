"""Resolve native sixth-jet difference as a radial polynomial density."""
from math import comb
from pathlib import Path
import json
import history_sixth_jet_probe as old

def radial(poly):
    out={}
    for (state,p),v in poly.items():
        if any(p[1:4]):continue
        out.setdefault(state,{})[p[0],p[4]]=v
    return out

def density(jets,n):
    radial_jets=[radial(p) for p in jets];out={}
    for k in range(n+1):
        a,b=radial_jets[n-k],radial_jets[k];factor=1j**n*(-1)**k*comb(n,k)
        for state,terms in a.items():
            if state not in b:continue
            for (r,s),v in terms.items():
                for (rr,ss),z in b[state].items():
                    key=(r+rr,s+ss);out[key]=out.get(key,0)+factor*v.conjugate()*z
    return out

if __name__=='__main__':
    a=old.more_jets(0);b=old.more_jets((1<<30)|(1<<31));result={}
    for n in (4,6):
        pa,pb=density(a,n),density(b,n)
        delta={key:pb.get(key,0)-pa.get(key,0) for key in set(pa)|set(pb)}
        kept={key:z for key,z in delta.items() if abs(z)>1e-8}
        result[str(n)]=dict(coefficients=[dict(r_power=k[0],s_power=k[1],real=z.real,imag=z.imag) for k,z in sorted(kept.items())],
            discarded_max=max((abs(z) for k,z in delta.items() if k not in kept),default=0))
    with Path(__file__).with_name('history_density_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))
