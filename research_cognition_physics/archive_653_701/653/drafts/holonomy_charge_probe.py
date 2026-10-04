"""653 entry probe: exact Fourier support of the actual 615 conditional weight.

No claim of complete Hamiltonian reconstruction or new particles. This tests
the specific identity with a fixed original 32-mode charge-twisted trace.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_subgroup_measure_source as old


def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():out[i+j]=out.get(i+j,F(0))+x*y
    return {k:v for k,v in out.items() if v}


def power(a,n):
    out={0:F(1)}
    for _ in range(n):out=mul(out,a)
    return out


def cos_square(q):
    if q==0:return {0:F(1)}
    return {-abs(q):F(1,4),0:F(1,2),abs(q):F(1,4)}


def evaluate(poly,theta):return sum(float(v)*np.exp(1j*k*theta) for k,v in poly.items())


def run():
    physical={0:F(1)}
    for q in old.Q:physical=mul(physical,cos_square(int(q)))
    auxiliary={}
    for k,c in enumerate(old.COEFFICIENTS):
        term=mul(power(cos_square(2),k),power(cos_square(3),8-k))
        for charge,v in term.items():auxiliary[charge]=auxiliary.get(charge,F(0))+c*v
    full=mul(physical,auxiliary)
    assert all(v>0 for poly in (physical,auxiliary,full) for v in poly.values())
    assert all(sum(poly.values())==1 for poly in (physical,auxiliary,full))
    assert sum(int(q) for q in old.Q)==0
    canonical_max=2*sum(max(int(q),0) for q in old.Q)
    assert canonical_max==36 and max(physical)==36 and max(auxiliary)==24 and max(full)==60
    assert full[60]==F(1,55*4**23)
    errors=[]
    for theta in (0.,.07,.17,.31,.67,1.2,2.4):
        errors.extend((abs(evaluate(physical,theta)-old.physical(theta)),
                       abs(evaluate(auxiliary,theta)-old.measure(theta)[0]),
                       abs(evaluate(full,theta)-old.full(theta))))
    assert max(errors)<2e-14
    # Explicit exact 32-mode occupation counting, independent of cos products.
    counts={0:1}
    for q in old.Q:
        for _ in range(2):
            nxt={}
            for n,c in counts.items():
                nxt[n]=nxt.get(n,0)+c;nxt[n+int(q)]=nxt.get(n+int(q),0)+c
            counts=nxt
    assert {q:F(c,2**32) for q,c in counts.items()}==physical
    variance={name:str(sum(F(q*q)*p for q,p in poly.items()))
              for name,poly in [('physical',physical),('auxiliary',auxiliary),('full',full)]}
    assert variance=={'physical':'60','auxiliary':'24','full':'84'}
    outside=sum(p for q,p in full.items() if abs(q)>canonical_max)
    return dict(date='2026-10-02',status='entry probe; round653 incomplete',
        fixed_original_CAR_modes=32,original_charge_range=[-36,36],
        auxiliary_fourier_range=[-24,24],full_fourier_range=[-60,60],
        full_highest_coefficient=str(full[60]),outside_original_charge_weight=str(outside),
        outside_original_charge_weight_float=float(outside),charge_variances=variance,
        numerical_reconstruction_error=float(max(errors)),
        positive_auxiliary_charge_support_size=len(auxiliary),
        coefficients={name:{str(k):str(v) for k,v in sorted(poly.items())}
                      for name,poly in [('physical',physical),('auxiliary',auxiliary),('full',full)]},
        dependency_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('research_note_615.md','joint_subgroup_measure_source.py')},
        limitation='Conditional fixed-holonomy fermion factor, before gauge/boson integration. No full Gauss trace or continuum impossibility claim.')


if __name__=='__main__':
    result=run()
    with (HERE/'holonomy_charge_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('coefficients','dependency_hashes')}))
