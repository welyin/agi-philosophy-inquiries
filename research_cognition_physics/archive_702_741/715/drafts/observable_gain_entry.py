"""715 entry: original loop readout under affine and bounded nonlinear gain."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_holonomy_readout_scale as old
TARGET=HERE/'observable_gain_entry_results.json'


def run():
    rng=np.random.default_rng(7151)
    n=4;m=3;eta=.6;eps=.73;sigma=.12;step=1e-5
    data=[[old.original.sample(rng) for _ in range(n)] for _ in range(m)]
    def mean(links):return np.mean([old.f(old.multiply(row)) for row in links])
    F=mean(data)
    gradient=n/eps*np.exp(-2*sigma)/m**2*sum(old.gamma(old.multiply(row)) for row in data)
    checks=[];scale=[]
    for z in (1.,2.,8.,32.):
        def amps(links):
            response=eta*np.tanh(z*mean(links))
            return np.sqrt(np.array([1+response,1-response])/2)
        h=eta*np.tanh(z*F)
        hp=eta*z/(np.cosh(z*F)**2)
        exact=gradient*hp*hp/(4*(1-h*h))
        finite=0.
        for l in range(m):
            for e in range(n):
                for kind,dim,b in (('weak',3,old.BW),('abelian',1,old.B0)):
                    for a in range(dim):
                        plus=[row.copy() for row in data];minus=[row.copy() for row in data]
                        plus[l][e]=old.shift(data[l][e],kind,a,step)
                        minus[l][e]=old.shift(data[l][e],kind,a,-step)
                        derivative=(amps(plus)-amps(minus))/(2*step)
                        finite+=b/eps*np.exp(-2*sigma)*float(derivative@derivative)
        assert abs(exact-finite)<1e-9
        checks.append(dict(gain=z,injection=exact,direct_Kraus_derivative_error=abs(exact-finite)))
        norm=eta**2*old.CSTAR/(4*.25)*z*z
        scale.append(dict(gain=z,affine_minimum_effect=(1-eta*z)/2,
                          bounded_nonlinear_minimum_effect=(1-eta*np.tanh(z))/2,
                          fixed_area_point_two_five_length_one_exact_norm=norm))
    # Algebraic uniform bound of the nonlinear derivative factor.
    t=np.linspace(0,1,1001)
    assert np.max((1-t)**2/(1-eta**2*t))<=1+1e-14
    maximizer=(np.eye(3),np.diag([1j,-1j]),1.) if old.BW/4>=9*old.B0 else (
        np.eye(3),np.eye(2),np.exp(1j*np.pi/6))
    assert abs(old.f(maximizer))<1e-13
    assert abs(old.gamma(maximizer)-old.CSTAR)<1e-13
    names=('research_note_592.md','research_note_639.md','research_note_701.md',
           'research_note_713.md','research_note_714.md','joint_holonomy_readout_scale.py')
    return dict(entry_round=715,new_formal_round=False,original_loop_mean=float(F),
                actual_original_Kraus_checks=checks,conditional_gain_family=scale,
                norm_is_worst_state_not_original_Gibbs_expectation=True,
                divergent_gain_is_hypothesis_not_original_renormalization_result=True,
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
