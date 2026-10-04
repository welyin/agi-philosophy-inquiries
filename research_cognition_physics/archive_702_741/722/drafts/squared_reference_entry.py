"""722 entry: test squared-reference flow before claiming energy preservation."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_velocity_record_energy as old
TARGET=HERE/'squared_reference_entry_results.json'


def run():
    # For unnormalised tent a=1-|eta|: jumps of (exp(i eta²/4)a)'.
    n=256;t,w=np.polynomial.legendre.leggauss(n)
    nodes=np.r_[(t-1)/2,(t+1)/2];weights=np.r_[w/2,w/2]
    a=1-abs(nodes);f=np.exp(.25j*nodes**2)*a
    rows=[]
    threshold=(33/4)/np.sin(.25)
    for k in (40.,80.,160.):
        actual=np.dot(weights,f*np.exp(1j*k*nodes))  # xi negative
        leading=-(2*np.exp(.25j)*np.cos(k)-2)/k**2
        remainder=abs(actual-leading)
        lower=np.sin(.25)/k**2
        assert k>=threshold and remainder<=(33/4)/k**3+1e-13
        assert abs(actual)>=lower
        # Normalized tent propagated at unit flow-Schrodinger parameter.
        amplitude=np.sqrt(1.5)*abs(actual)/np.sqrt(4*np.pi)
        xi=-2*k;cmin=old.M/(1+.5**2/6)
        rows.append(dict(k=k,abs_Fourier=float(abs(actual)),certified_asymptotic_lower=float(lower),
                         remainder_times_k_cubed=float(remainder*k**3),
                         log_weighted_density=float(-2*cmin*xi+2*np.log(amplitude))))
    errs=[]
    for p in (np.array([.3,.4,.1,-.2,.5]),np.array([1.2,.7,-.3,.2,-.4])):
        F=old.original.F(p);x=p/np.sqrt(F);R=np.linalg.norm(x[:4]);z=x[4]
        A=1+z*z/6
        # xi=(w/M)(A log R+R²/12), w=1; X_T xi=1.
        derivative=((A/R+R/6)/old.M)*(F*R)
        errs.append(abs(derivative-1))
    assert max(errs)<1e-14
    names=('joint_velocity_record_energy.py','research_note_721.md','research_note_652.md')
    return dict(entry_round=722,formal_round=False,rows=rows,
                Fourier_remainder_constant=33/4,lower_threshold=float(threshold),
                original_flow_coordinate_error=max(errs),
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names},
                scope='Exact scalar flow-coordinate and tent Fourier tail. Full Gauss-state lift, original Hardy/form estimate and squared-reference instrument conclusion require next formal audit.')


if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))

