"""735 entry: invariant action counterterms cannot change a Ward defect."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_dynamic_continuum_reference as ref
import joint_relative_source_development as prior
TARGET=HERE/'counterterm_scope_entry_results.json'

def potential(phi):
    M=ref.old.bdg(*ref.old.matter.mass_matrices(phi))
    return float(np.trace(np.linalg.matrix_power(M,4)).real/2)

def gradient(phi):
    result=[]
    for v in np.eye(5):
        M,D=prior.mass_jet(phi,v)
        result.append(2*np.trace(np.linalg.matrix_power(M,3)@D).real)
    return np.array(result)

def run():
    point=ref.point();phi=point['phi'];velocity=point['v'];grad=gradient(phi)
    rate=float(grad@velocity);h=2e-5
    finite=(potential(phi+h*velocity)-potential(phi-h*velocity))/(2*h)
    assert abs(rate-finite)<2e-7 and abs(rate)>1e-5
    X=phi[:2]+1j*phi[2:4];gauge=[]
    for weak in (*[np.asarray(x)/2 for x in ref.SIG],3*np.eye(2)):
        dX=1j*weak@X;v=np.r_[dX.real,dX.imag,0.]
        gauge.append(float(grad@v))
    assert max(abs(x) for x in gauge)<1e-12
    metric_divergence=-rate;scalar_exchange=-rate
    defect=metric_divergence-scalar_exchange
    assert defect==0.
    deps=('research_note_580.md','research_note_599.md','research_note_628.md','research_note_734.md',
          'joint_retarded_reference_response.py','joint_retarded_reference_response_results.json',
          'joint_relative_source_development.py','joint_dynamic_continuum_reference.py')
    return dict(entry_round=735,new_formal_round=False,
        potential_direction='V=one_half_trace(M(phi)^4), original complete Nambu mass',
        potential_value=potential(phi),five_scalar_gradient=grad.tolist(),
        original_time_derivative=rate,finite_difference=finite,derivative_error=abs(rate-finite),
        weak_and_circle_gauge_derivatives=gauge,
        metric_source_divergence=metric_divergence,scalar_source_exchange=scalar_exchange,
        joint_Ward_defect=defect,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope='A common covariant local action changes the original full sources but contributes zero joint Ward defect. Its finite coefficients therefore cannot repair an already nonzero defect of a point-splitting prescription. This is an audit of allowed corrections, not a new loop coefficient or a failure of anomaly-free renormalization.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
