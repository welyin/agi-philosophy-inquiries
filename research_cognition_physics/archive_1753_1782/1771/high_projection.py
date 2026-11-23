"""Evaluate the independent colour/weak tensor bound at the matched unified scale."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def calculate():
    m = runpy.run_path(str(HERE/'check.py'))
    i = runpy.run_path(str(HERE/'independent_check.py'))
    cfg = m['CFG']
    q,mt,g5,a3,aq = m['matching'](m['mass_spectrum'](True))
    y,pt,_ = m['charged_point'](a3)
    hi = m['run_yukawa'](y,cfg['flavour_scale_GeV'],q,aq,pt['flow']['beta'],cfg['step_high']/2)
    up,d,e = hi[0,:,:3],hi[1,:,:3],hi[2,:,:3].T
    ul,sv,_ = np.linalg.svd(up)
    symmetric = (ul*sv)@ul.T
    np.testing.assert_allclose(symmetric@symmetric.conj().T,up@up.conj().T,atol=1e-14)
    sym_distance = np.linalg.norm((up-up.T)/2)
    assert sym_distance>1e-10
    joint,dim = i['minimum'](symmetric,e,d,np.eye(3),np.eye(3))
    weights = np.array([sum(np.delete(sv,j)**2) for j in range(3)])
    weighted = np.sqrt(weights)[:,None]*(ul.conj().T@(d-e))
    np.testing.assert_allclose(joint, math.sqrt(2/7)*np.linalg.norm(weighted),rtol=3e-11)
    eps=q/(5*g5*cfg['cutoff_for_adjoint_expansion_GeV'])
    perturbation=25*eps**2*cfg['aggregate_quadratic_coefficient_budget']*np.sqrt(sum(weights))
    protected=math.sqrt(2/7)*max(0,np.linalg.norm(weighted)-perturbation)
    assert dim==3 and 0<protected<joint
    return dict(Q_GeV=q,MT_GeV=mt,distance_to_symmetric_up_in_frozen_frame=sym_distance,
                relative_up_frame_distance=sym_distance/np.linalg.norm(up),
                physical_tensor_norm_times_MT=joint,physical_tensor_norm_per_GeV=joint/mt,
                bounded_mixed_insertion_floor_per_GeV=protected/mt)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=json.loads(json.dumps(calculate()),parse_float=lambda v:float(f'{float(v):.10g}'))
    path=HERE/'high_projection_results.json'
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1771 matched scale: symmetric embedding and irreducible joint d=5 floor checked")
    else:
        with path.open('x',encoding='utf8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
