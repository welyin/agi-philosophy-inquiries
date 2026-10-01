"""657 entry: finite necessary tests of a declared two-slice reflection kernel.

E is reflected as a scalar, interchanging the two time slices with no internal
twist. Finite positive Gram samples do not prove reflection positivity, a
transfer semigroup, full S9 integration, or identification with physical CAR.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_spatial_auxiliary_geometry as model


def run():
    f=model.frame(1);rng=np.random.default_rng(657);count=10;rows=[]
    base=np.zeros((6,10));base[:,0]=1
    pf0=model.old.pfaffian(model.auxiliary(f,base))
    for mode,width in [('plane',.1),('plane',.7),('plane',2.),('S9',.1),('S9',.5),('S9',1.)]:
        e=np.zeros((count,3,10))
        if mode=='plane':
            th=rng.uniform(-width,width,(count,3));e[:,:,0]=np.cos(th);e[:,:,1]=np.sin(th)
        else:
            e=np.eye(10)[0]+width*rng.normal(size=(count,3,10));e/=np.linalg.norm(e,axis=2)[:,:,None]
        matrix=np.zeros((count,count),complex)
        for i in range(count):
            for j in range(count):
                # Compute both orientations; do not impose Hermiticity in input.
                matrix[i,j]=model.old.pfaffian(model.auxiliary(f,np.vstack((e[i],e[j]))))/pf0
        herm=float(np.max(abs(matrix-matrix.conj().T)))
        assert herm<3e-13 and min(matrix.diagonal().real)>0
        scale=np.sqrt(np.outer(matrix.diagonal().real,matrix.diagonal().real));gram=matrix/scale
        eig=np.linalg.eigvalsh((gram+gram.conj().T)/2)
        assert eig[0]>-3e-11
        rows.append(dict(configuration_family=mode,width=width,samples=count,
            Hermiticity_residual=herm,minimum_diagonal=float(min(matrix.diagonal().real)),
            normalized_Gram_eigenvalues=eig.tolist(),
            maximum_offdiagonal=float(np.max(abs(gram-np.eye(count))))))
    return dict(date='2026-10-02',status='657 entry; incomplete',rows=rows,
        reflection='swap t=0 and t=1; E transforms as scalar; no internal twist',
        original_background='nx=3, nt=2 antiperiodic, other two spatial directions length1, m0=1, kappa=1, identity gauge',
        dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in
            ('joint_spatial_auxiliary_geometry.py','joint_subgroup_measure_source.py','research_note_656.md')},
        scope='Only six finite Gram tests using actual signed Pfaffian. No proof of full-kernel positivity, arbitrary time extent, original S9 integration or common physical state.')


if __name__=='__main__':
    result=run()
    with (HERE/'reflection_kernel_probe_results.json').open('x',encoding='utf8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=result['status'],cases=len(result['rows']),
        minimum_sampled_eigenvalue=min(r['normalized_Gram_eigenvalues'][0] for r in result['rows']))))
