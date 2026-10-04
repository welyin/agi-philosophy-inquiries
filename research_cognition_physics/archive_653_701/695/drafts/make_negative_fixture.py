"""Freeze rational full-group links near a resolved original-model sign witness."""
import json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as b
GRID=4096
def hermitian_input(u):
    a=-1j*(u-np.eye(len(u)))@np.linalg.inv(u+np.eye(len(u)))
    a=(a+a.conj().T)/2
    return dict(real=np.rint(a.real*GRID).astype(int).tolist(),imag=np.rint(a.imag*GRID).astype(int).tolist())
links=[]
for mu in (0,1):
    row=[]
    for i in range(4):
        c,w,z=b.prior.group(67351+8*mu+i,.9)
        row.append(dict(color=hermitian_input(c),weak=hermitian_input(w),
                        u1_cayley=int(round(float(z.imag/(1+z.real))*GRID))))
    links.append(row)
rng=np.random.default_rng(69500);e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
chart=np.rint(e[:,1:]/(1+e[:,:1])*GRID).astype(int).tolist()
out=dict(grid=GRID,links=links,negative_candidate_E_stereographic=chart,
         defining_data_are_the_integers_not_float_seed=True,
         exploratory_origin=dict(full_group_seed=67351,gauge_scale=.9,E_seed=69500))
target=HERE/'negative_auxiliary_fixture.json'
with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,indent=2)+'\n')
print('Rational695 fixture saved')
