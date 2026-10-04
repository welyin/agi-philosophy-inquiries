"""711 entry: finite raw gauge samples do not determine continuum winding.

The continuum ball, smooth map and orientation are explicit comparison inputs.
This is not an instanton simulation or a new numbered research round.
"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_fermion_gauss_completion as matter
TARGET=HERE/'topology_sampling_entry_results.json'
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
CENTER=np.array([.37,.43,.59]);RADIUS=.17


def profile(s):
    s=np.asarray(s);a=np.exp(-1/s);b=np.exp(-1/(1-s))
    h=a/(a+b);dh=h*(1-h)*(1/s**2+1/(1-s)**2)
    return np.pi*(1-h),-np.pi*dh


def gauge_map(x):
    offset=np.asarray(x)-CENTER;r=float(np.linalg.norm(offset))
    if r>=RADIUS:return np.eye(2,dtype=complex)
    if r==0:return -np.eye(2,dtype=complex)
    f,_=profile(r/RADIUS)
    return np.cos(f)*np.eye(2)+1j*np.sin(f)*np.einsum('i,ijk->jk',offset/r,PAULI)


def run():
    vertices=[np.array(v,float) for v in product((0,1),repeat=3)]
    pairs=[(i,j) for i in range(8) for j in range(i+1,8) if np.sum(abs(vertices[i]-vertices[j]))==1]
    matrices=[gauge_map(x) for x in vertices]
    node_error=max(float(np.linalg.norm(g-np.eye(2))) for g in matrices)
    links=[matrices[i].conj().T@matrices[j] for i,j in pairs]
    link_error=max(float(np.linalg.norm(u-np.eye(2))) for u in links)
    original_error=max(float(np.linalg.norm(matter.representation(np.eye(3),u,1)-np.eye(32))) for u in links)
    assert node_error==link_error==original_error==0
    inside=matter.representation(np.eye(3),gauge_map(CENTER),1)
    contrast=float(np.linalg.norm(inside-np.eye(32)));assert contrast==8
    windings=[]
    for n in (48,96,144):
        x,w=np.polynomial.legendre.leggauss(n);s=(x+1)/2
        f,df=profile(s);integral=float(np.sum(w/2*(-2/np.pi)*df*np.sin(f)**2))
        windings.append(dict(nodes=n,winding=integral,error_from_exact_one=abs(integral-1)))
    assert windings[-1]['error_from_exact_one']<2e-12
    # A naive radial contraction f -> t f is discontinuous at its center for
    # 0<t<1; opposite-direction limits differ by this exact Frobenius norm.
    gaps=[dict(t=t,opposite_direction_limit_gap=float(2*np.sqrt(2)*abs(np.sin(np.pi*t)))) for t in (.25,.5,.75)]
    dependencies=('research_note_598.md','research_note_629.md','research_note_673.md','research_note_710.md',
                  'joint_fermion_gauss_completion.py')
    return dict(entry_round=711,new_formal_round=False,finite_nodes=8,finite_edges=len(pairs),
        empty_ball_radius=RADIUS,closest_node_distance=min(float(np.linalg.norm(x-CENTER)) for x in vertices),
        node_sample_error=node_error,pure_gauge_link_error=link_error,original_32_mode_link_error=original_error,
        hidden_center_representation_contrast=contrast,winding_quadrature=windings,naive_contraction_gaps=gaps,
        exact_winding_proof='-(2/pi) [f/2-sin(2f)/4]_{pi}^{0}=1; orientation declared',
        scope='Finite raw spatial gauge samples versus continuum pi3 class; no four-dimensional instanton charge, continuum no-go or new physical configuration sector inferred.',
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in dependencies})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
