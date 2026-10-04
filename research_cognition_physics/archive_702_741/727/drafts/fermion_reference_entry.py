"""727 entry: a valid alternative reference and its inherited Weyl restriction."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_fermion_gauss_completion as matter
TARGET=HERE/'fermion_reference_entry_results.json'


def run():
    rng=np.random.default_rng(727)
    occupied=np.arange(0,64,2);empty=np.arange(1,64,2)
    def measure(h,delta):
        mean=float(np.trace(h[np.ix_(occupied,occupied)]).real)
        normal=float(np.sum(abs(h[np.ix_(empty,occupied)])**2))
        creation=float(np.sum(abs(delta[np.ix_(empty,empty)])**2)/2)
        annihilation=float(np.sum(abs(delta[np.ix_(occupied,occupied)])**2)/2)
        return dict(mean=mean,normal_variance=normal,pair_variance=creation+annihilation,
                    total_variance=normal+creation+annihilation)
    rows=[];det_error=0.
    for _ in range(5):
        phis=rng.normal(size=(2,5))*.35
        h=np.zeros((64,64),complex);delta=np.zeros_like(h)
        for v,p in enumerate(phis):
            hh,dd=matter.mass_matrices(p);h[v*32:(v+1)*32,v*32:(v+1)*32]=hh
            delta[v*32:(v+1)*32,v*32:(v+1)*32]=dd
        C=matter.gauge.group_exp(rng.normal(size=8),3)
        W=matter.gauge.group_exp(rng.normal(size=3),2)
        R=matter.representation(C,W,np.exp(.43j))
        det_error=max(det_error,float(abs(np.linalg.det(R[::2,::2])-1)))
        h[:32,32:]=.3*R;h[32:,:32]=.3*R.conj().T
        result=measure(h,delta);assert abs(result['mean'])<1e-12 and result['total_variance']<1e-24
        rows.append(result)
    assert det_error<1e-12
    # The same original internal content with604's explicit spin-direction edges.
    sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
    eta=np.r_[-np.ones(6),np.ones(3),np.ones(3),-np.ones(2),np.ones(1),np.ones(1)]
    eps=.7;directions=[]
    for axis,sigma in zip(('x','y','z'),(sx,sy,sz)):
        h=np.zeros((64,64),complex);J=-1j*np.kron(np.diag(eta),sigma)/(2*eps)
        h[:32,32:]=J;h[32:,:32]=J.conj().T
        result=measure(h,np.zeros_like(h));result['axis']=axis;directions.append(result)
    expected=16/eps**2
    assert abs(sum(r['total_variance'] for r in directions)-expected)<1e-12
    assert directions[0]['total_variance']>0 and directions[2]['total_variance']==0
    deps=('research_note_598.md','research_note_604.md','research_note_614.md','research_note_642.md',
          'research_note_709.md','research_note_720.md','research_note_726.md',
          'joint_fermion_gauss_completion.py')
    return dict(entry_round=727,new_formal_round=False,tests_run=2,failures=0,errors=0,
                result=dict(nodes=2,physical_modes=64,occupied_modes=32,
                    fixed_spin_determinant_Gauss_error=det_error,original_mass_and_scalar_link_rows=rows,
                    original_Weyl_edge_rows=directions,direction_sum_variance=float(expected),
                    empty_reference_failure_is_not_all_reference_failure=True,
                    state_not_claimed_ground_thermal_or_autonomously_prepared=True),
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
                scope='Entry application of known CAR filling and720 direction restriction. Complete one-spin filling is a Gauss reference with zero original mass and spin-scalar edge action. It fails transverse604 Weyl edges. No classification of all references or full chiral continuum claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
