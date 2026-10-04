"""708 entry: inherited hyperboloid linearity ties nested centroid and collective mass.
An algebraic entry, not a new formal round. Actual normal-state insufficiency
of707 radial trace is analytic; samples are only coefficient/coordinate checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_geodesic_spatial_block as old
TARGET=HERE/'nested_mass_entry_results.json'
R=old.R

def merge(a,b):
    na,A,x=a;nb,B,y=b
    S=A*old.hyper(x)+B*old.hyper(y)
    C=np.sqrt(-S@old.ETA@S)/R
    return na+nb,float(C),old.klein(S)

def chain(items):
    value=items[0]
    for p in items[1:]:value=merge(value,p)
    return value

def run():
    rng=np.random.default_rng(708);points=.55*rng.normal(size=(4,5))
    leaves=[(1,1.,p) for p in points]
    direct=sum(old.hyper(p) for p in points)
    A=float(np.sqrt(-direct@old.ETA@direct)/R);m=old.klein(direct)
    trees=[chain(leaves),merge(merge(leaves[0],leaves[1]),merge(leaves[2],leaves[3])),
           merge(merge(leaves[0],leaves[2]),merge(leaves[1],leaves[3])),
           chain(leaves[::-1])]
    errors=[max(abs(t[1]-A),np.linalg.norm(t[2]-m)) for t in trees]
    hm,dm=old.matter.mass_matrices(m)
    masses=[old.matter.mass_matrices(p) for p in points]
    hbar=sum(h for h,d in masses)/4;dbar=sum(d for h,d in masses)/4
    masserr=max(np.linalg.norm(hbar-A/4*hm),np.linalg.norm(dbar-A/4*dm))
    # A real orthogonal four-site CAR transform with a normalized uniform first row.
    O=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])/2
    full=np.zeros((128,128),complex);delta=np.zeros_like(full)
    for i,(h,d) in enumerate(masses):
        full[32*i:32*i+32,32*i:32*i+32]=h
        delta[32*i:32*i+32,32*i:32*i+32]=d
    U=np.kron(O,np.eye(32))
    hc=U@full@U.T;dc=U@delta@U.T
    modeerr=max(np.linalg.norm(hc[:32,:32]-hbar),np.linalg.norm(dc[:32,:32]-dbar))
    # Same707 midpoint m0 and internal orientation u, different discarded radius.
    m0=np.array([.22,-.11,.08,.14,.36]);u=np.array([.4,.3,-.2,.1,.8]);u/=np.linalg.norm(u)
    third=np.array([-.1,.2,.04,-.13,-.27]);witness=[]
    for r in (.3,1.3):
        x,y=old.endpoints(m0,r*u);pair=merge((1,1.,x),(1,1.,y))
        assert np.linalg.norm(pair[2]-m0)<1e-14
        assert abs(pair[1]-2*np.cosh(r/R))<1e-14
        triple=merge(pair,(1,1.,third))
        witness.append(dict(discarded_radius=r,pair_amplitude=pair[1],
            pair_inverse_cosh=float(1/np.cosh(r/R)),third_merge_midpoint=triple[2].tolist(),
            third_merge_total_amplitude=triple[1]))
    future=np.linalg.norm(np.array(witness[0]['third_merge_midpoint'])-np.array(witness[1]['third_merge_midpoint']))
    # Dropping amplitude after each pair leads to pairing-dependent nested midpoints.
    bare01=old.midpoint(old.midpoint(points[0],points[1]),old.midpoint(points[2],points[3]))
    bare02=old.midpoint(old.midpoint(points[0],points[2]),old.midpoint(points[1],points[3]))
    naive=float(np.linalg.norm(bare01-bare02))
    assert max(*errors,masserr,modeerr)<1e-12 and future>1e-3 and naive>1e-4
    deps=('research_note_707.md','joint_geodesic_spatial_block.py','research_note_642.md',
          'research_note_598.md')
    return dict(date='2026-10-02',entry_round=708,latest_formal_round=707,new_formal_round=False,
        aggregate=dict(count=4,amplitude=A,mass_factor=A/4,tree_errors=errors,
            collective_Dirac_Majorana_error=float(masserr),full128mode_uniform_block_error=float(modeerr),
            retained_to_difference_mass_norm=float(np.sqrt(np.linalg.norm(hc[:32,32:])**2+np.linalg.norm(dc[:32,32:])**2))),
        same707_coarse_coordinates=witness,future_midpoint_difference=float(future),
        discarding_amplitude_pairing_difference=naive,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(fixed_equal_site_weights=True,amplitude_is_existing_relative_information=True,
            no_complete_H_or_instrument_history_closure=True,no_new_formal_scientific_round=True),
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False))
