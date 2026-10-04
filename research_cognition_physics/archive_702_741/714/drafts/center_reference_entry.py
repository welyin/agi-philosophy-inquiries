"""714 entry: original center-related signal states and equal energy means."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_holonomy_readout_scale as loops
import joint_vertex_shared_evolution as car
old=loops.original
TARGET=HERE/'center_reference_entry_results.json'
CENTER=(np.eye(3),-np.eye(2),1.)


def scalar_action(g,p):
    x=p[:2]+1j*p[2:4];y=g[2]**3*g[1]@x
    return np.r_[y.real,y.imag,p[4]]


def run():
    rng=np.random.default_rng(7141);dims=(4,3,3)
    sites=list(product(*(range(n) for n in dims)))
    def nxt(v,mu):
        result=list(v);result[mu]=(result[mu]+1)%dims[mu];return tuple(result)
    links={(v,mu):old.sample(rng) for v in sites for mu in range(3)}
    changed={key:old.product(CENTER,g) if key[1]==0 and key[0][0]==0 else g for key,g in links.items()}
    def plaquette(data,v,mu,nu):
        return loops.multiply((data[v,mu],data[nxt(v,mu),nu],
                               old.inverse(data[nxt(v,nu),mu]),old.inverse(data[v,nu])))
    errors=[]
    for v in sites:
        for mu,nu in ((0,1),(0,2),(1,2)):
            a=plaquette(links,v,mu,nu);b=plaquette(changed,v,mu,nu)
            errors.append(float(np.linalg.norm(old.matter.representation(*a)-old.matter.representation(*b))))
    flip_errors=[]
    for y,z in product(range(dims[1]),range(dims[2])):
        a=loops.multiply(links[(x,y,z),0] for x in range(dims[0]))
        b=loops.multiply(changed[(x,y,z),0] for x in range(dims[0]))
        flip_errors.append(abs(loops.f(a)+loops.f(b)))
    scalar_errors=[]
    for _ in range(12):
        p=rng.normal(size=5)*.3;q=rng.normal(size=5)*.3;g=old.sample(rng)
        cq=scalar_action(CENTER,q)
        def mean(link):
            return .5*(old.original.distance_squared(p,scalar_action(link,q))+
                       old.original.distance_squared(p,scalar_action(link,cq)))
        scalar_errors.append(float(abs(mean(g)-mean(old.product(CENTER,g)))))
    assert max(errors+flip_errors+scalar_errors)<1e-12
    h,d=old.matter.mass_matrices(np.array([.4,-.3,.2,.1,.35]))
    vac={0:1.};bv=car.quadratic(vac,h,d)
    vacuum_mean=car.old.inner(vac,bv)
    nonstationary=float(car.old.inner(bv,bv).real)
    assert abs(vacuum_mean)<1e-14 and nonstationary>1e-5
    previous=json.loads((ARCHIVE/'joint_holonomy_readout_scale_results.json').read_text('utf8'))
    row=previous['actual_state_energy']['tilted_Haar_checks'][1]
    names=('research_note_598.md','research_note_617.md','research_note_623.md',
           'research_note_637.md','research_note_713.md','joint_holonomy_readout_scale.py')
    return dict(entry_round=714,new_formal_round=False,periodic_grid=list(dims),original_modes_per_node=32,
        full_plaquettes_checked=len(errors),maximum_original_representation_plaquette_error=max(errors),
        loop_sign_error=max(flip_errors),original_scalar_edge_sign_average_error=max(scalar_errors),
        original_CAR_vacuum_mass_mean_real=float(vacuum_mean.real),
        nonzero_mass_action_norm_squared=nonstationary,
        inherited_tilt=3.,opposite_means=[row['signal_mean'],-row['signal_mean']],
        eta_point_six_record_TV=.6*row['signal_mean'],
        complete_H_mean_equality_is_analytic=True,geometry_source_mean_equality_is_analytic=True,
        actual_control_implementation_or_zero_preparation_work_not_proved=True,
        reference_is_explicit_normal_Gauss_state_not_Gibbs=True,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
