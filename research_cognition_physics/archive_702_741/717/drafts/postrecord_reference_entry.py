"""717 entry: actual full-CAR pinching, energy and commutator moments."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_smooth_mode_contract as model
car=model.car
entry=model.entry
TARGET=HERE/'postrecord_reference_entry_results.json'


def inner(a,b):
    return sum(value.conjugate()*b.get(k,0) for k,value in a.items())


def subtract(a,b):
    return car.old.add(dict(a),b,-1)


def run():
    rng=np.random.default_rng(7171)
    points=np.array([[.4,-.3,.2,.1,.35],[-.2,.15,.1,-.25,.3],
                     [.22,.34,-.12,.08,-.27]])
    links=[model.group.sample(rng) for _ in range(3)]
    mass,d,hop=model.coefficients(points,links,np.array([.21,.27,.18]))
    h=mass+hop
    f,_=model.profile(np.array([.7,1.1,.5]),
                      np.array([[1.,.2j],[.7+.3j,.15],[-.3,.4j]]))
    g=np.zeros_like(f);g[62]=1.
    g-=f*np.vdot(f,g);g/=np.linalg.norm(g)
    dh,dd=model.difference_coeff(h,d,f)
    psi=(f+g)/np.sqrt(2)
    state={1<<int(i):psi[i] for i in np.flatnonzero(abs(psi)>1e-14)}
    k1=entry.n_apply(state,f)
    k0=subtract(state,k1)
    H=lambda s:car.quadratic(s,h,d)
    delta=lambda s:car.quadratic(s,dh,dd)
    comm=lambda s:subtract(H(entry.n_apply(s,f)),entry.n_apply(H(s),f))
    samples=[state,{0:1.},{(1<<3)|(1<<63):1.}]
    errors=[]
    for s in samples:
        errors.append(car.difference(comm(s),delta(entry.reflect(s,f))))
        errors.append(car.difference(entry.reflect(delta(s),f),
                                     {i:-v for i,v in delta(entry.reflect(s,f)).items()}))
    pre=float(inner(state,delta(state)).real)
    post=float(sum(inner(k,delta(k)).real for k in (k0,k1)))
    pre2=float(inner(delta(state),delta(state)).real)
    post2=float(sum(inner(delta(k),delta(k)).real for k in (k0,k1)))
    energy=float(inner(state,H(state)).real)
    energy_post=float(sum(inner(k,H(k)).real for k in (k0,k1)))
    probability=float(inner(k1,k1).real)
    drift=sum(inner(k,comm(k)) for k in (k0,k1))
    error=max(*errors,abs(post),abs(pre2-post2),abs(energy_post-energy-pre),abs(drift))
    assert error<1e-12 and abs(pre)>.001
    assert abs(sum(inner(k,k).real for k in (k0,k1))-1)<1e-13
    names=('research_note_624.md','research_note_633.md','research_note_716.md',
           'joint_smooth_mode_contract.py','round716_drafts/sterile_mode_entry.py')
    return dict(entry_round=717,new_formal_round=False,full_original_modes=96,
                pure_sterile_input_on_full_CAR=True,selected_probability=probability,
                original_pointwise_energy=energy,postrecord_pointwise_energy=energy_post,
                first_readout_energy_difference=pre,immediate_second_difference=post,
                pre_readout_D2=pre2,post_readout_D2=post2,
                postrecord_first_occupation_drift_abs=float(abs(drift)),
                full_sparse_CAR_identity_max_error=float(error),
                bosonic_marginal_preservation_proved_analytically=True,
                selected_branches_not_claimed_to_preserve_bosonic_marginals=True,
                subsequent_H_evolution_not_claimed_stationary=True,
                not_a_full_Gibbs_or_continuum_simulation=True,
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
