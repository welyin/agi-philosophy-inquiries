"""716 entry: original64-mode CAR projection and its exact active subspace."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_vertex_shared_evolution as car
import joint_region_energy_gluing as group
TARGET=HERE/'sterile_mode_entry_results.json'


def c_apply(state,f,creation=False):
    out={}
    for i in np.flatnonzero(abs(f)>1e-14):
        term=car.old.create(state,int(i)) if creation else car.old.annihilate(state,int(i))
        out=car.old.add(out,term,f[i] if creation else f[i].conjugate())
    return out


def n_apply(state,f):
    return c_apply(c_apply(state,f),f,True)


def reflect(state,f):
    # The inherited sparse addition mutates its first argument.
    return car.old.add(dict(state),n_apply(state,f),-2)


def run():
    points=[np.array([.4,-.3,.2,.1,.35]),np.array([-.2,.15,.1,-.25,.3])]
    h=np.zeros((64,64),complex);d=np.zeros_like(h)
    for v,p in enumerate(points):
        sl=slice(32*v,32*v+32)
        h[sl,sl],d[sl,sl]=car.matter.mass_matrices(p)
    rng=np.random.default_rng(7161);g=group.sample(rng)
    t=.23*np.exp(-.24)
    hop=t*car.matter.representation(*g)
    h[:32,32:]=hop;h[32:,:32]=hop.conj().T
    f=np.zeros(64,complex);f[30]=np.sqrt(.6);f[62]=np.sqrt(.4)*np.exp(.4j)
    P=np.outer(f,f.conj());Q=np.eye(64)-P;U=np.eye(64)-2*P
    dh=(U@h@U-h)/2;dd=(U@d@U.T-d)/2
    dh2=-(P@h@Q+Q@h@P);dd2=-(P@d@Q.T+Q@d@P.T)
    coeff_error=float(max(np.linalg.norm(dh-dh2),np.linalg.norm(dd-dd2)))
    states=[{0:1.},{1<<30:1.},{(1<<4)|(1<<63):1.},
            car.old.add({1<<30:1/np.sqrt(2)},{1<<62:1j/np.sqrt(2)})]
    errors=[]
    for state in states:
        rbr=reflect(car.quadratic(reflect(state,f),h,d),f)
        expected=car.old.add(rbr,car.quadratic(state,h,d),-1)
        actual={k:2*v for k,v in car.quadratic(state,dh,dd).items()}
        errors.append(car.difference(expected,actual))
        errors.append(car.difference(n_apply(n_apply(state,f),f),n_apply(state,f)))
    # The exact coefficient support is spanned by f, Qhf and Q Delta conjugate(f).
    span=np.column_stack((f,Q@h@f,Q@d@f.conj()))
    v,s,_=np.linalg.svd(span,full_matrices=False);v=v[:,s>1e-12]
    hs=v.conj().T@dh@v;ds=v.conj().T@dd@v.conj()
    support_error=float(max(np.linalg.norm(dh-v@hs@v.conj().T),
                            np.linalg.norm(dd-v@ds@v.T)))
    n=v.shape[1];B=np.zeros((2**n,2**n),complex)
    for col in range(2**n):
        for row,value in car.quadratic({col:1.},hs,ds).items():B[row,col]=value
    norm=float(np.max(abs(np.linalg.eigvalsh(B))))
    bound=float(np.linalg.norm(Q@h@f)+np.linalg.norm(Q@d@f.conj()))
    R=np.zeros_like(h)
    R[:32,:32]=car.matter.representation(*group.sample(rng))
    R[32:,32:]=car.matter.representation(*group.sample(rng))
    gauss=float(np.linalg.norm(R@f-f))
    if max(coeff_error,support_error,gauss,*errors)>=1e-12:
        diagnostic=dict(coefficient_error=coeff_error,support_error=support_error,
                        Gauss_error=gauss,CAR_errors=errors)
        if not (HERE/'first_action_diagnostic.json').exists():
            with (HERE/'first_action_diagnostic.json').open('x',encoding='utf8') as out:
                json.dump(diagnostic,out,indent=2)
        print(diagnostic)
    assert max(coeff_error,support_error,gauss,*errors)<1e-12
    assert norm<=bound+1e-12 and n<=3
    names=('research_note_598.md','research_note_633.md','research_note_666.md',
           'research_note_667.md','research_note_715.md','joint_vertex_shared_evolution.py')
    return dict(entry_round=716,new_formal_round=False,total_original_modes=64,
                actual_sterile_mode_indices=[30,62],CAR_action_error=max(errors),
                exact_coefficient_identity_error=coeff_error,Gauss_mode_error=gauss,
                active_single_particle_dimension=int(n),exact_active_Fock_dimension=2**n,
                support_reconstruction_error=support_error,full_CAR_operator_norm=norm,
                original_coefficient_norm_bound=bound,
                eigenvalues_of_actual_active_difference=np.linalg.eigvalsh(B).tolist(),
                normal_and_pairing_ranks=[int(np.linalg.matrix_rank(dh,tol=1e-12)),
                                          int(np.linalg.matrix_rank(dd,tol=1e-12))],
                comparison_does_not_compute_original_Gibbs=True,
                no_continuum_propagation_or_autonomous_device_claim=True,
                dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
