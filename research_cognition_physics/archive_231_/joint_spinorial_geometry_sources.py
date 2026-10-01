"""663: original full-mass continuum symbol, spin lift and metric sources.

Coefficient identities on a declared continuum background, not the missing
full lattice/chiral/Gauss limit or quantized Einstein constraints.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_chiral_source_matching as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_spinorial_geometry_sources_results.json'
GAMMA=old.kinetic_matrices()
MASS=old.mass()


def norm(a):return float(np.linalg.norm(a,2))


def exp_sym(a):
    w,v=np.linalg.eigh(a);return (v*np.exp(w))@v.T


def sqrt_spd(g):
    w,v=np.linalg.eigh(g);assert min(w)>0
    return (v*np.sqrt(w))@v.T


def rotation_generator(a):
    assert norm(a+a.T)<2e-12
    return sum(.25*a[i,j]*(GAMMA[i]@GAMMA[j]) for i in range(3) for j in range(3))


def spin(r):
    assert max(norm(r.T@r-np.eye(3)),abs(np.linalg.det(r)-1))<2e-12
    cosine=np.clip((np.trace(r)-1)/2,-1,1);theta=float(np.arccos(cosine))
    if theta<1e-7:
        a=.5*(r-r.T)*(1+theta**2/6)
    else:
        assert theta<2.5  # Local continuous branch only, no global spin assertion.
        a=theta/(2*np.sin(theta))*(r-r.T)
    s=rotation_generator(a);e,v=np.linalg.eigh(1j*s)
    u=(v*np.exp(-1j*e))@v.conj().T
    assert norm(u.conj().T@u-np.eye(64))<2e-12
    return u


def geometric_pullback(g,f):
    gp=f.T@g@f;e=sqrt_spd(g);ep=sqrt_spd(gp)
    r=e@f@np.linalg.inv(ep)
    return gp,r,spin(r).conj().T


def clifford(g,p):
    covector=np.linalg.solve(sqrt_spd(g).T,p)
    return sum(x*c for x,c in zip(covector,GAMMA))


def local_symbols(g,p,n,m,hn,hm,N=1.2,M=.7):
    cp=clifford(g,p);cn=clifford(g,n);cm=clifford(g,m)
    h=N*(cp+MASS)-.5j*cn
    inverse=np.linalg.inv(g);v=inverse@(N*m-M*n)
    div=N*np.trace(inverse@hm)-M*np.trace(inverse@hn)
    bracket=(-v@p+.5j*div)*np.eye(64)+.25j*(cn@cm-cm@cn)
    return h,bracket


def full_source_and_frame_check():
    g=exp_sym(np.array([[.25,.11,-.04],[.11,-.18,.07],[-.04,.07,.09]]))
    a=np.array([[.2,.14,-.06],[.14,-.12,.03],[-.06,.03,.08]])
    b=np.array([[-.1,.08,.02],[.08,.18,-.05],[.02,-.05,.03]])
    f=exp_sym(a)@exp_sym(b)
    gp,r,u=geometric_pullback(g,f)
    p=np.array([.41,-.28,.33]);n=np.array([.2,-.13,.07]);m=np.array([-.11,.31,.16])
    hn=np.array([[.13,.04,-.02],[.04,-.1,.05],[-.02,.05,.06]])
    hm=np.array([[-.07,.03,.01],[.03,.09,-.04],[.01,-.04,.15]])
    h,c=local_symbols(g,p,n,m,hn,hm)
    hp,cp=local_symbols(gp,f.T@p,f.T@n,f.T@m,f.T@hn@f,f.T@hm@f)
    errors=dict(mass=norm(u@MASS@u.conj().T-MASS),
        local_lapse_symbol=norm(hp-u@h@u.conj().T),
        complete_local_bracket=norm(cp-u@c@u.conj().T))
    assert max(errors.values())<2e-12
    # Actual same geometric parameter in old and new spinor identifications.
    dg=np.array([[.27,-.19,.09],[-.19,-.16,.14],[.09,.14,.06]])
    step=2e-5
    def family(t):
        gt=g+t*dg;gtp,_,ut=geometric_pullback(gt,f)
        ht=clifford(gt,p)+MASS;htp=clifford(gtp,f.T@p)+MASS
        return ht,htp,ut
    h0,hp0,u0=family(0.);plus=family(step);minus=family(-step)
    dh=(plus[0]-minus[0])/(2*step);dhp=(plus[1]-minus[1])/(2*step)
    du=(plus[2]-minus[2])/(2*step);connection=du@u0.conj().T
    bare=u0@dh@u0.conj().T;correction=connection@hp0-hp0@connection
    source_error=norm(dhp-bare-correction)
    missing=norm(dhp-bare)
    assert source_error<2e-8 and missing>1e-3
    return dict(original_CAR_modes=32,Nambu_coefficients=64,
        general_constant_metric=g.tolist(),nonorthogonal_map=f.tolist(),
        polar_rotation=r.tolist(),symbol_dictionary_errors=errors,
        geometric_source_step=step,full_geometric_source_error=source_error,
        omitted_spin_dictionary_derivative_error=missing,
        commutator_correction_norm=norm(correction),
        half_density_Jacobian=float(np.sqrt(np.linalg.det(f))),
        background_spin_structure_and_common_metric_are_inputs=True,
        not_a_lattice_continuum_or_Gauss_state_calculation=True)


def composition_and_metric_check():
    a=np.diag([1.,-1.,0.]);b=np.array([[0.,1.,0.],[1.,0.,0.],[0.,0.,0.]])
    def skew(x):return (x-x.T)/2
    comm=a@b-b@a
    # Vector fields xi=A x, eta=B x: [xi,eta]=(BA-AB)x.
    ka=-rotation_generator(skew(a));kb=-rotation_generator(skew(b))
    kbracket=-rotation_generator(skew(b@a-a@b))
    defect=ka@kb-kb@ka-kbracket
    predicted=-rotation_generator(comm)
    assert norm(defect-predicted)<1e-13 and abs(norm(defect)-1)<1e-13
    rows=[];g0=np.eye(3)
    for t in (.2,.1,.05):
        f=exp_sym(t*a);f2=exp_sym(t*b)
        g1,r1,u1=geometric_pullback(g0,f)
        g2,r2,u2=geometric_pullback(g1,f2)
        gt,rt,ut=geometric_pullback(g0,f@f2)
        _,rfrozen,ufrozen=geometric_pullback(g0,f2)
        composed=u2@u1
        full_error=norm(composed-ut)
        frozen_error=norm(ufrozen@u1-ut)
        jac_error=abs(np.sqrt(np.linalg.det(f2))*np.sqrt(np.linalg.det(f))-
                      np.sqrt(np.linalg.det(f@f2)))
        assert max(norm(g2-gt),norm(r1@r2-rt),full_error,jac_error)<3e-12
        assert frozen_error>.3*t*t
        rows.append(dict(parameter=t,updated_metric_error=norm(g2-gt),
            rotation_cocycle_error=norm(r1@r2-rt),spin_cocycle_error=full_error,
            frozen_metric_spin_error=frozen_error,half_density_cocycle_error=jac_error))
    return dict(linear_vector_matrices=[a.tolist(),b.tolist()],
        fixed_metric_Lie_bracket_defect_norm=norm(defect),analytic_defect_error=norm(defect-predicted),
        finite_composition_rows=rows,local_spin_sign_branch=True,
        counterexample_only_to_fixed_metric_endomorphism_composition=True,
        moving_metric_identification_does_not_supply_Einstein_dynamics=True)


def run():
    deps=('joint_chiral_source_matching.py','joint_fermion_gauss_completion.py',
        'research_note_375.md','research_note_600.md','research_note_604.md','research_note_611.md',
        'research_note_612.md','research_note_619.md','research_note_662.md',
        'round663_drafts/spinorial_local_source_probe_results.json')
    return dict(date='2026-10-02',round=663,tests_run=2,failures=0,errors=0,
        full_source_and_frame=full_source_and_frame_check(),
        composition_and_metric=composition_and_metric_check(),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original full mass and continuum Clifford symbol share local-lapse spinorial displacement, metric/frame identification and geometric source. Explicit fixed-metric composition failure is repaired by a moving-metric local spin cocycle. Given continuum geometry and spin structure; no actual lattice/chiral-state limit, regulated field-commutator anomaly calculation, ADM closure or quantum GR.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
