"""714: original charged cut operator, common sources and thermal constraints.

Numerics evaluate original coefficient functions, not the full interacting
Gibbs state. Thermal trace statements are proved in research_note_714.md.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_region_energy_gluing as old
import joint_full_spatial_metric as metric
import joint_holonomy_readout_scale as loops

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_center_thermal_reference_results.json'
CENTER=(np.eye(3),-np.eye(2),1.)


def act(g,p):
    x=g[2]**3*g[1]@(p[:2]+1j*p[2:4])
    return np.r_[x.real,x.imag,p[4]]


def vectors(phi,links):
    return np.stack([metric.target_log(phi[0],act(g,q))
                     for g,q in zip(links,phi[1:])])


def gradient_difference(phi,links,shape,sigma=.12,eps=.73):
    changed=links.copy();changed[0]=old.product(CENTER,changed[0])
    u=np.exp(sigma)*vectors(phi,links)
    v=np.exp(sigma)*vectors(phi,changed)
    K=metric.original.metric(phi[0]);inverse=np.linalg.inv(shape)
    gram0=u@K@u.T;gram1=v@K@v.T
    direct=float(eps/2*np.sum(inverse*(gram1-gram0)))
    d=v[0]-u[0]
    expanded=float(eps/2*(2*np.sum(inverse[0]*(u@K@d))+inverse[0,0]*(d@K@d)))
    return direct,expanded,gram1-gram0


def charged_cut_check():
    rng=np.random.default_rng(7143)
    Rcenter=old.matter.representation(*CENTER)
    P=(np.eye(32)-Rcenter)/2
    assert np.linalg.matrix_rank(P)==16
    assert np.linalg.norm(P@P-P)<1e-13
    scalar_errors=[];fermion_errors=[];mass_errors=[]
    for _ in range(18):
        phi=rng.normal(size=(4,5))*.35
        links=[old.sample(rng) for _ in range(3)]
        A=rng.normal(size=(3,3));shape=np.eye(3)+A@A.T
        shape/=np.linalg.det(shape)**(1/3)
        direct,expanded,_=gradient_difference(phi,links,shape)
        scalar_errors.append(abs(direct-expanded))
        h1,d1=old.matter.mass_matrices(phi[0]);h2,d2=old.matter.mass_matrices(phi[1])
        h=np.zeros((64,64),complex);h[:32,:32]=h1;h[32:,32:]=h2
        hc=h.copy()
        t=.31+.17j
        R=old.matter.representation(*links[0])
        RC=old.matter.representation(*old.product(CENTER,links[0]))
        h[:32,32:]=t*R;h[32:,:32]=(t*R).conj().T
        hc[:32,32:]=t*RC;hc[32:,:32]=(t*RC).conj().T
        predicted=np.zeros_like(h)
        predicted[:32,32:]=-2*t*P@R
        predicted[32:,:32]=predicted[:32,32:].conj().T
        fermion_errors.append(float(np.linalg.norm(hc-h-predicted)))
        mass_errors.append(float(np.linalg.norm(hc[:32,:32]-h[:32,:32])+
                                 np.linalg.norm(hc[32:,32:]-h[32:,32:])))
        assert np.linalg.norm(d1)+np.linalg.norm(d2)>0
    assert max(scalar_errors)<1e-12 and max(fermion_errors+mass_errors)<1e-12
    return dict(samples=18,negative_center_modes=16,total_original_modes_per_node=32,
                scalar_full_cross_term_error=max(scalar_errors),
                fermion_two_node_coefficient_error=max(fermion_errors),
                all_onsite_mass_and_pairing_coefficients_retained=True,
                onsite_mass_change=max(mass_errors))


def source_check():
    rng=np.random.default_rng(7144)
    phi=rng.normal(size=(4,5))*.31;links=[old.sample(rng) for _ in range(3)]
    S=np.array([[.3,.2,-.1],[.2,-.2,.07],[-.1,.07,-.1]])
    sigma=.12;t=.27;eps=.73;step=2e-5
    shape=metric.shape_exp(S,t);inverse=np.linalg.inv(shape)
    delta,_,gram=gradient_difference(phi,links,shape,sigma,eps)
    ds=(gradient_difference(phi,links,shape,sigma+step,eps)[0]-
        gradient_difference(phi,links,shape,sigma-step,eps)[0])/(2*step)
    dt=(gradient_difference(phi,links,metric.shape_exp(S,t+step),sigma,eps)[0]-
        gradient_difference(phi,links,metric.shape_exp(S,t-step),sigma,eps)[0])/(2*step)
    predicted=float(-eps/2*np.sum((inverse@S)*gram))
    assert abs(ds-2*delta)<1e-8 and abs(dt-predicted)<1e-8
    return dict(original_cut_gradient_difference=delta,conformal_source=2*delta,
                conformal_finite_difference_error=abs(ds-2*delta),
                shape_source=predicted,shape_finite_difference_error=abs(dt-predicted),
                thermal_reference_response_not_included_in_this_coefficient_check=True)


def physical_witness_check():
    p=np.array([.4,-.3,.2,.1,.35])
    phi=np.repeat(p[None,:],4,axis=0)
    links=[loops.IDENTITY]*3
    S=np.array([[.3,.2,-.1],[.2,-.2,.07],[-.1,.07,-.1]])
    shape=metric.shape_exp(S,.27)
    delta,expanded,_=gradient_difference(phi,links,shape)
    assert delta>0 and abs(delta-expanded)<1e-13
    h,d=old.matter.mass_matrices(p)
    assert np.linalg.norm(d)>0
    rows=[]
    for m in (0.,.1,.36752608040636575,.9):
        eta=.6;beta=2.;q=eta*m
        prob=np.array([(1+q)/2,(1-q)/2])
        reversed_prob=prob[::-1]
        direct=float(np.sum(reversed_prob*np.log(reversed_prob/prob)))
        bound=float(2*q*np.arctanh(q))
        assert abs(direct-bound)<1e-13 and bound>=2*q*q-1e-14
        rows.append(dict(hypothetical_thermal_mean=m,record_TV=q,
                         beta_times_work_lower_bound=bound,work_lower_bound_at_beta_two=bound/beta))
    return dict(original_uniform_scalar_configuration=p.tolist(),
                strict_cut_gradient_witness_per_cut_node=delta,
                physical_normal_Gauss_witness_obtained_analytically_by_orbit_localization=True,
                point_configuration_not_claimed_as_a_normal_state=True,
                binary_bound_calibration=rows,
                table_means_not_computed_from_original_Gibbs=True)


def run():
    names=('research_note_589.md','research_note_598.md','research_note_603.md',
           'research_note_623.md','research_note_638.md','research_note_703.md',
           'research_note_713.md','round714_drafts/center_reference_entry.md')
    return dict(round=714,tests_run=3,failures=0,errors=0,
                charged_cut=charged_cut_check(),sources=source_check(),witness=physical_witness_check(),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original full finite graph, positive fixed geometry, exact Gibbs trace theorems analytic. '
                      'Numerics check original charged coefficients and sources; no interacting thermal mean, '
                      'autonomous center operation, uniform spatial continuum cost, or GR derivation computed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
