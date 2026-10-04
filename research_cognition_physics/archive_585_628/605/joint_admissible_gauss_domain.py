"""605: admissible configurations, Gauss thermal states and source domains.
Whole-model statements are proved in the note. Rotor and Fourier diagnostics
test domain/positivity mechanisms; they are not the full matter partition sum.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_quotient_gauge_completion as gauge
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_admissible_gauss_domain_results.json'

def roughness(C,W,z):
    return float(np.linalg.norm(np.eye(32)-matter.representation(C,W,z),2))

def quotient_domain_check():
    rng=np.random.default_rng(605)
    rows=[]
    for _ in range(12):
        C=gauge.group_exp(rng.normal(size=8)*.12,3)
        W=gauge.group_exp(rng.normal(size=3)*.12,2);z=np.exp(1j*rng.normal()*.08)
        A=gauge.group_exp(rng.normal(size=8),3)
        B=gauge.group_exp(rng.normal(size=3),2)
        u=roughness(C,W,z)
        conjugate=roughness(A@C@A.conj().T,B@W@B.conj().T,z)
        lift=roughness(np.exp(2j*np.pi/3)*C,-W,np.exp(1j*np.pi/3)*z)
        rows.append(max(abs(u-conjugate),abs(u-lift)))
    assert max(rows)<1e-12
    epsilon=.02
    inside=roughness(np.eye(3),np.eye(2),np.exp(.0005j))
    outside=roughness(np.eye(3),np.eye(2),np.exp(1j*np.pi/12))
    p=gauge.parameters()
    finite_bad_potential=float(gauge.potential(np.eye(3),np.eye(2),np.exp(1j*np.pi/12),p))
    assert inside<epsilon<outside and np.isfinite(finite_bad_potential)
    return dict(epsilon_input=epsilon,inside_norm=inside,outside_norm=outside,
        original_magnetic_potential_at_bad_point=finite_bad_potential,
        max_quotient_and_conjugation_error=max(rows),sample_count=len(rows),
        full_Gauss_support_is_analytic_not_inferred_from_sampling=True)

def rotor(N,kappa=0):
    j=np.arange(-N//2,N//2);h=2*np.pi/N;theta=h*j
    H=(2*np.eye(N)-np.roll(np.eye(N),1,axis=0)-np.roll(np.eye(N),-1,axis=0))/h**2
    H+=np.diag(kappa*(1-np.cos(theta)))
    allowed=abs(j)<N//6
    return H,allowed,h

def gibbs(H,beta):
    E,V=np.linalg.eigh(H)
    weight=np.exp(-beta*(E-E[0]));Z0=weight.sum()
    prob=weight/Z0
    rho=(V*prob)@V.conj().T
    free=float(E[0]-np.log(Z0)/beta)
    return rho,free,float(prob@E),E

def penalty_check():
    H,allowed,h=rotor(120,.4);beta=.8
    bad=~allowed;HD=H[np.ix_(allowed,allowed)]
    rhoD,FD,ED,_=gibbs(HD,beta)
    rho0,F0,_,_=gibbs(H,beta)
    rows=[];lastF=F0-1
    for strength in (0.,5.,50.,500.,5000.):
        rho,F,E,spec=gibbs(H+strength*np.diag(bad),beta)
        p=float(np.trace(rho[np.ix_(bad,bad)]))
        bound=None if strength==0 else (FD-F0)/strength
        assert 0<p<1 and F>lastF-1e-10 and F<=FD+1e-10
        if bound is not None: assert p<=bound+1e-11
        rows.append(dict(penalty=strength,bad_probability=p,free_energy=F,
            free_energy_gap_to_restricted=FD-F,variational_probability_bound=bound))
        lastF=F
    delta=.001
    derivative=(gibbs(H+(5+delta)*np.diag(bad),beta)[1]-
                gibbs(H+(5-delta)*np.diag(bad),beta)[1])/(2*delta)
    assert abs(derivative-rows[1]['bad_probability'])<1e-7
    post=rho0[np.ix_(allowed,allowed)];post/=np.trace(post)
    distance=.5*float(np.sum(abs(np.linalg.eigvalsh(post-rhoD))))
    assert distance>.05
    return dict(N=120,beta=beta,rows=rows,restricted_free_energy=FD,
        restricted_mean_energy=ED,derivative_residual=abs(derivative-rows[1]['bad_probability']),
        postselected_vs_restricted_trace_distance=distance,
        full_model_penalty_limit_proved_in_note=True)

def sharp_projection_check():
    rows=[]
    for N in (96,192,384):
        H,allowed,h=rotor(N);beta=.7
        rho,_,_,_=gibbs(H,beta)
        conditional=rho[np.ix_(allowed,allowed)]
        prob=float(np.trace(conditional));conditional/=prob
        HD=H[np.ix_(allowed,allowed)]
        _,_,energyD,_=gibbs(HD,beta)
        energyP=float(np.trace(HD@conditional))
        assert abs(prob-allowed.mean())<1e-12
        rows.append(dict(N=N,step=h,allowed_probability=prob,
            postselected_kinetic_energy=energyP,restricted_thermal_energy=energyD,
            step_times_postselected_energy=h*energyP))
    assert all(rows[i+1]['postselected_kinetic_energy']>1.7*rows[i]['postselected_kinetic_energy'] for i in range(2))
    # Interval half-width pi/3: boundary contribution tends to 1/theta0.
    assert abs(rows[-1]['step_times_postselected_energy']-3/np.pi)<.04
    assert abs(rows[-1]['restricted_thermal_energy']-rows[-2]['restricted_thermal_energy'])<.003
    return dict(rows=rows,continuum_sharp_projection_not_in_form_domain=True,
        infinite_energy_is_proved_by_H1_jump_not_numerical_extrapolation=True)

def temporal_weight_check():
    theta0=.4;kappa=1.5;ns=np.arange(65)
    def coefficients(count):
        x,w=np.polynomial.legendre.leggauss(count)
        theta=theta0*x
        weight=np.exp(kappa*(np.cos(theta)-1))
        return theta0/(2*np.pi)*(np.cos(ns[:,None]*theta)@(w*weight))
    coarse=coefficients(256);fine=coefficients(512)
    error=float(np.max(abs(coarse-fine)))
    assert error<1e-13
    negative=np.where(fine< -1e-5)[0]
    assert len(negative)>0 and fine[0]>0
    first=int(negative[0])
    return dict(theta0=theta0,kappa=kappa,quadrature_error=error,
        first_negative_mode=first,first_negative_coefficient=float(fine[first]),
        zero_mode=float(fine[0]),minimum_mode=int(np.argmin(fine)),
        minimum_coefficient=float(fine.min()),
        scope='specific temporal one-plaquette weight, not all admissible Hamiltonians')

def run():
    data=dict(quotient_domain=quotient_domain_check(),penalty=penalty_check(),
              projection=sharp_projection_check(),temporal_weight=temporal_weight_check())
    deps=('research_note_569.md','research_note_598.md','research_note_603.md','research_note_604.md',
          'joint_fermion_gauss_completion.py','joint_quotient_gauge_completion.py')
    return dict(round=605,tests_run=4,failures=0,errors=0,**data,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_Gauss_Gibbs_cannot_have_hard_admissibility=True,
            spatial_restricted_form_preserves_finite_graph_thermal_existence=True,
            penalty_and_readout_source_bounds_are_analytic=True,
            rotor_is_diagnostic_not_full_model=True,
            no_Euclidean_chiral_measure_or_SM_completion=True,
            spatial_constraint_not_full_four_dimensional_admissibility=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(dict(round=605,tests=4,all_passed=True,temporal_weight=r['temporal_weight'])))
