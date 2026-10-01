"""545: common-vacuum relations for scalar and neutral-fermion spectral weights.

Exact finite tree matrices, plus diagnostics of round 544 parameter trajectories.
Not a pole calculation, collider fit, or a no-go theorem for the running family.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_singlet_common_mass_rg as rg

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_singlet_spectrum_and_weights_results.json'


def scalar_spectrum(lh,p,ls,ratio):
    matrix=2*np.array([[lh,p*math.sqrt(ratio)],[p*math.sqrt(ratio),ls*ratio]])
    vals,vectors=np.linalg.eigh(matrix);weights=vectors[0,:]**2
    assert min(vals)>0
    schur=2*(lh-p*p/ls)
    return dict(squared_mass_over_h_squared=vals.tolist(),Higgs_weights=weights.tolist(),
        static_Schur_mass_squared_over_h_squared=schur,
        static_inverse_response=float(np.linalg.inv(matrix)[0,0]),
        low_momentum_kinetic_coefficient=1+p*p/(ls*ls*ratio))


def neutral_spectrum(c,z,ratio):
    """h=1; U Takagi convention from eigenvectors of M†M."""
    S=float(np.vdot(c,c).real);Y=np.column_stack((c,1j*c))/math.sqrt(2)
    M=np.zeros((5,5),complex);M[:3,3:]=Y/math.sqrt(2);M[3:,:3]=M[:3,3:].T
    M[3:,3:]=math.sqrt(z*ratio)*np.eye(2)
    vals,U=np.linalg.eigh(M.conj().T@M)
    heavy=U[:,-2:]@U[:,-2:].conj().T
    return vals,heavy[:3,:3],M,Y


def row_from_state(a):
    q,S,z,lh,p,ls,x,y=map(float,a);_,A,Q=rg.invariants(a);ratio=Q/A
    assert ratio>0
    scalar=scalar_spectrum(lh,p,ls,ratio)
    mn2=z*ratio+S/2;alpha=(S/2)/mn2
    return dict(q=q,S=S,ratio=float(ratio),scalar=scalar,
        neutral_heavy_squared_over_h_squared=mn2,active_weight_in_heavy_pair=alpha,
        top_squared_over_h_squared=q/2,neutral_over_top_squared_ratio=mn2/(q/2),
        maximum_scalar_over_four_top_squared_ratio=max(scalar['squared_mass_over_h_squared'])/(2*q))


def run():
    checks=[];rng=np.random.default_rng(545);error=0.;projector_error=0.
    for _ in range(18):
        c=rng.normal(size=3)+1j*rng.normal(size=3);c*=math.sqrt(.37)/np.linalg.norm(c)
        z=.4;ratio=.6;S=float(np.vdot(c,c).real);mn2=z*ratio+S/2
        vals,P,M,Y=neutral_spectrum(c,z,ratio)
        alpha=S/(2*z*ratio+S)
        expected=alpha*np.outer(c.conj(),c)/S
        error=max(error,float(max(abs(vals-np.array([0,0,0,mn2,mn2])))))
        projector_error=max(projector_error,float(np.max(abs(P-expected))))
        assert np.allclose(P,expected,atol=2e-14)
        assert abs(np.trace(P).real-alpha)<2e-14
        assert np.linalg.norm(Y@Y.T)<1e-14
    assert error<1e-14 and projector_error<1e-14
    checks.append('exact_Takagi_spectrum_and_two_heavy_state_active_projector_for_complex_flavours')

    rows=[]
    for r in (F(2),F(7,3),F(4)):
        for alpha in (F(1,100),F(1,5),F(4,5)):
            q=F(3,10);S=alpha*q;T=3*q+r*S;z=T/(2*r);ratio=r*(q-S)/T
            lh=(3*q*q+r*S*S)/T;p=S;ls=T/r
            beta=3*(q-S)/T
            assert z*ratio+S/2==q/2
            assert (S/2)/(z*ratio+S/2)==alpha
            assert beta==3*(1-alpha)/(3+r*alpha)
            zh=r*(1-alpha)/(r+3)
            assert zh==ratio/(1+ratio)
            spec=scalar_spectrum(*map(float,(lh,p,ls,ratio)))
            assert np.allclose(spec['squared_mass_over_h_squared'],[float(2*q*beta),float(2*q)],atol=2e-15)
            assert np.allclose(spec['Higgs_weights'],[float(zh),float(1-zh)],atol=2e-14)
            rows.append(dict(r=str(r),alpha=str(alpha),light_scalar_over_four_top_squared=str(beta),
                             light_Higgs_weight=str(zh)))
    checks.append('common_bare_boundary_exact_scalar_neutral_top_and_weight_relations')

    # Resource/tolerance parameters here are diagnostic contracts, not empirical bounds.
    r=F(7,3);epsilon=F(1,100);beta_min=3*(1-epsilon)/(3+r*epsilon)
    for alpha in (F(1,1000),F(1,500),epsilon):
        beta=3*(1-alpha)/(3+r*alpha)
        assert beta>=beta_min
        assert alpha==3*(1-beta)/(3+r*beta)
        assert r*(1-alpha)/(r+3)<F(1,2)
    for r2 in (F(2),F(3),F(4)):
        a=F(1,100)
        assert (r2*(1-a)/(r2+3)>F(1,2))==(a<(r2-3)/(2*r2))
    checks.append('exact_active_weight_light_mass_tradeoff_and_Higgs_dominance_criterion')

    # Generic spectral sum rules distinguish static inverse response from an eigenmass.
    sum_error=0.;inv_error=0.
    for _ in range(20):
        lh=.2+rng.random();ls=.3+rng.random();p=.5*math.sqrt(lh*ls)*rng.random();ratio=.05+2*rng.random()
        spec=scalar_spectrum(lh,p,ls,ratio);e=np.array(spec['squared_mass_over_h_squared']);w=np.array(spec['Higgs_weights'])
        schur=spec['static_Schur_mass_squared_over_h_squared']
        sum_error=max(sum_error,abs(float(w@e)-2*lh))
        inv_error=max(inv_error,abs(float(w@(1/e))-1/schur))
        assert np.isclose(w@e,2*lh,rtol=2e-14)
        assert np.isclose(w@(1/e),1/schur,rtol=2e-14)
        assert e[0]<=schur<=2*lh+1e-14
    checks.append('full_Higgs_spectral_sum_and_inverse_sum_rules')

    # Exact inverse-propagator expansion; small momentum is needed as well as stable V.
    a,b,c=F(3),F(1),F(5);P=F(1,100)
    exact=P-a-b*b/(P-c);leading=(1+b*b/(c*c))*P-(a-b*b/c)
    remainder=b*b*P*P/(c*c*(c-P))
    assert exact-leading==remainder
    assert 1+b*b/(c*c)>1
    checks.append('momentum_dependent_Schur_complement_and_nontrivial_kinetic_normalization')

    # Both endpoints have different qualifications; do not assign a unique degenerate residue.
    q=.3;r=rg.R;S=0.;T=3*q
    degenerate=2*q*np.eye(2)
    assert np.linalg.eigvalsh(degenerate).tolist()==[2*q,2*q]
    for theta in (.2,.7):
        U=np.array([[math.cos(theta),-math.sin(theta)],[math.sin(theta),math.cos(theta)]])
        assert np.allclose(U.T@degenerate@U,degenerate)
    cvec=np.array([math.sqrt(q),0,0],complex)
    vals,P,_,_=neutral_spectrum(cvec,T/(2*r),0.)
    assert np.allclose(vals,[0,0,0,q/2,q/2],atol=1e-14)
    assert abs(np.trace(P).real-1)<1e-14
    checks.append('zero_portal_degenerate_scalar_and_zero_Majorana_block_endpoints')

    saved=json.loads(rg.TARGET.read_text('utf8'));running=[]
    for original in saved['examples']:
        if original['u']<18:continue
        a=[original['state'][k] for k in ('q','S','z','lambda_H','p','lambda_s','x','y')]
        row=row_from_state(a);row.update(q0=original['q0'],u=original['u'])
        S=row['S'];cvec=np.array([math.sqrt(S),0,0],complex)
        vals,proj,_,_=neutral_spectrum(cvec,a[2],row['ratio'])
        assert np.allclose(vals[-2:],row['neutral_heavy_squared_over_h_squared'],atol=1e-14)
        assert abs(np.trace(proj).real-row['active_weight_in_heavy_pair'])<1e-14
        running.append(row)
    balanced=running[0]
    assert balanced['scalar']['Higgs_weights'][0]<.5
    assert balanced['active_weight_in_heavy_pair']>.4
    assert abs(balanced['neutral_over_top_squared_ratio']-1)>.1
    assert abs(balanced['scalar']['static_Schur_mass_squared_over_h_squared']-
               balanced['scalar']['squared_mass_over_h_squared'][0])>.1
    checks.append('same_frozen_RG_trajectories_change_bare_relations_but_retain_joint_observable_conditions')

    deps=('research_note_543.md','research_note_544.md','joint_singlet_common_mass_rg.py',
          'joint_singlet_common_mass_rg_results.json','research_round_544_checks.json')
    return dict(round=545,tests_run=len(checks),failures=0,errors=0,checks=checks,
        neutral_spectrum_max_error=error,active_projector_max_error=projector_error,
        scalar_sum_max_error=sum_error,scalar_inverse_sum_max_error=inv_error,
        rational_common_boundary_cases=rows,
        diagnostic_tolerance_not_experimental=dict(r='7/3',epsilon='1/100',
            light_scalar_squared_over_four_top_squared_lower=str(beta_min)),
        running_spectra=running,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(common_bare_relations_require_same_scale_and_0_less_S_less_q=True,
          exact_tree_spectral_weights_not_pole_observables=True,
          frozen_running_trajectories_not_rescanned=True,physical_decoupling_completed=False,
          precision_data_fit=False,universal_running_family_excluded=False,unification_completed=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','neutral_spectrum_max_error','active_projector_max_error')},ensure_ascii=False))
