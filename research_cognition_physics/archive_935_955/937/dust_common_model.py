"""937: positive material reference in the existing Einstein--gauge background.

Classical two-derivative physics. The lattice is a collocation check, not a
fundamental space. The moving-point test is not a global constrained solution.
"""
from pathlib import Path
from fractions import Fraction as Q
import sys, json, hashlib, argparse
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STAGE = HERE.parent
TARGET = HERE/'dust_common_model_results.json'
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import Layout, ResearchRuntime


def maximum(a):
    return float(np.max(np.abs(a)))


def fraction_check():
    # Future dust: C=-5, Ca=(3,0,0), q=I, root=4. All statements exact.
    c, a, r, volume = Q(-5), Q(3), Q(4), Q(1)
    lapse, shift = -c/r, a/r
    rho = r*r/(volume*(-c))
    un, ua = c/r, a/r
    assert -un*un+ua*ua == -1
    assert -lapse*lapse+shift*shift == -1
    assert rho*un*un == -c/volume
    assert rho*un*ua == -a/volume
    assert c+r*(-un) == 0 and a-r*ua == 0
    # Algebra of the two distinct even, square-zero Grassmann bilinears.
    def mul(x, y):
        z = [Q(0)]*4
        for i in range(4):
            for j in range(4):
                if not i & j:
                    z[i | j] += x[i]*y[j]
        return z
    u, v = Q(3,7), Q(5,11)
    cn = [c,u,v,Q(0)]
    rad = mul(cn,cn);rad[0] -= a*a
    z = rad.copy();z[0] -= r*r
    z2 = mul(z,z)
    root = [r*(i == 0)+z[i]/(2*r)-z2[i]/(8*r**3) for i in range(4)]
    assert mul(root,root) == rad
    mixed = -root[3]
    assert mixed == a*a/r**3*u*v and mixed > 0
    return dict(C=str(c), Ca=[str(a),'0','0'], root=str(r),
        lapse=str(lapse), shift=[str(shift),'0','0'], proper_density=str(rho),
        material_H_is_negative_root=True, positive_normal_energy=str(-c),
        exact_total_constraints_zero=True, exact_unit_timelike=True,
        exact_mixed_mass_four_leg=str(mixed),
        mass_bilinear_amplitudes=[str(u),str(v)],
        scope='Local constraint/source algebra; not a global moving-dust PDE solution.')


def moving_source_check():
    # Complete variation includes explicit inverse-metric dependence of D.
    qi = np.array([[1.3,.1,.03],[.1,.9,-.02],[.03,-.02,1.1]])
    a = np.array([.37,-.24,.18]);c = -1.7
    dc = .11;da = np.array([.05,-.03,.02])
    dq = np.array([[.04,.02,0],[.02,-.03,.01],[0,.01,.05]])
    def hh(cc,aa,gg):return -np.sqrt(cc*cc-aa@gg@aa)
    d = float(a@qi@a);r = np.sqrt(c*c-d)
    predicted = -c/r*dc+(qi@a)@da/r+(a@dq@a)/(2*r)
    metric_term = float((a@dq@a)/(2*r))
    fd=[]
    for e in (.01,.005):
        value=(hh(c+e*dc,a+e*da,qi+e*dq)-hh(c-e*dc,a-e*da,qi-e*dq))/(2*e)
        fd.append(dict(step=e,value=float(value),error=abs(float(value-predicted))))
    assert 3.9 < fd[0]['error']/fd[1]['error'] < 4.1
    assert abs(metric_term)>1e-5
    volume=1/np.sqrt(np.linalg.det(qi));density=r*r/(volume*(-c))
    lapse=-c/r;shift=qi@a/r;un=c/r;ua=a/r
    g00=-lapse*lapse+shift@np.linalg.solve(qi,shift)
    # Compare physical stress, solved dust constraints and lapse/shift directly.
    assert abs(g00+1)<1e-14
    assert abs(density*un*un+c/volume)<1e-14
    assert maximum(density*un*ua+a/volume)<1e-14
    assert abs(c+np.sqrt(r*r+d))<1e-14
    assert maximum(a-r*ua)<1e-14
    u,v=.23,.31
    mixed=d/r**3*u*v
    source=[]
    for e in (.01,.005):
        def h(s,t):return hh(c+s*u+t*v,a,qi)
        val=(h(e,e)-h(e,-e)-h(-e,e)+h(-e,-e))/(4*e*e)
        source.append(dict(step=e,value=float(val),error=abs(float(val-mixed))))
    assert 3.7<source[0]['error']/source[1]['error']<4.3
    return dict(non_dust_C=c, inverse_metric=qi.tolist(), non_dust_Ca=a.tolist(),
        density=density, root=r, predicted_full_variation=float(predicted),
        explicit_metric_variation=metric_term, variation_checks=fd,
        mixed_mass_coefficient=float(mixed), mixed_mass_checks=source,
        material_g00=float(g00), scope='Nonzero momentum local algebra, not a solved global background.')


def background_checks(old):
    geo=old.geo;rows=[]
    for n in (16,24):
        z,psi,tensor,info,color=old.completed(n,1.)
        tau=np.sqrt(z['tau2']);vol=psi**6;metric=psi[...,None,None]**4*np.eye(3)
        inverse=psi[...,None,None]**-4*np.eye(3)
        curvature=-8*psi**-5*geo.laplace(psi)
        shear=psi[...,None,None]**-2*tensor
        rho0=.5*psi**-12*z['pKp']+.5*psi**-4*z['B']+z['U']+psi**-8*z['Y']
        matter_C=vol*rho0
        def canonical(kk,tt):
            # K=-L_n q/2, as required by the inherited div(A)=-M convention.
            kup=psi[...,None,None]**-8*kk
            return -.5*vol[...,None,None]*(kup-tt*inverse)
        def gravity_C(pp):
            trace=np.einsum('...ij,...ij->...',metric,pp)
            p2=psi**8*np.sum(pp*pp,axis=(-1,-2))
            return 2/vol*(p2-.5*trace*trace)-.5*vol*curvature
        k0=shear+tau/3*metric;p0=canonical(k0,tau)
        c0=gravity_C(p0)+matter_C
        momentum=sum(geo.derivative(tensor[..., :,j],j) for j in range(3))+z['mom']
        k=old.source.tangent_matrix(z)
        gauss=old.source.gauss(z,k,z['p'],z['f']['E'],z['f']['E0'])
        gauss_error=maximum(gauss)
        assert gauss_error<1e-9 and maximum(color['G'])<1e-14
        assert color['electric']>0 and color['magnetic']>0
        for delta in (.01,.03):
            tp=tau+delta;rhoD=(tp*tp-tau*tau)/3
            knew=shear+tp/3*metric;pnew=canonical(knew,tp)
            cnew=gravity_C(pnew)+matter_C
            pdust=vol*rhoD
            shift_error=maximum(pnew-p0-vol[...,None,None]*delta/3*inverse)
            identity=maximum(cnew-c0+pdust)
            total=maximum(cnew+pdust)
            k2=psi**-8*np.sum(knew*knew,axis=(-1,-2))
            einstein=curvature-k2+tp*tp-2*(rho0+rhoD)
            assert shift_error<1e-13 and identity<1e-13
            assert total<3e-8 and maximum(einstein)<3e-8
            assert maximum(momentum)<1e-9
            assert np.max(cnew)<0 and np.min(pdust)>0
            # Body at rest: H=C. Source Hessian in Ca is q^ab / r, not zero.
            r=-cnew
            sensitivity=maximum(inverse/r[...,None,None])
            rows.append(dict(N=n,trace_shift=delta,old_tau=float(tau),new_tau=float(tp),
                dust_density=float(rhoD),dust_energy_integral=float(z['dx']**3*np.sum(pdust)),
                old_constraint_density_residual=maximum(c0),new_total_density_residual=total,
                new_Einstein_constraint_residual=maximum(einstein),
                canonical_trace_shift_error=shift_error,exact_shift_identity_residual=identity,
                unchanged_momentum_constraint_residual=maximum(momentum),
                unchanged_electroweak_Gauss_residual=gauss_error,color_Gauss_residual=maximum(color['G']),
                color_electric_coefficient=color['electric'],color_magnetic_coefficient=color['magnetic'],
                minimum_reference_momentum=float(np.min(pdust)),
                maximum_rest_momentum_source_hessian=sensitivity,
                at_rest_mixed_mass_four_leg_coefficient=0.,
                space_and_matter_canonical_data_preserved=True))
    assert rows[0]['maximum_rest_momentum_source_hessian']>2.9*rows[1]['maximum_rest_momentum_source_hessian']
    return rows


def run():
    exact=fraction_check();moving=moving_source_check()
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        rows=background_checks(old)
    oldpaths=[STAGE.parent/'archive_742_763/753/joint_reference_constraint_strata.py',
              STAGE.parent/'archive_554_584/572/joint_gauss_einstein_initial_data.py',
              STAGE/'research_note_873.md',STAGE/'936/native_source_quantization_results.json']
    paths=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/dust_candidate_decision.md']+oldpaths
    return dict(round=937,date='2026-10-07',all_scientific_checks_passed=True,
        exact_reference_and_source_check=exact,moving_local_source_check=moving,
        original_full_field_backgrounds=rows,
        new_positive_reference_matter_is_explicit_physical_input=True,
        classical_two_derivative_and_zero_fermion_background_scope=True,
        full_parent_action_retains_fermions_but_not_fully_quantized=True,
        original873_clock_is_different_physical_reference=True,
        dust_low_density_limit_not_claimed_uniform=True,
        classical_constraints_not_quantum_Ward_certificate=True,
        no_new_spatial_dimension_claim=True,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
