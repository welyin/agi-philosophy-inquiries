"""703: one fixed confining anchor for geometric thermal-source limits.

Full finite-graph/Gauss theorem is analytic. Group1 calibrates original589
coefficient inequalities, including non-diagonal shape. Group2 retains702's
declared radial/neutral-Majorana sector with actual volume-dependent kinetic
and potential, and evaluates source jets without finite-difference fitting.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_full_spatial_metric as geometry
import joint_gibbs_preserving_transfer as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_geometry_thermal_limit_results.json'


def muljet(x,y):
    return [x[0]@y[0],x[1]@y[0]+x[0]@y[1],
            x[2]@y[0]+2*x[1]@y[1]+x[0]@y[2]]


def expjet(h,g,c,t):
    """Taylor scaling/squaring of value/first/second jets, complex h allowed."""
    scale=max(0,int(np.ceil(np.log2(max(1.,8*t*np.linalg.norm(h,np.inf))))))
    x=[-t*m/(2**scale) for m in (h,g,c)]
    term=[np.eye(len(h),dtype=complex),np.zeros_like(h,dtype=complex),np.zeros_like(h,dtype=complex)]
    total=[m.copy() for m in term]
    for k in range(1,25):
        term=[m/k for m in muljet(term,x)]
        total=[a+b for a,b in zip(total,term)]
    assert max(np.linalg.norm(m,np.inf) for m in term)<1e-17
    for _ in range(scale):total=muljet(total,total)
    return total


def powerjet(s,n):
    out=[np.eye(len(s[0]),dtype=complex),np.zeros_like(s[0]),np.zeros_like(s[0])]
    while n:
        if n&1:out=muljet(out,s)
        n//=2
        if n:s=muljet(s,s)
    return out


def statejets(p):
    z=np.array([np.trace(m).real for m in p]);rho=p[0]/z[0]
    first=p[1]/z[0]-rho*z[1]/z[0]
    second=p[2]/z[0]-rho*z[2]/z[0]-2*first*z[1]/z[0]
    logjets=[np.log(z[0]),z[1]/z[0],z[2]/z[0]-(z[1]/z[0])**2]
    return [rho,first,second],logjets


def original_geometry_anchor():
    rv,rs,rn=.12,.35,.18
    lower=float(np.exp(-rn-max(6*rv,2*rv+2*rs)));eta=lower/4
    rng=np.random.default_rng(703);rows=[]
    for k in range(18):
        raw=rng.normal(size=(3,3));raw=(raw+raw.T)/2;raw-=np.trace(raw)*np.eye(3)/3
        raw*=rs/np.linalg.norm(raw,2)
        shape=geometry.shape_exp(raw)
        # Original edge/face means of positive vertex scales; all in declared box.
        u=rng.uniform(-rv,rv,size=9);psi=np.exp(u)
        pe=(psi[0]+psi[1:4])/2
        pf=np.array([np.mean(psi[[0,1,2,4]]),np.mean(psi[[0,1,3,5]]),np.mean(psi[[0,2,3,6]])])
        lapse=float(np.exp(rng.uniform(-rn,rn)))
        ce=lapse*shape/(pe[:,None]*pe[None,:])
        cg=lapse*np.linalg.inv(shape)*(pe[:,None]*pe[None,:])
        cb=lapse*geometry.gram(shape)/(pf[:,None]*pf[None,:])
        candidates=[lapse*np.exp(-6*u[0]),lapse*np.exp(6*u[0]),
            np.linalg.eigvalsh(ce)[0],np.linalg.eigvalsh(cg)[0],np.linalg.eigvalsh(cb)[0]]
        actual=float(min(candidates));assert actual>=lower-1e-13
        assert np.linalg.eigvalsh(ce-eta*np.eye(3))[0]>0
        rows.append(dict(case=k,minimum_original_coefficient=actual,
            anchor_remainder_margin=actual-eta,
            electric_cross_size=float(np.linalg.norm(ce-np.diag(np.diag(ce))))))
    return dict(log_volume_radius=rv,log_shape_norm_radius=rs,log_lapse_radius=rn,
        analytic_uniform_coefficient_lower=lower,anchor_fraction=eta,rows=rows,
        bound_is_analytic_over_whole_box_not_inferred_from_samples=True,
        original_non_diagonal_electric_gradient_and_magnetic_coefficients=True,
        not_a_quantum_Gauss_thermal_simulation=True)


def volume_thermal_sources():
    K,W,B,effects=old.interacting_radial_diagnostic()
    eta=.15;A=eta*(K+W);beta=1.2;u0=.03
    def family(u):return np.exp(-6*u)*K+np.exp(6*u)*W+B
    H=family(u0);G=-6*np.exp(-6*u0)*K+6*np.exp(6*u0)*W
    C=36*np.exp(-6*u0)*K+36*np.exp(6*u0)*W
    target=expjet(H,G,C,beta);rhot,logt=statejets(target)
    exact=old.exp_h(H,beta);assert np.linalg.norm(target[0]-exact)<2e-12
    e,v=np.linalg.eigh(H);weight=np.exp(-beta*e);Z=float(weight.sum())
    gt=v.conj().T@G@v;ct=v.conj().T@C@v
    z=beta*(e[None,:]-e[:,None])/2
    sinhc=np.divide(np.sinh(z),z,out=np.ones_like(z),where=abs(z)>1e-10)
    kernel=np.exp(-beta*(e[:,None]+e[None,:])/2)*sinhc
    first=-beta*float(weight@np.diag(gt).real)/Z
    contact=-beta*float(weight@np.diag(ct).real)/Z
    no_contact=beta**2*float(np.sum(kernel*abs(gt)**2))/Z-first**2
    exactsecond=no_contact+contact
    assert abs(first-logt[1])<2e-11 and abs(exactsecond-logt[2])<2e-10
    O=effects[0]@effects[0]
    target_p=[float(np.trace(O@r).real) for r in rhot]
    rows=[]
    for n in (4,8,16,32,64):
        a=beta/n;E=old.exp_h(A,a/2)
        middle=expjet(H-A,G,C,a)
        step=[E@m@E for m in middle]
        vals=np.linalg.eigvalsh(step[0]);assert min(vals)>0
        p=powerjet(step,n);rho,logs=statejets(p)
        errors=[old.norm1(rho[j]-rhot[j]) for j in range(3)]
        record_errors=[abs(float(np.trace(O@rho[j]).real)-target_p[j]) for j in range(3)]
        rows.append(dict(slices=n,positive_step_min=float(min(vals)),
            Gibbs_jet_trace_errors=errors,log_partition_jet_errors=[abs(logs[j]-logt[j]) for j in range(3)],
            original_record_probability_jet_errors=record_errors,
            Hessian=float(logs[2])))
    assert rows[-1]['Gibbs_jet_trace_errors'][2]<rows[0]['Gibbs_jet_trace_errors'][2]/100
    assert rows[-1]['log_partition_jet_errors'][2]<rows[0]['log_partition_jet_errors'][2]/100
    # Direct thermal finite differences independently check the source convention.
    fd=[]
    for h in (.0004,.0002):
        zp=np.log(np.trace(old.exp_h(family(u0+h),beta)).real)
        zm=np.log(np.trace(old.exp_h(family(u0-h),beta)).real)
        fd.append(dict(step=h,first_error=abs((zp-zm)/(2*h)-first),
            second_error=abs((zp-2*logt[0]+zm)/h**2-exactsecond)))
    assert fd[-1]['second_error']<2e-5
    # Uniform complex coercivity for this original conformal diagnostic.
    rreal,rimag=.1,.06;coercivity=float(np.exp(-6*rreal)*np.cos(6*rimag)-eta)
    assert coercivity>0
    c=max(0.,-float(np.linalg.eigvalsh(B)[0]));ZA=float(np.exp(-beta*np.linalg.eigvalsh(A)).sum())
    complex_rows=[]
    for u in (.025+.04j,-.06+.03j):
        n=8;a=beta/n;E=old.exp_h(A,a/2);h=family(u)-A
        middle=expjet(h,np.zeros_like(h),np.zeros_like(h),a)[0]
        p=np.linalg.matrix_power(E@middle@E,n)
        norm=float(np.linalg.svd(p,compute_uv=False).sum());bound=float(np.exp(beta*c)*ZA)
        assert norm<bound
        complex_rows.append(dict(real=u.real,imag=u.imag,trace_norm=norm,uniform_bound=bound))
    return dict(dimension=len(H),eta=eta,beta=beta,log_scale=u0,
        original_volume_powers=dict(kinetic=-6,onsite=6,canonical_local_mass=0),
        target_log_partition_jets=logt,target_record_probability_jets=target_p,
        source_Hessian_contact=contact,incorrect_Hessian_if_contact_omitted=no_contact,
        independent_Duhamel_Hessian_error=abs(exactsecond-logt[2]),finite_difference_check=fd,
        rows=rows,complex_coercivity_margin=coercivity,complex_checks=complex_rows,
        finite_radial_neutral_diagnostic_scope_inherited_702=True,
        full_graph_general_geometry_theorem_not_claimed_numerically_proved=True)


def run():
    deps=('research_note_589.md','research_note_590.md','research_note_603.md',
        'research_note_623.md','research_note_624.md','research_note_655.md','research_note_662.md',
        'research_note_665.md','research_note_702.md','joint_full_spatial_metric.py',
        'joint_gibbs_preserving_transfer.py','round703_drafts/bounded_source_entry.md')
    return dict(date='2026-10-02',round=703,tests_run=2,failures=0,errors=0,
        original_geometry_anchor=original_geometry_anchor(),thermal_geometry_jets=volume_thermal_sources(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_full_fixed_graph_analytic_theorem=True,
            fixed_anchor_is_new_finite_step_representation=True,
            actual_geometry_source_and_contact_thermal_limits=True,
            source_derivatives_not_inferred_from_trace_distance_alone=True,
            all_interactions_and_original_Gauss_preserved=True,
            real_time_source_derivatives_not_proved_for_the_approximation=True,
            quantum_geometry_continuum_and_GR_not_derived=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=703,groups=2,last=result['thermal_geometry_jets']['rows'][-1],all_checks_passed=True)))
