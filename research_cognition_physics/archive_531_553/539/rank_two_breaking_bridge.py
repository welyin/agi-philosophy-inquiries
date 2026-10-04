"""Round 539: rank-two breaking, common moments, full matrix RG and dim-5 EFT.

One-loop MS dimension-four running plus tree common-scale matching.
Charged Yukawas other than top are zero. No finite one-loop matching or pole fit.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import protected_pair_rg as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'rank_two_breaking_bridge_results.json'
LOOP=16*math.pi**2
VEV=246.0  # Declared tree conversion, not a computed running vev.
TARGET_MASSES=np.sqrt([2.513e-3,7.49e-5])  # eV; NO, m1=0, NuFIT 6.0 IC24.
HIST,BOUND=old.inputs()
XSTAR=HIST['matched_inverse_couplings']
MU=BOUND['matching_scale_GeV']
END=math.log(HIST['MZ_GeV']/MU)


def construct(a,b,M,S,U=None):
    """C=a*u*u^T+b*v*v^T, exact total Yukawa norm S in a Dirac basis."""
    assert a>=b>0 and S>=M*(a+b)
    if U is None: U=np.eye(3,dtype=complex)
    q=(S+math.sqrt(max(0,S*S-(M*(a+b))**2)))/2
    alpha=math.sqrt(q/(a+b))
    c=alpha*(math.sqrt(a)*U[:,0]+1j*math.sqrt(b)*U[:,1])
    d=M/(2*alpha)*(math.sqrt(a)*U[:,0]-1j*math.sqrt(b)*U[:,1])
    return np.column_stack([c,d])


def sharp_min(a,b,M,r,T):
    B=M*(a+b); w=T/(3+r)
    if B>T/r: return None
    S=max(B,w); x=(T-r*S)/3
    return (3*x*x+r*(S*S-2*M*M*a*b))/T


def rk4(fun,state,t0,t1,steps):
    z=np.array(state,complex); dt=(t1-t0)/steps
    for j in range(steps):
        t=t0+j*dt
        k1=fun(t,z);k2=fun(t+dt/2,z+dt*k1/2)
        k3=fun(t+dt/2,z+dt*k2/2);k4=fun(t+dt,z+dt*k3)
        z+=dt*(k1+2*k2+2*k3+k4)/6
    assert np.all(np.isfinite(z))
    return z


def full_rhs(t,z):
    yt,lam=z[:2].real;Y=z[2:8].reshape(3,2);B=z[8:12].reshape(2,2)
    H=Y.conj().T@Y; S=np.trace(H).real;Q=np.trace(H@H).real
    gy,gw,gc=old.gauge_squared(t,XSTAR)
    dy=yt*(4.5*yt*yt+S-17*gy/12-9*gw/4-8*gc)
    dl=24*lam*lam-(3*gy+9*gw)*lam+3/8*(gy*gy+2*gy*gw+3*gw*gw)
    dl+=4*(3*yt*yt+S)*lam-2*(3*yt**4+Q)
    dY=1.5*Y@H+(3*yt*yt+S-.75*gy-2.25*gw)*Y
    dB=H.T@B+B@H
    return np.r_[dy,dl,dY.ravel(),dB.ravel()]/LOOP


def eft_rhs(t,z,quartic_factor=4):
    yt,lam=z[:2].real; gy,gw,gc=old.gauge_squared(t,XSTAR)
    dy=yt*(4.5*yt*yt-17*gy/12-9*gw/4-8*gc)
    dl=24*lam*lam-(3*gy+9*gw)*lam+3/8*(gy*gy+2*gy*gw+3*gw*gw)
    dl+=12*yt*yt*lam-6*yt**4
    # z[2] integrates ln of the multiplicative C5 factor.
    return np.array([dy,dl,quartic_factor*lam-3*gw+6*yt*yt])/LOOP


def trajectory(logpars,decoupling=1e6,steps=400,quartic_factor=4):
    a,b,M=np.exp(logpars); x=BOUND['T']/(3+BOUND['r'])
    S=(BOUND['T']-3*x)/BOUND['r']; Y=construct(a,b,M,S)
    H=Y.conj().T@Y; Q=np.trace(H@H).real
    lam=(3*x*x+BOUND['r']*Q)/BOUND['T']
    initial=np.r_[math.sqrt(x),lam,Y.ravel(),np.array([[0,1],[1,0]]).ravel()]
    td=math.log(decoupling/MU)
    high=rk4(full_rhs,initial,0,td,steps)
    Yd=high[2:8].reshape(3,2); Bd=high[8:12].reshape(2,2)
    C=Yd@np.linalg.solve(Bd,Yd.T)/M
    assert np.linalg.norm(C-C.T)<1e-25
    masses=M*np.linalg.svd(Bd,compute_uv=False)
    low=rk4(lambda t,z:eft_rhs(t,z,quartic_factor),[high[0],high[1],0],td,END,steps)
    factor=math.exp(low[2].real); Clow=factor*C
    light=np.linalg.svd(Clow,compute_uv=False)*VEV**2/2*1e9
    error_bound=2*BOUND['r']*M*M*a*b/BOUND['T']
    return dict(high_C5_singular_values_GeV_inverse=[a,b],high_Majorana_GeV=M,
        high_S=S,high_Q=Q,high_lambda=lam,high_yukawa_small_to_large_norm=float(np.linalg.norm(Y[:,1])/np.linalg.norm(Y[:,0])),
        analytic_decrease_of_high_quartic_bound=error_bound,
        threshold_GeV=decoupling,heavy_masses_at_threshold_GeV=masses.tolist(),
        relative_heavy_splitting=float((masses[0]-masses[1])/decoupling),
        relative_generated_Majorana_diagonal=float(np.linalg.norm(np.diag(Bd))/np.linalg.norm(Bd)),
        high_to_threshold_C5_matrix=C,low_C5_matrix=Clow,low_masses_eV=light.tolist(),
        C5_EFT_running_factor=factor,top_at_MZ=low[0].real,lambda_at_MZ=low[1].real,
        # Matching-scale diagnostic only: dim-6 running is not implemented.
        active_heavy_mixing_trace_at_threshold=float(VEV**2/(2*M*M)*np.linalg.norm(Yd@np.linalg.inv(Bd))**2))


def shoot(decoupling,steps=400):
    protected=old.trajectory(BOUND['gw_squared'],decoupling,steps)
    pars=np.log(np.r_[2*TARGET_MASSES*1e-9/VEV**2,protected['Majorana_at_matching_GeV']])
    def residual(p):
        row=trajectory(p,decoupling,steps)
        return np.r_[np.log(np.array(row['low_masses_eV'][:2])/TARGET_MASSES),
            math.log(math.sqrt(np.prod(row['heavy_masses_at_threshold_GeV']))/decoupling)]
    for iteration in range(8):
        err=residual(pars)
        if max(abs(err))<2e-11:break
        h=1e-4;J=np.column_stack([(residual(pars+h*np.eye(3)[j])-residual(pars-h*np.eye(3)[j]))/(2*h) for j in range(3)])
        pars-=np.linalg.solve(J,err)
    assert max(abs(residual(pars)))<2e-11
    return pars,iteration+1


def run():
    checks=[];rng=np.random.default_rng(539);errors=[]
    for k in range(20):
        U,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        a,b,M=0.14,0.07,1.2; S=M*(a+b)*(1+k/4)
        Y=construct(a,b,M,S,U); H=Y.conj().T@Y; MR=M*np.array([[0,1],[1,0]])
        C=Y@np.linalg.solve(MR,Y.T);target=a*np.outer(U[:,0],U[:,0])+b*np.outer(U[:,1],U[:,1])
        errors.append(float(np.linalg.norm(C-target)))
        assert np.allclose(C,target,atol=1e-14)
        assert abs(np.trace(H).real-S)<1e-13
        assert abs(np.linalg.det(H).real-M*M*a*b)<1e-13
        assert abs(np.trace(H@H).real-(S*S-2*M*M*a*b))<1e-13
    checks.append('rank_two_complex_construction_and_exact_common_moment_identity')

    for k in range(30):
        Y=rng.normal(size=(3,2))+1j*rng.normal(size=(3,2));M=1.7
        C=Y@Y.T/M;ab=np.linalg.svd(C,compute_uv=False)[:2];H=Y.conj().T@Y
        assert np.trace(H).real+1e-12>=M*sum(ab)
        assert abs(np.linalg.det(H).real-M*M*np.prod(ab))<2e-12
    checks.append('arbitrary_two_column_degenerate_basis_necessity_without_real_vector_assumption')

    T,r=1.7,2.3; cases=[]
    for M in (1,4,8):
        a,b=.06,.02; exact=sharp_min(a,b,M,r,T)
        lower=M*(a+b)
        if exact is None:continue
        sample=np.linspace(lower,T/r,10001);x=(T-r*sample)/3
        vals=(3*x*x+r*(sample**2-2*M*M*a*b))/T
        assert min(vals)>=exact-1e-13 and min(vals)-exact<1e-8
        cases.append(dict(M=M,sharp_minimum=exact,trace_floor=lower))
    assert sharp_min(.06,.02,10,r,T) is None
    checks.append('piecewise_sharp_Higgs_minimum_and_empty_domain')

    c=np.array([1,2j,3],complex);v=np.array([1j,2,1],complex);Y=np.outer(c,v)
    for k in range(8):
        A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));M=A+A.T+5*np.eye(3)
        sv=np.linalg.svd(Y@np.linalg.solve(M,Y.T),compute_uv=False)
        assert sv[1]<1e-13*sv[0]
    Y=construct(.01,.004,1,.4);H=Y.conj().T@Y;B=np.array([[0,1],[1,0]])
    assert np.linalg.norm(np.diag(H.T@B+B@H))>1e-4
    checks.append('Majorana_only_rank_one_obstruction_and_broken_pair_RG_not_kept_degenerate')

    rows=[];convergence=[];solver=[];saved_pars=[]
    for dec in (1e6,1e9):
        pars,it=shoot(dec,400);fine=trajectory(pars,dec,800);coarse=trajectory(pars,dec,400)
        rel=float(np.max(abs(np.array(fine['low_masses_eV'][:2])/TARGET_MASSES-1)))
        assert rel<2e-8 and fine['low_masses_eV'][2]<1e-16
        thresh=abs(math.sqrt(np.prod(fine['heavy_masses_at_threshold_GeV']))/dec-1)
        assert thresh<2e-10
        rows.append(fine);saved_pars.append(pars)
        solver.append(dict(threshold=dec,iterations=it,low_mass_relative_residual=rel,threshold_relative_residual=thresh))
        maxrel=max(abs(np.array(fine['low_masses_eV'][:2])/coarse['low_masses_eV'][:2]-1))
        assert maxrel<2e-8
        convergence.append(dict(threshold=dec,mass_step_halving_relative_change=float(maxrel),
            quartic_step_halving_change=abs(fine['lambda_at_MZ']-coarse['lambda_at_MZ'])))
    checks.append('two_nonzero_mass_targets_shot_through_full_matrix_RG_and_dim5_running')
    checks.append('same_parameters_step_halving_and_common_running_mass_threshold')

    controls=[]
    for dec,pars,row in zip((1e6,1e9),saved_pars,rows):
        wrong=trajectory(pars,dec,800,1)
        ratio=np.array(wrong['low_masses_eV'][:2])/row['low_masses_eV'][:2]
        assert min(ratio)>1.01 and abs(ratio[0]-ratio[1])<1e-12
        protected=old.trajectory(BOUND['gw_squared'],dec,800)
        assert abs(protected['lambda_at_MZ']-row['lambda_at_MZ'])<1e-9
        controls.append(dict(threshold=dec,wrong_quartic_coefficient_mass_ratio=ratio.tolist(),
            quartic_change_from_protected_538=row['lambda_at_MZ']-protected['lambda_at_MZ'],
            generated_mass_diagonal=row['relative_generated_Majorana_diagonal']))
    checks.append('wrong_Weinberg_quartic_coefficient_changes_nonzero_masses')
    checks.append('small_breaking_preserves_large_common_Higgs_boundary_numerically_not_exactly')

    def tidy(v):
        if isinstance(v,np.ndarray):
            if np.iscomplexobj(v):return dict(real=tidy(v.real.tolist()),imaginary=tidy(v.imag.tolist()))
            return tidy(v.tolist())
        if isinstance(v,(float,np.floating)):return float(format(float(v),'.13g'))
        if isinstance(v,list):return [tidy(x) for x in v]
        if isinstance(v,dict):return {k:tidy(x) for k,x in v.items()}
        return v
    deps=('research_note_538.md','protected_pair_rg.py','protected_pair_rg_results.json',
          'joint_gauge_matter_constraints_results.json','joint_yukawa_higgs_matching_results.json',
          'unified_physics_condition_ledger_538.md')
    return tidy(dict(round=539,tests_run=len(checks),failures=0,errors=0,checks=checks,
        construction_max_matrix_error=max(errors),sharp_minimum_examples=cases,
        target_mass_input=dict(source='NuFIT 6.0 (2024), Table 1, IC24 with SK, NO; declared m1=0',
            source_url='https://arxiv.org/html/2410.05380',delta_m31_eV2=2.513e-3,delta_m21_eV2=7.49e-5,
            m3_m2_eV=TARGET_MASSES.tolist(),tree_vev_conversion_GeV=VEV,
            not_current_global_fit=True),solver=solver,step_halving=convergence,
        trajectories=rows,controls=controls,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(common_spectral_boundary_exact_at_high_scale=True,
            degenerate_heavy_boundary_not_reimposed_during_RG=True,
            dimension_four_full_matrix_one_loop_and_dim5_EFT_running=True,
            common_scale_tree_matching_near_degenerate_pair=True,
            mass_splittings_used_as_targets_not_predictions=True,
            nonzero_charged_flavor_PMNS_or_pole_fit_completed=False,
            dimension_six_running_and_finite_one_loop_matching_included=False,
            spacetime_GR_or_full_unification_derived=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','solver','controls')},ensure_ascii=False))
