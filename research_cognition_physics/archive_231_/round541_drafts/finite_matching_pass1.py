"""Round 541: genuine finite matching in the exact lepton-number pair.

Zhang--Zhou 2107.12133 v3 Eqs(125),(126),(128), in positive heavy mass basis.
Dimension-four one-loop flow with the displayed matching through m^2/M^2;
not a complete dimension-six flow or a calculation of the spectral mass boundary.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import protected_pair_rg as old
import joint_higgs_flow_obstruction as obstruction

HERE=Path(__file__).resolve().parent
TARGET=HERE/'protected_pair_finite_matching_results.json'
LOOP=16*math.pi**2
HIST,BOUND=old.inputs()
XSTAR=HIST['matched_inverse_couplings'];MU=BOUND['matching_scale_GeV']
IR=173.34
MASS_TARGET=-(131.55**2)/2  # Historical MS convention conversion, not a prediction.


def finite_matching(yt,lam,S,gw,xi):
    dl=(-lam*S+xi*((5*gw/18-4*lam/3)*S+S*S))/LOOP
    dy=-yt*S*(.25+xi/3)/LOOP
    eff=xi+S/LOOP*(2-xi/2-xi*xi/3)
    return yt+dy,lam+dl,eff


def inverse_mass(S,eff):
    """Rationalized root near zero; keep cancellation diagnostic in Decimal."""
    with localcontext() as ctx:
        ctx.prec=60
        e=Decimal(str(S))/Decimal(str(LOOP)); target=Decimal(str(eff))
        A=-e/3;B=1-e/2;C=2*e-target
        if e==0:return target,Decimal(0),None
        disc=(B*B-4*A*C).sqrt()
        root=-2*C/(B+disc);other=(-B-disc)/(2*A)
        recovered=root+e*(2-root/2-root*root/3)
        assert abs(recovered-target)<Decimal('1e-55')
        return root,recovered,other


def low_rhs(t,z):
    y,l=z[:2];gy,gw,gc=old.gauge_squared(t,XSTAR)
    dy=y*(4.5*y*y-17*gy/12-9*gw/4-8*gc)
    dl=24*l*l-(3*gy+9*gw)*l+3/8*(gy*gy+2*gy*gw+3*gw*gw)+12*y*y*l-6*y**4
    dm=-1.5*gy-4.5*gw+12*l+6*y*y
    return np.array([dy,dl,dm])/LOOP


def low_flow(yt,lam,td,steps):
    z=np.array([yt,lam,0.]);end=math.log(IR/MU);h=(end-td)/steps
    for j in range(steps):
        t=td+j*h;k1=low_rhs(t,z);k2=low_rhs(t+h/2,z+h*k1/2)
        k3=low_rhs(t+h/2,z+h*k2/2);k4=low_rhs(t+h,z+h*k3)
        z+=h*(k1+2*k2+2*k3+k4)/6
    assert np.all(np.isfinite(z)) and z[0]>=0 and z[1]>0
    return z


def trajectory(M,steps=600):
    x=BOUND['T']/(3+BOUND['r']);S0=(BOUND['T']-3*x)/BOUND['r']
    lam0=(3*x*x+BOUND['r']*S0*S0)/BOUND['T'];td=math.log(M/MU)
    high,_,_=old.evolve([math.sqrt(x),math.sqrt(S0),lam0,0.,1.],0,td,XSTAR,True,steps)
    yt,ss,lam=high[:3];S=ss*ss;gw=old.gauge_squared(td,XSTAR)[1]
    # Fit only the threshold mass parameter; all dimensionless high inputs stay fixed.
    xi=0.
    for iterations in range(20):
        yplus,lplus,_=finite_matching(yt,lam,S,gw,xi)
        low=low_flow(yplus,lplus,td,steps)
        eff=MASS_TARGET/(M*M*math.exp(low[2]))
        root,recovered,other=inverse_mass(S,eff)
        new=float(root)
        if abs(new-xi)<1e-15:xi=new;break
        xi=new
    else:raise AssertionError('mass boundary iteration failed')
    assert abs(xi)<.01
    yplus,lplus,_=finite_matching(yt,lam,S,gw,xi);low=low_flow(yplus,lplus,td,steps)
    eff=MASS_TARGET/(M*M*math.exp(low[2]));root,recovered,other=inverse_mass(S,eff)
    assert abs(float(root)-xi)<1e-15
    control=low_flow(yt,lam,td,steps)
    zero_y,zero_l,zero_eff=finite_matching(yt,lam,S,gw,0.)
    zero_low=low_flow(zero_y,zero_l,td,steps)
    return dict(threshold_GeV=M,Majorana_at_high_GeV=M*math.exp(-high[3]),
        S_at_threshold=S,lambda_before=lam,lambda_after=lplus,
        delta_lambda=lplus-lam,delta_top=yplus-yt,
        low_top=low[0],low_lambda=low[1],low_mass_squared_target_GeV2=MASS_TARGET,
        low_mass_running_factor=math.exp(low[2]),
        mass_ratio_UV_at_threshold=str(root),mass_ratio_EFT_at_threshold=str(recovered),
        large_other_root_outside_hierarchy=str(other),
        inverse_mass_equation_decimal_precision=60,
        mass_counterterm_input_GeV2=float(root)*M*M,
        leading_heavy_mass_shift_GeV2=2*M*M*S/LOOP,
        retained_light_over_leading_heavy_mass_ratio=abs(eff)/(2*S/LOOP),
        zero_UV_mass_control_low_GeV2=zero_eff*M*M*math.exp(zero_low[2]),
        tree_matching_control_low_lambda=control[1],
        low_quartic_change_from_tree_matching=low[1]-control[1],
        iterations=iterations+1)


def run():
    checks=[];rng=np.random.default_rng(541);matrix_error=0.
    for _ in range(20):
        c=rng.normal(size=3)+1j*rng.normal(size=3);c*=math.sqrt(.4)/np.linalg.norm(c)
        Y=np.column_stack([c,1j*c,np.zeros(3)])/math.sqrt(2);G=Y.conj().T@Y
        P=np.trace(G@G.T);Q=np.trace(G@G)
        matrix_error=max(matrix_error,float(abs(P)))
        assert abs(P)<1e-15 and abs(Q-.16)<1e-14
        assert np.linalg.norm(Y@Y.T)<1e-15
        # Positive-heavy-mass basis is preserved by real orthogonal rotations.
        O,_=np.linalg.qr(rng.normal(size=(3,3)));Gr=(Y@O).conj().T@(Y@O)
        assert abs(np.trace(Gr@Gr.T))<1e-15
    checks.append('protected_transpose_invariant_zero_but_Hermitian_trace_nonzero')

    # Reduce the source formula with exact rational arithmetic, including xi terms.
    l,S,xi,gw=F(3,10),F(2,5),F(-1,200),F(2,5)
    source=((5*gw*xi-18*l-24*l*xi)/18)*S+xi*S*S
    reduced=-l*S+xi*((F(5,18)*gw-F(4,3)*l)*S+S*S)
    assert source==reduced
    assert float(finite_matching(1.,float(l),float(S),float(gw),float(xi))[1]-float(l))<0
    checks.append('source_finite_quartic_reduction_in_same_normalization')

    # Strict endpoint barriers, using the already certified historical gauge enclosure.
    _,bmax,_,_=obstruction.exact_gauge_enclosures(HIST,BOUND)
    assert bmax<F(431,100)
    # Conservative gY²<.162, gw²<.425, gc²<1.489, verified from frozen endpoints.
    gs=old.gauge_squared(math.log(HIST['MZ_GeV']/MU),XSTAR)
    assert max(old.gauge_squared(0,XSTAR)[0],gs[0])<.162 and gs[1]<.425 and gs[2]<1.489
    Amax=F(17,6)*F(162,1000)+F(9,2)*F(425,1000)+16*F(1489,1000)
    assert Amax<27
    Smax=F(739,1000);r=F(str(BOUND['r']));T=F(str(BOUND['T']))
    assert T/r<Smax and T/3<3
    assert F(431,200)-5*Smax<0
    assert -24+F(431,100)+18<0
    checks.append('analytic_S_top_and_quartic_barriers_before_the_real_threshold')

    rho=F(1,100);jump=(Smax+rho*((F(5,18)*F(425,1000)+F(4,3))*Smax+Smax*Smax))/157
    assert jump<F(1,200) and jump<F(881,30000)
    refined=F(17009,100000)-F(3,2)*F(1,200)
    assert refined==F(16259,100000) and refined>F(12604,100000)
    assert Smax*(F(1,4)+rho/3)/157<1
    checks.append('strict_finite_jump_bound_is_insufficient_for_round540_repair_budget')

    for s in (.001,.3,.738):
        for eff in (-1e-6,-1e-10,-1e-16):
            root,recovered,other=inverse_mass(s,eff)
            assert abs(root)<Decimal('.01') and other>100
            assert abs(float(recovered)-eff)<abs(eff)*1e-12
    checks.append('stable_mass_inverse_and_extraneous_large_root_exclusion')

    rows=[];steps=[]
    for M in (1e4,1e6,1e9):
        coarse=trajectory(M,300);fine=trajectory(M,600)
        delta=max(abs(coarse[k]-fine[k]) for k in ('low_top','low_lambda','lambda_after'))
        assert delta<1e-9 and fine['low_lambda']>.16259
        assert fine['zero_UV_mass_control_low_GeV2']>0 and fine['mass_counterterm_input_GeV2']<0
        assert fine['delta_lambda']<0 and fine['delta_top']<0
        rows.append(fine);steps.append(dict(M=M,max_step_halving_change=delta))
    checks.append('joint_top_Higgs_finite_matching_and_declared_mass_boundary_solutions')
    checks.append('fixed_parameters_step_halving_and_zero_UV_mass_counterexample')

    deps=('research_note_540.md','joint_higgs_flow_obstruction.py','joint_higgs_flow_obstruction_results.json',
        'protected_pair_rg.py','protected_pair_rg_results.json','joint_yukawa_higgs_matching_results.json',
        'unified_physics_condition_ledger_540.md')
    return dict(round=541,tests_run=len(checks),failures=0,errors=0,checks=checks,
        transpose_invariant_max_error=matrix_error,
        analytic_global_bounds=dict(S_max=str(Smax),top_squared_max=3,quartic_before_max=1,
            abs_mass_ratio_max=str(rho),exact_negative_jump_upper=str(jump),
            negative_jump_upper=float(jump),safe_jump_upper=.005,
            low_quartic_floor=.16259,required_round540_budget=float(F(881,30000))),
        mass_input=dict(source='Buttazzo 1307.3536 Eq(56), historical MS m_B=131.55 GeV',
            conversion='V=-(m_B^2/2) HdagH + lambda(HdagH)^2; project m^2=-m_B^2/2',
            low_scale_GeV=IR,mass_squared_GeV2=MASS_TARGET,
            spectral_high_scale_mass_boundary_imposed=False),
        trajectories=rows,step_halving=steps,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(exact_lepton_number_pair_only=True,small_broken_nondegenerate_extension_proved=False,
            finite_one_loop_matching_for_displayed_couplings_included=True,
            truncated_inverse_heavy_mass_expansion=True,full_dimension_six_running_included=False,
            mass_boundary_is_reverse_input_not_prediction=True,all_top_threshold_choices_under_declared_bounds=True,
            physical_global_fit_or_full_unification_completed=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','analytic_global_bounds','step_halving')},ensure_ascii=False))
