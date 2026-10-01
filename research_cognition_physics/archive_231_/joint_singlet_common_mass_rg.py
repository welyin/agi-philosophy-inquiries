"""544: common UV masses and a protected pair in an explicit singlet extension.

One-loop unbroken MS full-theory flow. Fixed-scale running tree potential only;
no decoupling, Coleman-Weinberg vacuum, pole fit or universal hierarchy claim.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from protected_pair_rg import inputs,gauge_squared,LOOP

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_singlet_common_mass_rg_results.json'
OLD,BOUND=inputs()
T,R=BOUND['T'],BOUND['r']
XSTAR=OLD['matched_inverse_couplings']
MU=BOUND['matching_scale_GeV']


def initial(q0):
    q0=np.asarray(q0,dtype=float); S=(T-3*q0)/R
    assert np.all(S>=-1e-15) and np.all(q0>=0)
    return np.array([q0,S,np.full_like(q0,T/(2*R)),(3*q0*q0+R*S*S)/T,
                     S,np.full_like(q0,T/R),np.ones_like(q0),np.ones_like(q0)])


def beta_numerator(a,g):
    """L d/dln(mu), masses represented by x=-mH²/C, y=-mS²/C."""
    q,S,z,lh,p,ls,x,y=a; gy,gw,gc=g
    B=3*gy+9*gw; C=3/8*(gy*gy+2*gy*gw+3*gw*gw)
    return [2*q*(4.5*q+S-17*gy/12-9*gw/4-8*gc),
        S*(6*q+5*S+z-B/2),2*z*(5*z+S),
        24*lh*lh-B*lh+C+4*(3*q+S)*lh-2*(3*q*q+S*S)+2*p*p,
        p*(12*lh+6*ls+8*p+6*q+2*S+4*z-B/2)-4*z*S,
        18*ls*ls+8*ls*z+8*p*p-8*z*z,
        (12*lh+6*q+2*S-B/2)*x+2*p*y,(6*ls+4*z)*y+8*p*x]


def rhs(u,a):
    return -np.array(beta_numerator(a,gauge_squared(-u,XSTAR)))/LOOP


def flow(q0,u,steps=800,mass_scale=1.):
    a=initial(q0);a[6:]*=mass_scale;du=u/steps
    for j in range(steps):
        t=j*du;k1=rhs(t,a);k2=rhs(t+du/2,a+du*k1/2)
        k3=rhs(t+du/2,a+du*k2/2);k4=rhs(t+du,a+du*k3)
        a+=du*(k1+2*k2+2*k3+k4)/6
    assert np.all(np.isfinite(a))
    return a


def invariants(a):
    q,S,z,lh,p,ls,x,y=a
    return lh*ls-p*p,ls*x-p*y,lh*y-p*x


def diagnostics(a):
    q,S,z,lh,p,ls,x,y=map(float,a);D,A,B=invariants(a)
    out=dict(state=dict(zip(('q','S','z','lambda_H','p','lambda_s','x','y'),map(float,a))),
             determinant=float(D),H_numerator=float(A),singlet_numerator=float(B),
             strict_double_vev=bool(lh>0 and ls>0 and D>0 and A>0 and B>0))
    if out['strict_double_vev']:
        ratio=B/A;h2=A/D;s2=B/D
        massmat=2*np.array([[lh,p*math.sqrt(ratio)],[p*math.sqrt(ratio),ls*ratio]])
        vals=np.linalg.eigvalsh(massmat)
        assert min(vals)>0
        out.update(h_squared_over_common_scale=float(h2),s_squared_over_common_scale=float(s2),
                   squared_vev_ratio=float(ratio),Majorana_squared_over_h_squared=float(z*ratio),
                   radial_squared_over_h_squared=vals.tolist(),formal_Schur_lambda=float(D/ls),
                   tadpole_residual=float(max(abs(lh*h2+p*s2-x),abs(p*h2+ls*s2-y))))
    return out


def run():
    checks=[]
    # Exact conversion of source quartic/mass coefficients, after YN=k/2.
    q,S,z,lh,p,ls,x,y,gy,gw,gc=map(F,('0.2','0.3','0.4','0.25','0.15','0.7','0.8','0.6','0.16','0.32','0.9'))
    Hs,Ps,Ls=2*lh,2*p,2*ls;trN2=z/2;trN4=z*z/8;trNG=z*S/4
    scalar_source=[(Hs*(12*Hs-9*gw-3*gy+12*q+4*S)+F(9,4)*gw*gw+
        F(3,2)*gw*gy+F(3,4)*gy*gy-12*q*q+Ps*Ps-4*S*S)/2,
        (Ps*(6*Hs-F(9,2)*gw-F(3,2)*gy+6*q+4*Ps+3*Ls+2*S+8*trN2)-32*trNG)/2,
        (Ls*(9*Ls+16*trN2)+4*Ps*Ps-128*trN4)/2,
        (6*Hs-F(9,2)*gw-F(3,2)*gy+6*q+2*S)*x+Ps*y,
        (3*Ls+8*trN2)*y+4*Ps*x]
    exact=[24*lh*lh-(3*gy+9*gw)*lh+F(3,8)*(gy*gy+2*gy*gw+3*gw*gw)+
           4*(3*q+S)*lh-2*(3*q*q+S*S)+2*p*p,
           p*(12*lh+6*ls+8*p+6*q+2*S+4*z-F(3,2)*gy-F(9,2)*gw)-4*z*S,
           18*ls*ls+8*ls*z+8*p*p-8*z*z,
           (12*lh+6*q+2*S-F(3,2)*gy-F(9,2)*gw)*x+2*p*y,
           (6*ls+4*z)*y+8*p*x]
    assert scalar_source==exact
    assert np.allclose(beta_numerator(list(map(float,(q,S,z,lh,p,ls,x,y))),
        list(map(float,(gy,gw,gc))))[3:],list(map(float,exact)),atol=1e-14)
    checks.append('exact_source_to_canonical_quartic_and_mass_normalization')

    rng=np.random.default_rng(544);matrix_error=0.;spectrum_error=0.
    for _ in range(12):
        c=rng.normal(size=3)+1j*rng.normal(size=3);c*=math.sqrt(.31)/np.linalg.norm(c)
        Y=np.column_stack((c,1j*c))/math.sqrt(2);G=Y.conj().T@Y;K=math.sqrt(.37)*np.eye(2)
        z=.37;S=.31;k=math.sqrt(z)
        # Symmetric K beta: self term plus G.T K + K G; G+G.T=S I.
        dK=5*z*K+G.T@K+K@G
        dY=(3*.2+S-.75*.16-2.25*.32)*Y+1.5*Y@G+.5*Y@(K.conj().T@K)
        expected=(3*.2+2.5*S+.5*z-.75*.16-2.25*.32)*Y
        matrix_error=max(matrix_error,float(np.max(abs(dK-(5*z+S)*K))),
                         float(np.max(abs(dY-expected))))
        inv=np.linalg.inv(K);zero=Y@inv@Y.T
        tangent=dY@inv@Y.T+Y@inv@dY.T-Y@inv@dK@inv@Y.T
        assert np.linalg.norm(zero)<1e-14 and np.linalg.norm(tangent)<1e-14
        h=.7;s=1.1;M=np.zeros((5,5),complex)
        M[:3,3:]=Y*h/math.sqrt(2);M[3:,:3]=M[:3,3:].T;M[3:,3:]=s*K
        squared=np.linalg.eigvalsh(M.conj().T@M)
        target=np.array([0,0,0,z*s*s+S*h*h/2,z*s*s+S*h*h/2])
        spectrum_error=max(spectrum_error,float(np.max(abs(squared-target))))
        assert np.allclose(squared,target,atol=1e-14)
    assert matrix_error<1e-14 and spectrum_error<1e-14
    checks.append('complex_pair_matrix_closure_zero_Weinberg_and_independent_fermion_box_spectrum')

    # Hard Majorana and the two odd couplings have identically zero RHS at zero.
    YN=np.eye(2)*.3;M=np.zeros((2,2));G=np.array([[.2,.2j],[-.2j,.2]])
    hard=M@G.T+G@M+4*np.trace(M@YN)*YN+12*M@YN@YN
    oddS=0.-96*np.trace(M@np.linalg.matrix_power(YN,3))
    oddHS=0.-16*np.trace(M@YN@G)
    assert np.max(abs(hard))==0 and oddS==oddHS==0
    # A nonzero hard mass would source an odd term: not a manual projection.
    assert -96*np.trace(np.eye(2)@np.linalg.matrix_power(YN,3))!=0
    checks.append('hard_mass_and_odd_zero_branch_is_tangent_not_manually_reset')

    # Exact balanced-point derivative in arbitrary rational r,w,gauge values.
    rr,w,gy,gw=map(F,('2','0.3','0.16','0.32'))
    ls=(3+rr)*w/rr;z=ls/2;B=3*gy+9*gw;Cg=F(3,8)*(gy*gy+2*gy*gw+3*gw*gw)
    blh=34*w*w-B*w+Cg;bp=w*(28*w+6*ls-B/2)
    dx=-(22*w-B/2);dy=-(8*ls+8*w)
    dG=-blh+w*dy+bp-w*dx
    assert dG==6*(1-1/rr)*w*w-Cg
    w=T/(3+R);g=gauge_squared(0,XSTAR);Cg=3/8*(g[0]**2+2*g[0]*g[1]+3*g[1]**2)
    coefficient=6*(1-1/R)*w*w-Cg;F0=3*w/R
    assert coefficient>0
    slope=coefficient/(LOOP*F0)
    tiny=flow(w,1e-4,16);tiny_d=diagnostics(tiny)
    assert tiny_d['strict_double_vev']
    assert abs(tiny_d['squared_vev_ratio']/1e-4-slope)<2e-7
    checks.append('exact_balanced_point_lifting_and_positive_small_scale_double_condensation')

    end=math.log(MU/173.34);rows=[];convergence=0.
    for q0,u in ((w,1.),(w,10.),(w,end),(.25,end),(.5,end),(T/3,end)):
        a=flow(q0,u,800);fine=flow(q0,u,1600)
        err=float(np.max(abs(a-fine)));convergence=max(convergence,err)
        assert err<2e-10
        row=diagnostics(fine);row.update(q0=float(q0),u=float(u),mu_GeV=MU*math.exp(-u),
                                       step_doubling_max_error=err)
        assert row['strict_double_vev'] and row['tadpole_residual']<1e-13
        rows.append(row)
    checks.append('coupled_eight_variable_flow_step_doubling_and_vacuum_residuals')

    connected=rows[4]
    assert connected['state']['S']>0 and connected['state']['p']>0
    assert connected['squared_vev_ratio']>R/3
    assert connected['Majorana_squared_over_h_squared']>T/6
    checks.append('nonzero_portal_running_counterexample_to_extending_bare_scale_bounds')

    a=flow(w,end,800);scaled=flow(w,end,800,7.)
    assert np.max(abs(scaled[:6]-a[:6]))==0
    assert np.allclose(scaled[6:],7*a[6:],rtol=2e-14)
    assert math.isclose(diagnostics(a)['squared_vev_ratio'],diagnostics(scaled)['squared_vev_ratio'],rel_tol=2e-14)
    checks.append('common_mass_rescaling_cannot_tune_the_dimensionless_hierarchy')

    # Exploratory grid, explicitly not a continuum proof or exclusion certificate.
    qs=np.linspace(0,T/3,257);grid=[]
    for u in (1.,5.,10.,15.,end):
        a=flow(qs,u,800);D,Fh,Gs=invariants(a)
        mask=(D>0)&(Fh>0)&(Gs>0)&(a[3]>0)&(a[5]>0)
        ratios=Gs[mask]/Fh[mask]
        assert np.all(np.isfinite(ratios)) and np.all(ratios>0)
        grid.append(dict(u=float(u),samples=len(qs),stable_double_samples=int(mask.sum()),
          max_squared_vev_ratio=float(max(ratios)),min_H_numerator=float(min(Fh)),
          max_Majorana_squared_over_h_squared=float(max(a[2,mask]*ratios)),
          no_continuum_exclusion_claim=True))
    assert rhs(0,initial(T/3))[1]==rhs(0,initial(T/3))[4]==0
    checks.append('finite_grid_phase_diagnostics_and_exact_decoupled_endpoint_control')

    deps=('research_note_543.md','joint_singlet_origin_and_vacuum.py',
        'joint_singlet_origin_and_vacuum_results.json','protected_pair_rg.py',
        'joint_gauge_matter_constraints_results.json','joint_yukawa_higgs_matching_results.json')
    return dict(round=544,tests_run=len(checks),failures=0,errors=0,checks=checks,
        source='https://arxiv.org/html/1608.00087',source_equations='HTML 53-63; v2 PDF 39-49',
        matrix_closure_max_error=matrix_error,fermion_spectrum_max_error=spectrum_error,
        common_boundary=dict(T=T,r=R,mu_GeV=MU,normalized_negative_masses=[1,1]),
        balanced_lifting=dict(w=w,gauge_box=Cg,L_times_G_prime=coefficient,
           squared_vev_ratio_slope=slope,tiny_u=1e-4,tiny_ratio=tiny_d['squared_vev_ratio']),
        examples=rows,max_step_doubling_error=convergence,finite_grid=grid,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(singlet_is_extra_input=True,common_MS_UV_mass_boundary_is_extra_matching_input=True,
          one_loop_full_theory_only=True,unbroken_parameter_running_tree_potential_only=True,
          wavefunction_improved_Coleman_Weinberg_vacuum=False,physical_decoupling_and_pole_fit=False,
          all_family_hierarchy_obstruction_proved=False,all_physics_derived=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','balanced_lifting','max_step_doubling_error')},ensure_ascii=False))
