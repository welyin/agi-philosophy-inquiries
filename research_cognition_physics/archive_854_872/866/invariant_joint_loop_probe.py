"""Working 866: a positive joint class-function readout on the original G.
No formal866 count. Finite graph instrument; no continuum CP lift asserted.
SU(3) Haar integration uses the Weyl torus formula, and the U(1) cover
integrates quotient-invariant functions. SU(2) is a spectator.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'invariant_joint_loop_probe_results.json'
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(HERE.parent/'865'))
from research_layout import Layout,ResearchRuntime
import relational_loop_fusion_matching as original

ETA_R=1/30
ETA_A=1/320
LOWER=1-18*ETA_R-64*ETA_A

def menu_grid(N):
    t=2*np.pi*np.arange(N)/N
    a,b,z=np.meshgrid(t,t,t,indexing='ij')
    u,v=np.exp(1j*a),np.exp(1j*b);w=np.exp(-1j*(a+b))
    chi=u+v+w
    weight=abs((u-v)*(u-w)*(v-w))**2/(6*N**3)
    return np.exp(-2j*z)*chi,abs(chi)**2-1,weight

def integrate(U,dU,T,N,theta=.37):
    cy,ay,weight=menu_grid(N)
    phase=np.exp(-2j*theta);color=np.trace(U);dc=np.trace(dU)
    x=phase*color;dx=phase*dc
    a=abs(color)**2-1;da=2*np.real(np.conj(color)*dc)
    K=1+2*ETA_R*np.real(x*np.conj(cy))+ETA_A*a*ay
    dK=2*ETA_R*np.real(dx*np.conj(cy))+ETA_A*da*ay
    assert float(np.min(K))>=LOWER-1e-14
    mean=lambda f:np.sum(weight*K*f)
    naive=abs(cy/ETA_R)**2
    corrected=1+ay/ETA_A
    variance=float(mean(naive).real-abs(x)**2)
    variance_exact=1/ETA_R**2-1+(ETA_A/ETA_R**2-1)*a
    derivative=float(np.sum(weight*dK*naive)-2*np.real(np.conj(x)*dx))
    derivative_exact=float((ETA_A/ETA_R**2-1)*da)
    errors=[abs(np.sum(weight)-1),abs(mean(1)-1),
        abs(mean(cy/ETA_R)-x),abs(mean(corrected)-(1+a)),
        abs(variance-variance_exact),abs(derivative-derivative_exact)]
    assert max(errors)<2e-10
    color_dirichlet=0.
    for generator in T:
        du=1j*generator@U;dc_i=np.trace(du);dx_i=phase*dc_i
        da_i=2*np.real(np.conj(color)*dc_i)
        dk=2*ETA_R*np.real(dx_i*np.conj(cy))+ETA_A*da_i*ay
        color_dirichlet+=float(np.sum(weight*dk*dk/(4*K)))
    dx0=-2j*x
    dk0=2*ETA_R*np.real(dx0*np.conj(cy))
    abelian_dirichlet=float(np.sum(weight*dk0*dk0/(4*K)))
    tangent_dirichlet=float(np.sum(weight*dK*dK/(4*K)))
    assert 0<color_dirichlet<1.5 and 0<abelian_dirichlet<1.8 and tangent_dirichlet>0
    return dict(N=N,kernel_minimum_sample=float(np.min(K)),
        polynomial_moment_error=float(max(errors)),
        added_complex_variance=variance,added_complex_variance_exact=float(variance_exact),
        physical_tangent_of_variance=derivative,
        physical_tangent_of_variance_exact=derivative_exact,
        source_tangent_dirichlet=tangent_dirichlet,
        color_Dirichlet_per_link=color_dirichlet,
        U1_Dirichlet_per_link=abelian_dirichlet)

def run():
    assert abs(LOWER-.2)<1e-14
    c=original.inherited.constants()
    rows=[]
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        for L in (.4,.2):
            U,dU,_,_=original.matrices(old,c,L)
            lo=integrate(U,dU,old.T,48)
            hi=integrate(U,dU,old.T,72)
            diffs={key:abs(hi[key]-lo[key]) for key in (
                'color_Dirichlet_per_link','U1_Dirichlet_per_link','source_tangent_dirichlet')}
            assert max(diffs.values())<1e-12
            rows.append(dict(original_coordinate_side=L,lower_resolution=lo,higher_resolution=hi,
                             Dirichlet_quadrature_differences=diffs))
    return dict(date='2026-10-06',status='866_working_not_formal',
        formal_reports=865,cumulative_numbered_groups=3650,fresh_numbered_groups=0,
        all_diagnostic_checks_passed=True,
        input='fixed dimensionless finite-resolution kernel; not derived bulk dynamics',
        kernel='1 + 2 eta_R Re(chi_R(g) conj(chi_R(y))) + eta_A chi_A(g) chi_A(y)',
        eta_R=ETA_R,eta_A=ETA_A,analytic_uniform_lower_bound=LOWER,
        exact_group='(SU3 x SU2 x U1) / Z6',
        branchwise_Gauss_preservation='analytic class-function multiplication Kraus proof',
        independent_moment_check='SU3 Weyl integration and U1 cover quadrature',
        rows=rows,
        energy_bound_scope='old fixed finite graph electric form, hbar=1; not full continuum source energy',
        original_continuum_CP_instrument_constructed=False,
        full_original_dynamics_match_proved=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
