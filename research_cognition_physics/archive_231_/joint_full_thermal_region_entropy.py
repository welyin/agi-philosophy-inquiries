"""637: full-model regional entropy certificate and quotient heat counting.

The full H Gibbs spectrum is NOT numerically evaluated. Its regional entropy
existence follows from the original analytic coercive bound and Haar cutting.
Numerics check the actual quotient reference spectrum, geometric constants,
and entropy bounds on the already existing original-space loop states.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quotient_boundary_entropy as boundary
import joint_thermal_gauss_source as thermal
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_full_thermal_region_entropy_results.json'


def sequence_bound(alpha,k,N):
    n=np.arange(N+1,dtype=float)
    partial=float(np.sum((n+1)**k*np.exp(-alpha*n*n)))
    first=(N+2.)**k*np.exp(-alpha*(N+1)**2)
    ratio=((N+3.)/(N+2.))**k*np.exp(-alpha*(2*N+3))
    assert ratio<1
    tail=float(first/(1-ratio))
    return partial+tail,tail


def rectangle_tail(s,N):
    f,ft=sequence_bound(s/3,4,N)
    w,wt=sequence_bound(s/4,2,N)
    q,qt=sequence_bound(s,0,N)
    q=2*q-1;qt=2*qt
    # Union bound on the four label coordinates, before imposing Z6.
    return float(2*ft*f*w*q+f*f*wt*q+f*f*w*qt)


def product_upper(s):
    N=max(16,int(np.ceil(np.sqrt(160/s))))
    f,_=sequence_bound(s/3,4,N)
    w,_=sequence_bound(s/4,2,N)
    q,_=sequence_bound(s,0,N)
    return f*f*w*(2*q-1)


def quotient_heat(s,N=32):
    x=np.arange(N+1);a,b=np.meshgrid(x,x,indexing='ij')
    dim3=(a+1)*(b+1)*(a+b+2)/2
    C3=(a*a+b*b+a*b+3*a+3*b)/3
    ell=x;C2=ell*(ell+2)/4
    color=dim3**2*np.exp(-s*C3)
    weak=(ell+1)**2*np.exp(-s*C2)
    q=np.arange(-N,N+1);qw=np.exp(-s*q*q)
    theta=np.array([qw[q%6==r].sum() for r in range(6)])
    residue=(-2*(a+2*b)[...,None]-3*ell)%6
    direct=float(np.sum(color[...,None]*weak*theta[residue]))
    fourier=0j
    for k in range(6):
        z3=np.sum(color*np.exp(2j*np.pi*k*(a+2*b)/3))
        z2=np.sum(weak*np.exp(1j*np.pi*k*ell))
        z1=np.sum(qw*np.exp(1j*np.pi*k*q/3))
        fourier+=z3*z2*z1/6
    direct_product=float(color.sum()*weak.sum()*qw.sum())
    return direct,fourier,direct_product


def quotient_heat_check():
    rows=[]
    for s in (.3,.6,1.2):
        z,f,prod=quotient_heat(s)
        large,_,_=quotient_heat(s,40)
        tail=rectangle_tail(s,32)
        assert abs(z-f.real)<2e-12*max(1,z)
        assert abs(f.imag)<2e-12*max(1,z)
        assert abs(large-z)<=tail+2e-12*max(1,z)
        assert z>=1 and z<prod
        assert z+tail<=product_upper(s)*(1+1e-13)
        rows.append(dict(s=s,quotient_heat_box=z,Fourier_projection_real=float(f.real),
            Fourier_projection_imaginary=float(f.imag),direct_product_heat=prod,
            rigorous_label_tail_bound=tail,larger_box_difference=large-z,
            unrestricted_polynomial_upper=product_upper(s)))
    return dict(rows=rows,degeneracy_is_squared_original_representation_dimension=True,
                quotient_projection_not_replaced_by_one_sixth_of_product=True,
                floating_evaluation_not_interval_arithmetic=True)


def scalar_integral_upper(alpha):
    a,b,R0=thermal.constants()
    x=R0/np.sqrt(6);omega=8*np.pi**2/3
    core=omega*36*np.sqrt(6)/8*(np.sinh(4*x)/4-2*np.sinh(2*x)+3*x)
    tail=9*omega/(4*b*alpha*a)*np.exp(-alpha*a*np.exp(b*R0))
    return float(core+tail),float(core)


def region_check():
    eps=.73;sigma=.12;t=.8
    S=np.array([[.24,.31,-.12],[.31,-.07,.17],[-.12,.17,-.17]])
    shape=boundary.gluing.geometry.shape_exp(S,t)
    bs=np.array(boundary.gluing.geometry.PAR['b'])
    electric_min=float(bs.min()/eps*np.exp(-2*sigma)*np.linalg.eigvalsh(shape).min())
    scalar_coefficient=float(1/(2*eps**3*np.exp(6*sigma)))
    cT=min(electric_min,scalar_coefficient)
    w=float(eps**3*np.exp(6*sigma));u=w/2
    # Same 3-neighbour periodic cubic prescription, N=3; A is the x=0 slab.
    vertices=list(np.ndindex(3,3,3));A={v for v in vertices if v[0]==0}
    inside=cross=0
    for v in vertices:
        for direction in range(3):
            v2=list(v);v2[direction]=(v2[direction]+1)%3;v2=tuple(v2)
            if v in A and v2 in A:inside+=1
            if (v in A)!=(v2 in A):cross+=1
    assert (len(A),inside,cross)==(9,18,18)
    rows=[]
    for s in (.3,.6,1.2):
        eta=s/cT
        integral,core=scalar_integral_upper(eta*u)
        heat=(4*np.pi*s)**-2.5*np.exp(-2*s/3)*(1+s/9)
        zg=quotient_heat(s)[0]+rectangle_tail(s,32)
        zh=quotient_heat(s/2)[0]+rectangle_tail(s/2,32)
        logZA=(32*len(A)*np.log(2)+len(A)*np.log(heat*integral)
               +inside*np.log(zg)+cross*np.log(zh))
        assert np.isfinite(logZA) and cT>0
        rows.append(dict(reference_heat_s=s,entropy_energy_slope_eta=eta,
            scalar_integral_upper=integral,scalar_core_volume=core,
            log_regional_reference_trace_upper=float(logZA),
            physical_full_thermal_energy_not_computed=True))
    # Original actual loop states from 636: the open arc costs 2*cT*C.
    C=np.array([0.,3.,6.,8.]);dims=np.array([1,8,10,27])
    delta=np.array([-5,16,-20,9])/160
    plus=np.full(4,.25)+delta;minus=np.full(4,.25)-delta
    fixtures=[]
    for s in (.3,.6,1.2):
        z=quotient_heat(2*s)[0];tail=rectangle_tail(2*s,32)
        for label,p in (('plus',plus),('minus',minus)):
            H=boundary.entropy_parts(p,dims)['extended']
            bound=2*s*float(p@C)+np.log(z+tail)
            assert H<=bound+1e-12
            fixtures.append(dict(state=label,s=s,entropy=H,
                reference_arc_entropy_upper=float(bound)))
    continuity=[]
    Sp=boundary.entropy_parts(plus,dims)['extended']
    for mix in (.2,.05,.01,.001):
        p=(1-mix)*plus+mix*minus
        distance=float(np.sum(abs(p-plus))/2)
        actual=abs(boundary.entropy_parts(p,dims)['extended']-Sp)
        # eta*cT=distance, E=2*cT*<C>; valid all-representation upper bound.
        Hb=boundary.entropy([distance,1-distance])
        bound=4*distance*float(plus@C)+2*distance*np.log(product_upper(2*distance))+Hb
        assert actual<bound
        continuity.append(dict(mixing=mix,regional_trace_distance=distance,
            actual_entropy_difference=actual,all_label_energy_continuity_upper=float(bound)))
    assert all(continuity[k+1]['all_label_energy_continuity_upper']<
               continuity[k]['all_label_energy_continuity_upper']
               for k in range(len(continuity)-1))
    return dict(original_electric_lower_coefficient=electric_min,
        original_scalar_kinetic_coefficient=scalar_coefficient,
        common_cT=cT,original_volume_weight=w,confining_u=u,
        slab_vertices=len(A),internal_links=inside,cut_links=cross,
        full_CAR_modes_per_region=32*len(A),regional_bounds=rows,
        original_loop_entropy_checks=fixtures,loop_continuity_checks=continuity,
        no_original_Gibbs_spectrum_or_entropy_numerically_substituted=True)


def run():
    deps=('research_note_598.md','research_note_603.md','research_note_617.md',
          'research_note_636.md','joint_thermal_gauss_source.py',
          'joint_full_spatial_metric.py','joint_region_energy_gluing.py',
          'joint_quotient_boundary_entropy.py',
          'joint_quotient_boundary_entropy_results.json')
    return dict(round=637,tests_run=2,failures=0,errors=0,
        quotient_reference_heat=quotient_heat_check(),original_region=region_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_original_H_controls_comparison_energy_analytically=True,
            matching_half_Casimir_cut_weights=True,
            full_fixed_graph_finite_energy_implies_finite_regional_entropy=True,
            full_Gauss_Gibbs_and_original_bounded_injection_records_covered=True,
            fixed_graph_energy_bounded_trace_approximation_controls_entropy=True,
            comparison_Gibbs_is_not_original_physical_equilibrium=True,
            no_uniform_continuum_or_area_law_or_Newton_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
