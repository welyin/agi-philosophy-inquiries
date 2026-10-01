"""646: the actual612 free overlap spectrum in a common continuum window.

The same positive spectral measure carries a correlation and two energy
insertions. Nonzero-time convergence does not include coincident second
moments. Fixed physical spatial torus, no interacting Gauss/SM/GR claim.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_overlap_time_spectrum as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_overlap_continuum_sources_results.json'
PBOUND=2.
TAU=np.pi/16
TMAX=np.pi/2


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def data(size, p):
    a=2*np.pi/size
    m=old.parameters(a*np.asarray(p,float))
    m.update(a=a,p=np.asarray(p,float),B=np.sqrt(1+m['s2']),
             delta=m['b']/(2*np.sqrt(1+m['s2'])))
    assert a*PBOUND<=.5
    return m


def cut_moment(m,t,power,order=240):
    # E=E1+u^2; numerical tail bound is explicit, not a physical cutoff.
    v,w=nodes(order);u=(v+1)*np.sqrt(90.)/2
    shift=u*u;E=m['E1']+shift
    radical=np.sqrt(4*m['b']*np.sinh(m['E1']+shift/2)*np.sinh(shift/2))
    density=np.sinh(E)*radical/(2*np.pi*(np.sinh(E)**2-m['s2']))
    weight=w*np.sqrt(90.)/2*2*u*density
    result=float(np.dot(weight,(E/m['a'])**power*np.exp(-t*E/m['a'])))
    # For a|p|<=1/2: b<=1/8, s^2<=1/4, rho(E)<=sqrt(b)exp(-E/2).
    top=m['E1']+90;lam=.5+t/m['a']
    poly=top**power/lam
    if power>=1:poly+=power*top**(power-1)/lam**2
    if power>=2:poly+=power*(power-1)/lam**3
    tail=float(np.sqrt(m['b'])/m['a']**power*np.exp(-lam*top)*poly)
    return result,tail


def moment(m,t,power,order=240):
    e=m['E0']/m['a'];cut,tail=cut_moment(m,t,power,order)
    return float(m['w']*e**power*np.exp(-t*e)+cut),tail


def window_bound(a,power):
    L=np.log(2/(a*a*PBOUND**2))/a
    assert TAU*L>=power
    cut=a*a*PBOUND**2/4*L**power*np.exp(-TAU*L)
    deriv=(power*PBOUND**(power-1) if power else 0)+TMAX*PBOUND**power
    pole=a*a*PBOUND**(power+2)/4+deriv*a*a*PBOUND**3/6
    return float(pole+cut),float(cut),float(L)


def normalization_check():
    rows=[];identity_error=0.
    for size in (32,64,128):
        for p in ((1,0,0),(1,1,1)):
            m=data(size,p)
            mass,tail=cut_moment(m,0.,0,order=320)
            err=abs(mass-m['delta']);assert err<2e-14
            assert abs(m['w']+mass-.5)<2e-14 and tail<1e-19
            # Independent transformed density rho dE = rational(y) dy.
            for y in (.15,.7,2.1,7.):
                E=np.arccosh((m['B']**2+m['b']**2+y*y)/(2*m['b']))
                actual=old.density(E,m)*y/(m['b']*np.sinh(E))
                expected=2*m['b']/np.pi*y*y/((y*y+(m['B']-m['b'])**2)
                                                   *(y*y+(m['B']+m['b'])**2))
                identity_error=max(identity_error,float(abs(actual-expected)))
            # Original lattice Fourier data only at actual positive integer times.
            omega=2*np.pi*(np.arange(32768)+.5)/32768-np.pi
            ff=old.scalar_fourier(omega,m)
            n=size//32
            actual=np.mean(np.exp(1j*omega*n)*ff)
            expected,_=moment(m,n*m['a'],0)
            errf=float(abs(actual-expected));assert errf<2e-12
            rows.append(dict(spatial_sites_per_side=size,physical_momentum=list(p),
                a=m['a'],pole_physical_energy=m['E0']/m['a'],
                cut_physical_threshold=m['E1']/m['a'],cut_weight=m['delta'],
                total_weight_error=err,tail_bound=tail,integer_time=n,Fourier_error=errf))
    assert identity_error<2e-14
    return dict(rows=rows,rational_density_identity_error=identity_error,
                physical_torus_side=2*np.pi,original_kernel_at_allowed_momenta=True,
                spectral_right_limit_not_the_raw_lattice_equal_time_value=True)


def common_window_check():
    rows=[];max_quad=0.;max_source=0.
    for size in (32,64,128,256):
        a=2*np.pi/size;errs=np.zeros(3);ratios=np.zeros(3)
        for p in ((1,0,0),(1,1,0),(1,1,1)):
            m=data(size,p);r=float(np.linalg.norm(p))
            for t in (TAU,2*TAU,TMAX):
                assert abs(t/a-round(t/a))<1e-12
                for j in range(3):
                    val,tail=moment(m,t,j);fine,_=moment(m,t,j,order=320)
                    max_quad=max(max_quad,abs(val-fine))
                    target=.5*r**j*np.exp(-t*r)
                    errs[j]=max(errs[j],abs(val-target)+tail)
                    bound,cutbound,_=window_bound(a,j)
                    cut,ctail=cut_moment(m,t,j)
                    assert cut+ctail<=cutbound*(1+2e-12)+1e-300
                    assert abs(val-target)+tail<=bound
                    ratios[j]=max(ratios[j],(cut+ctail)/cutbound if cutbound>0 else 0)
            # One shared constant lapse on the reconstructed cyclic Hamiltonian.
            t=2*TAU;h=1e-4
            c0=moment(m,t,0)[0];cp=moment(m,t*(1+h),0)[0];cm=moment(m,t*(1-h),0)[0]
            fd1=(cp-cm)/(2*h);fd2=(cp-2*c0+cm)/h**2
            exact1=-t*moment(m,t,1)[0];exact2=t*t*moment(m,t,2)[0]
            max_source=max(max_source,abs(fd1-exact1),abs(fd2-exact2))
        rows.append(dict(spatial_sites_per_side=size,a=a,
            max_window_errors_orders012=errs.tolist(),
            rigorous_uniform_bounds_orders012=[window_bound(a,j)[0] for j in range(3)],
            cut_to_uniform_bound_ratios=ratios.tolist()))
    assert max_quad<3e-13 and max_source<5e-8
    return dict(rows=rows,physical_momentum_bound=PBOUND,physical_time_window=[TAU,TMAX],
                quadrature_difference=max_quad,constant_lapse_difference_error=max_source,
                full_tensor_metric_source_not_claimed=True,
                low_energy_matching_does_not_impose_a_new_hard_physical_cutoff=True)


def coincident_source_check():
    rows=[];last_lower=0.
    for size in (32,64,128,256,512,1024,2048):
        m=data(size,(1,0,0));values=[];tails=[];error=0.
        for j in range(3):
            value,tail=cut_moment(m,0.,j,order=320)
            second,_=cut_moment(m,0.,j,order=400)
            error=max(error,abs(value-second));values.append(value);tails.append(tail)
        lower=m['delta']*(m['E1']/m['a'])**2
        assert values[2]+tails[2]>=lower and lower>last_lower
        assert error<3e-11 and max(tails)<2e-15
        last_lower=lower
        rows.append(dict(spatial_sites_per_side=size,a=m['a'],
            cut_moments_orders012=values,coincident_second_moment_lower_bound=lower,
            quadrature_difference=error,tail_bounds_orders012=tails,
            total_right_limit_second_moment=values[2]+m['w']*(m['E0']/m['a'])**2,
            limiting_pole_second_moment=.5))
    assert rows[-1]['cut_moments_orders012'][1]<rows[0]['cut_moments_orders012'][1]
    return dict(rows=rows,
                second_moment_divergence_proved_by_threshold_lower_bound=True,
                raw_lattice_contact_prescription_not_replaced=True,
                compact_resolvent_exact_fixed_a_obstruction_from612_unchanged=True,
                no_full_Gauss_thermal_trace_or_stress_renormalization_solution=True)


def run():
    deps=('research_note_612.md','research_note_613.md','research_note_614.md',
          'research_note_615.md','research_note_645.md','joint_overlap_time_spectrum.py')
    return dict(round=646,tests_run=3,failures=0,errors=0,
        actual_spectral_dictionary=normalization_check(),
        common_correlation_energy_window=common_window_check(),
        coincident_source_boundary=coincident_source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(free_massless_frozen_original_overlap_slice_only=True,
            fixed_physical_3torus_and_refinement_are_inputs=True,
            uniform_correlation_and_two_energy_insertions_on_nonzero_time_window=True,
            coincident_second_energy_moment_not_convergent=True,
            no_interacting_chiral_Gauss_state_or_region_continuum_map=True,
            no_quantum_gravity_SM_or_GR_generation_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
