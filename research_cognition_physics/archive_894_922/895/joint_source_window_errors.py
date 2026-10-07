"""895: one finite-frequency window for the original mass/stress source pair.
Bounds concern one-loop quadratic source kernels on the declared flat reference.
They do not approximate the full interacting record process or prove feedback.
"""
from pathlib import Path
from functools import lru_cache
from math import comb, factorial
import json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
with ResearchRuntime(Layout()).installed():
    import joint_continuum_source_spectrum as scalar
    import joint_tensor_stress_spectrum as tensor
TARGET=HERE/'joint_source_window_errors_results.json'

@lru_cache(None)
def quad(n):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)/2,w/2

def rho(omega,p,m,d):
    return scalar.density(omega,m,d) if p==2 else tensor.spin2(omega,m,d)

def spectral_tail(z,cut,p,m,d,n):
    x,w=quad(n);om=cut/x
    return z**(p+2)*np.dot(w, rho(om,p,m,d)/(np.pi*om**(p+1)*(om**2-z*z))*cut/x**2)

def full_remainder(z,p,m,d,n):
    return scalar.finite_subtracted(z,m,d,n) if p==2 else tensor.tensor_dispersion(z,m,d,n)

def finite_remainder(z,cut,p,m,d,n):
    # Independent integration of each retained threshold interval.
    u,w=quad(n);ans=0j
    for mi,di in zip(m,d):
        if 2*mi>=cut:continue
        # x=2m/omega, integrate x from 2m/cut to 1.
        lo=2*mi/cut;x=lo+(1-lo)*u
        base=x*(1-x*x)**1.5/(1-(z*x/(2*mi))**2)
        if p==2:ans+=z**4*di/(16*np.pi**2)*(1-lo)*np.dot(w,base)
        else:ans+=z**6*di/(960*np.pi**2*mi*mi)*(1-lo)*np.dot(w,base*(3+2*x*x))
    return ans

def pulse_coefficients():
    c=np.zeros(17)
    for j in range(9):c[8+j]=4.**8*(-1)**j*comb(8,j)
    return c

def bernstein_derivative_bound(r):
    b=np.zeros(17);b[8]=4.**8/comb(16,8)
    for j in range(r):b=(16-j)*np.diff(b)
    return float(np.max(abs(b)))

def pulse_fourier(omega,T):
    # Compact pulse phi(t/T)=(4x(1-x))^8, C^7 globally.
    out=np.zeros_like(np.asarray(omega),dtype=complex);s=np.asarray(omega)*T
    low=abs(s)<40
    if np.any(low):
        x,w=quad(128);shape=(4*x*(1-x))**8
        out[low]=T*(np.exp(1j*np.outer(s[low],x))@(w*shape))
    if np.any(~low):
        z=1j*s[~low];v=0j;coeff=pulse_coefficients()
        for j in range(8,17):
            at0=factorial(j)*coeff[j]
            # symmetry phi(1-x)=phi(x), so endpoint derivatives differ by sign
            at1=(-1)**j*at0
            v=v+(-1)**j*(at1*np.exp(z)-at0)/z**(j+1)
        out[~low]=T*v
    return out

def noise_tail(cut,p,m,d,T,n):
    x,w=quad(n);om=cut/x;ft=pulse_fourier(om,T)
    return float(np.dot(w,rho(om,p,m,d)*abs(ft)**2/(2*np.pi)*cut/x**2))

def run():
    _,m,d,phi,_=scalar.original_data()
    constants={2:float(np.dot(d,m*m)/(4*np.pi**2)),4:float(d.sum()/(80*np.pi**2))}
    b=float(1/np.sqrt(6)/np.tanh(scalar.QSTAR/np.sqrt(6)));v=np.array([b,1.])
    rows=[]
    for p in (2,4):
        cp=constants[p]
        for cut in (1.,2.,4.):
            for z in (.07+.03j,.2+.1j):
                tail=spectral_tail(z,cut,p,m,d,256)
                tail2=spectral_tail(z,cut,p,m,d,512)
                bound=cp*abs(z)**(p+2)/(2*cut**2*(1-abs(z)**2/cut**2))
                full=full_remainder(z,p,m,d,384)
                finite=finite_remainder(z,cut,p,m,d,384)
                identity=abs(full-finite-tail)
                assert abs(tail)<=bound*(1+1e-12)
                assert identity<2e-12 and abs(tail-tail2)<2e-14
                rows.append(dict(channel_power=p,cutoff=cut,z_real=z.real,z_imag=z.imag,
                    source_kernel_tail_norm=float(abs(tail)),source_kernel_tail_bound=bound,
                    independent_full_minus_finite_error=float(identity),tail_quadrature_change=float(abs(tail-tail2))))
    T=80.;eps=.001;r=4;window=[]
    for cut in (1.,2.,4.):
        entry=dict(cutoff=cut,pulse_duration=T,source_amplitude=eps)
        for p in (2,4):
            cp=constants[p]
            # Rigorous Bernstein convex hull bounds on derivatives; no fit.
            fbound=bernstein_derivative_bound(p+2)/T**(p+2)
            dfbound=bernstein_derivative_bound(p+3)/T**(p+2)
            rtime=cp/(2*cut**2)*(fbound+dfbound)
            deriv_l1=bernstein_derivative_bound(r)*T**(1-r)
            nbound=cp*deriv_l1**2/(2*(2*r-p-1)*cut**(2*r-p-1))
            n1=noise_tail(cut,p,m,d,T,256);n2=noise_tail(cut,p,m,d,T,512)
            assert n2>=0 and n2<=nbound*(1+1e-12)
            entry[str(p)]=dict(causal_source_error_bound=eps*rtime,
                noise_quadratic_error_bound=eps**2*nbound,
                numerical_noise_tail=eps**2*n2,
                noise_tail_quadrature_change=eps**2*abs(n2-n1))
        entry['joint_seven_source_response_bound']=max(np.dot(v,v)*entry['2']['causal_source_error_bound'],entry['4']['causal_source_error_bound'])
        window.append(entry)
    # Exact retained covariance (positive measure) preserves the same matrix
    # directions; no independent fit of scalar, conformal and shear channels.
    omega=np.array([1.3]);rm=float(rho(omega,2,m,d)[0]);rt=float(rho(omega,4,m,d)[0])
    joint=np.zeros((7,7));joint[:2,:2]=rm*np.outer(v,v);joint[2:,2:]=rt*np.eye(5)
    null=np.r_[1.,-b,np.zeros(5)]
    assert abs(null@joint@null)<1e-12 and np.linalg.eigvalsh(joint).min()>-1e-12
    # Endpoint Fourier formula cross-check in a moderate oscillatory window.
    x,w=quad(512);om=np.array([.51,.8,1.2]);ft=pulse_fourier(om,T)
    direct=T*np.exp(1j*np.outer(om*T,x))@(w*(4*x*(1-x))**8)
    ferr=float(np.max(abs(ft-direct)));assert ferr<2e-12
    return dict(round=895,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3680,
        argument_scope='Original630/631 one-generation vacuum, flat constant background, zero spatial momentum, one-fermion-loop quadratic sources; common finite window preserves positive noise and causal subtracted response with explicit tail bounds.',
        original_Weyl_components=16,original_masses=np.sort(m).tolist(),
        positive_spectral_envelopes={str(k):v for k,v in constants.items()},
        mass_geometry_direction=[b,1.],complex_response_checks=rows,
        compact_pulse_window_checks=window,Fourier_endpoint_crosscheck_error=ferr,
        original_joint_null_residual=float(abs(null@joint@null)),
        finite_local_matching_coefficients_unchanged=True,
        same_positive_measure_used_for_noise_and_response=True,
        microscopic_minimum_scale_assumed=False,horizontal_transport_added=False,
        actual_feedback_error_bound_requires_stability_estimate=True,
        full_interacting_probability_remainder_proved=False,
        full_general_background_Ward_reproved=False,full_goal_completed=False)

if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
