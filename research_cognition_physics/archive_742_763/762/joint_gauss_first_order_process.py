"""762: first-order Gauss packet, matrix transport and actual source histories.

Three coefficient checks; no full graph or continuum dynamics simulation.
"""
import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import numpy as np
import joint_loop_scale_transport as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_first_order_process_results.json'


def fall(n,k):
    return factorial(n)//factorial(n-k)


def add(a,b,scale=1):
    out={k:v.copy() for k,v in a.items()}
    for k,v in b.items():
        out[k]=out.get(k,np.zeros_like(v))+scale*v
    return {k:v for k,v in out.items() if np.max(abs(v))>1e-13}


def scale(a,x):
    return {k:x*v for k,v in a.items() if np.max(abs(x*v))>1e-13}


def multiply(a,b):
    out={}
    for (q,p,e),A in a.items():
        for (r,s,f),B in b.items():
            key=(q+r,p+s,e+f)
            value=A@B
            out[key]=out.get(key,np.zeros_like(value))+value
    return out


def diff(a,axis):
    out={}
    for key,A in a.items():
        n=key[axis]
        if n:
            k=list(key);k[axis]-=1
            out[tuple(k)]=n*A
    return out


def poisson(a,b):
    # P(a,b)=a_p b_q-a_q b_p; matrix order matters.
    return add(multiply(diff(a,1),diff(b,0)),
               multiply(diff(a,0),diff(b,1)),-1)


def epsilon_shift(a,k):
    return {(q,p,e+k):A for (q,p,e),A in a.items()}


def weyl_to_ordered(a):
    # A q^r p^s epsilon^e -> sum A (-i)^s epsilon^(e+s)
    # binom(s,k) (r)_k / 2^k q^(r-k) d_q^(s-k).
    out={}
    for (q,p,e),A in a.items():
        for k in range(min(q,p)+1):
            key=(q-k,p-k,e+p)
            value=(-1j)**p*comb(p,k)*fall(q,k)/2**k*A
            out[key]=out.get(key,np.zeros_like(value))+value
    return out


def compose_ordered(a,b):
    out={}
    for (q,d,e),A in a.items():
        for (r,j,f),B in b.items():
            for k in range(min(d,r)+1):
                key=(q+r-k,d+j-k,e+f)
                value=comb(d,k)*fall(r,k)*(A@B)
                out[key]=out.get(key,np.zeros_like(value))+value
    return out


def small_order_residual(a,order=2):
    return max([float(np.max(abs(A))) for (q,d,e),A in a.items()
                if e-d<order]+[0.])


def gauss_packet_check():
    # One compact group direction calibrated with original hypercharge values:
    # vacuum Q=0 and a paired e_R state Q=-12. This is not a full Gauss graph.
    Q=np.diag([0.,-12.])
    A=np.array([[.8,.1-.2j],[.1+.2j,-.4]])
    N=np.array([[.2,.3+.07j],[.3-.07j,.9]])
    C=np.array([[.4,.11j],[-.11j,.7]])
    E=np.array([[1.,.1],[.1,.6]])
    F=np.array([[.3,.2-.15j],[.2+.15j,-.1]])
    c=.8; ell=.7
    theta=2*np.pi*np.arange(128)/128
    R=np.zeros((128,2,2),complex)
    R[:,0,0]=1
    R[:,1,1]=np.exp(-12j*theta)
    modes=np.fft.fftfreq(128,d=1/128)
    def derivative(values):
        return np.fft.ifft(1j*modes[:,None,None]*np.fft.fft(values,axis=0),axis=0)
    Rprime=derivative(R)
    Ftheta=R@F@R.conj().transpose(0,2,1)
    # Direct Weyl symmetrization of F(theta)*p_theta.
    direct=-1j*(Ftheta@Rprime+.5*derivative(Ftheta)@R)
    orbit=np.mean(R.conj().transpose(0,2,1)@direct,axis=0)
    expected_orbit=(Q@F+F@Q)/2
    orbit_error=float(np.max(abs(orbit-expected_orbit)))
    assert orbit_error<1e-12
    nonsymmetric_error=float(np.linalg.norm(F@Q-expected_orbit))
    assert nonsymmetric_error>.5

    nodes,weights=np.polynomial.hermite.hermgauss(80)
    weights=weights/np.sqrt(np.pi)
    L1=c*ell*N+c*C+E/(4*c)+expected_orbit
    rows=[]
    for eps in (.04,.02,.01,.005):
        n=np.sqrt(2*eps*c)*nodes
        density=weights*np.exp(ell*n)
        density/=density.sum()
        mean=float(density@n); second=float(density@(n*n))
        logarithmic_derivative=-n/(2*eps*c)+ell/2
        centered_p_square=eps**2*float(density@(logarithmic_derivative**2))
        actual=A+mean*N+second*C+centered_p_square*E+eps*orbit
        coefficient=(actual-A)/eps
        error=float(np.max(abs(coefficient-L1)))
        exact_remaining=eps*c*c*ell*ell*C
        assert np.max(abs(coefficient-L1-exact_remaining))<1e-12
        rows.append(dict(epsilon=eps,first_coefficient_error=error))
    return dict(group_jet_error=orbit_error,
                wrong_unsymmetrized_charge_error=nonsymmetric_error,
                density_gradient_and_momentum_variance_retained=True,rows=rows,
                scope='Local Gaussian times a compact U(1) orbit coefficient calibration, '
                      'using actual hypercharge values. Not original graph propagation.')


def transport_algebra_check(data):
    phi,x,F,B,M,vp,vm,psi=data
    I=np.eye(len(B),dtype=complex)
    h={(0,2,0):.5*I,(2,2,0):I/12,(4,0,0):3*I/8,(2,0,0):I/9}
    b={(0,0,0):B,(1,0,0):M}
    H=add(h,epsilon_shift(b,1))
    c={(0,0,0):M,(0,1,0):B,(1,0,0):M@M,
       (1,1,0):(M@B+B@M)/2,(0,2,0):B@B}
    leading=add(poisson(h,c),scale(add(multiply(b,c),multiply(c,b),-1),1j))
    correction=scale(add(poisson(b,c),poisson(c,b),-1),.5)
    theory=add(leading,epsilon_shift(correction,1))
    Hop=weyl_to_ordered(H);Cop=weyl_to_ordered(c)
    direct=epsilon_shift(scale(add(compose_ordered(Hop,Cop),
                                    compose_ordered(Cop,Hop),-1),1j),-1)
    error=small_order_residual(add(direct,weyl_to_ordered(theory),-1))
    assert error<2e-11
    wrong=add(leading,epsilon_shift(poisson(b,c),1))
    wrong_error=small_order_residual(add(direct,weyl_to_ordered(wrong),-1))
    assert wrong_error>1e-4

    # Original positive L_+(s(x)) through its second x5 jet.
    u=1+float(x@x)/6
    A=1+float(x[:4]@x[:4])/6
    ds=np.sqrt(2)*A/u**1.5
    dds=-np.sqrt(2)*A*x[4]/(2*u**2.5)
    e=.5+np.sin(phi[4])/4
    de=np.cos(phi[4])*ds/4
    dde=(np.cos(phi[4])*dds-np.sin(phi[4])*ds*ds)/4
    l0=np.sqrt(e);l1=de/(2*l0)
    l2=dde/(4*l0)-de*de/(8*l0**3)
    L={(0,0,0):l0*I,(1,0,0):l1*I,(2,0,0):l2*I}
    Lop=weyl_to_ordered(L)
    sandwich=compose_ordered(Lop,compose_ordered(Cop,Lop))
    leading_sandwich=weyl_to_ordered(multiply(multiply(L,c),L))
    instrument_error=small_order_residual(add(sandwich,leading_sandwich,-1))
    assert instrument_error<2e-11
    return dict(original_neutral_Fock_dimension=len(B),
                noncommuting_original_mass_and_source=float(np.linalg.norm(B@M-M@B)),
                direct_differential_operator_first_order_error=error,
                wrong_one_sided_matrix_Poisson_bracket_error=wrong_error,
                original_real_Kraus_first_order_star_term_error=instrument_error,
                scope='Independent differential-operator calculation of matrix-symbol identities. '
                      'Original neutral mass and original Kraus jets; scalar polynomial h is '
                      'an algebra calibration, not a solved original trajectory.')


def original_source_conditioning(data):
    phi,x,F,B,M,vp,vm,psi=data
    # Two invariant slice coordinates (Higgs radial mass coordinate, sterile).
    q0=np.array([np.linalg.norm(x[:4]),x[4]])
    covariance=np.array([[.6,.17],[.17,1.1]])
    root=np.linalg.cholesky(covariance)
    nodes,w=np.polynomial.hermite.hermgauss(48)
    grid=np.stack(np.meshgrid(nodes,nodes,indexing='ij'),axis=-1).reshape(-1,2)
    weights=np.outer(w,w).ravel()/np.pi
    ys=float(abs(old.old.matter.Y['s']))
    mu=.4*ys  # actual mixture .7*v+ + .3*v-
    assert abs(.7*np.vdot(vp,M@vp).real+.3*np.vdot(vm,M@vm).real-mu)<1e-13
    assert abs(.7*np.vdot(vp,M@M@vp).real+.3*np.vdot(vm,M@M@vm).real-ys**2)<1e-13

    def value(q):
        u=1+np.sum(q*q,axis=-1)/6
        s=np.sqrt(2)*q[...,1]/np.sqrt(u)
        return .5+np.sin(s)/4
    u=1+q0@q0/6
    k=np.array([0.,1.])
    gradient_s=np.sqrt(2)*(k/u**.5-q0[1]*q0/(6*u**1.5))
    hessian_s=np.sqrt(2)*(
        -(np.outer(k,q0)+np.outer(q0,k)+q0[1]*np.eye(2))/(6*u**1.5)
        +q0[1]*np.outer(q0,q0)/(12*u**2.5))
    s0=float(phi[4])
    grad=np.cos(s0)*gradient_s/4
    hess=(np.cos(s0)*hessian_s-np.sin(s0)*np.outer(gradient_s,gradient_s))/4
    e0=float(value(q0))
    p1=float(np.sum(covariance*hess))/2
    delta=float((covariance@grad)[1])/e0
    mean0=mu*q0[1];mean1=mu*delta
    variance0=(ys**2-mu**2)*q0[1]**2
    variance1=ys**2*covariance[1,1]+2*(ys**2-mu**2)*q0[1]*delta
    # Finite-difference cross-check of analytic derivatives; not used in prediction.
    d=1e-4
    numerical=np.empty((2,2))
    for i in range(2):
        for j in range(2):
            ei=np.eye(2)[i]*d;ej=np.eye(2)[j]*d
            numerical[i,j]=(value(q0+ei+ej)-value(q0+ei-ej)
                            -value(q0-ei+ej)+value(q0-ei-ej))/(4*d*d)
    derivative_error=float(np.max(abs(hess-numerical)))
    assert derivative_error<2e-8

    import sys
    sys.path.insert(0,str(HERE/'round762_drafts'))
    from resource_scale_entry import smooth_cutoff
    rows=[]
    for eps in (1/4096,1/8192,1/16384,1/32768):
        shift=np.sqrt(2*eps)*(grid@root.T)
        density=weights*np.prod(smooth_cutoff(shift/.2)**2,axis=1)
        density/=density.sum()
        q=q0+shift
        assert np.min(q[:,0])>0
        effects=value(q)
        p=float(density@effects)
        post=density*effects/p
        mean=float(post@(mu*q[:,1]))
        second=float(post@(ys**2*q[:,1]**2))
        variance=second-mean**2
        errors=[abs((p-e0)/eps-p1),abs((mean-mean0)/eps-mean1),
                abs((variance-variance0)/eps-variance1)]
        assert max(errors)<1e-4
        rows.append(dict(epsilon=eps,probability=p,conditional_mean=mean,
                         conditional_source_variance=variance,
                         first_coefficient_errors=errors))
    assert max(rows[-1]['first_coefficient_errors'])<max(rows[0]['first_coefficient_errors'])
    assert abs(mean1)>1e-3 and variance0>0
    return dict(original_s=s0,covariance=covariance.tolist(),
                probability_first_coefficient=p1,
                conditional_Majorana_energy_first_coefficient=mean1,
                conditional_variance_leading=variance0,
                conditional_variance_first_coefficient=variance1,
                derivative_crosscheck_error=derivative_error,rows=rows,
                scope='Original invariant h/s slice, real field effect and true mixed sterile '
                      'CAR source. Half-density Gaussian preparation, compact cutoff. Initial '
                      'coefficient quadrature only, not full graph dynamics or Einstein fluctuations.')


def run():
    data=old.original_data()
    result=dict(round=762,tests_run=3,failures=0,errors=0,
                gauss_initial_jet=gauss_packet_check(),
                transport_and_instrument=transport_algebra_check(data),
                original_conditional_source=original_source_conditioning(data))
    deps=['research_note_525.md','research_note_574.md','research_note_758.md',
          'research_note_761.md','joint_loop_scale_transport.py',
          'round762_drafts/resource_scale_entry.py',
          'round762_drafts/research_note_762_working.md']
    result['dependency_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps}
    result['scope']='Conditional fixed-graph first-order observable and finite-history expansion '
    result['scope']+='for the explicit 761 scaling, with original Gauss packet data and all B '
    result['scope']+='retained analytically. No first-order trace-norm state theorem, uniform '
    result['scope']+='spatial refinement, epsilon=1 accuracy, continuum renormalization or quantum GR.'
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))

