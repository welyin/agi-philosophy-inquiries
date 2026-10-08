"""Independent 1037 review of 1038; no imports from the author's science code."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'review_1037_results.json'
I = np.eye(2, dtype=complex)
P = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)
FROZEN = [HERE.parent/'research_note_1038.md'] + [HERE/x for x in
    ['moving_direction_transport.py', 'moving_direction_transport_results.json',
     'drafts/moving_direction_transport_derivation.md', 'selection_audit.md',
     'input_dependency_update.md', 'NEXT.md', 'review.md', 'verify_round1038.py']]


def sigma(v):
    return sum(a*b for a,b in zip(v,P))


def psqrt(a, inverse=False):
    w,v = np.linalg.eigh(a)
    assert w.min() > 0
    return (v*(1/np.sqrt(w) if inverse else np.sqrt(w)))@v.conj().T


def polar_wigner(p, n, xi, mass=1.):
    # Canonical rest boost = positive square root of the mass-shell matrix.
    p=np.asarray(p); e=np.sqrt(mass*mass+p@p)
    rest=psqrt((e*I+sigma(p))/mass)
    a=np.cosh(xi/2)*I+np.sinh(xi/2)*sigma(n)
    product=a@rest
    # Left polar decomposition; does not construct transformed momentum or use
    # the closed Wigner formula being reviewed.
    return psqrt(product@product.conj().T, inverse=True)@product


def closed_wigner(p,n,xi,mass=1.):
    p=np.asarray(p);v=p/(np.sqrt(mass*mass+p@p)+mass);t=np.tanh(xi/2)
    a=1+t*np.dot(n,v);b=t*np.cross(n,v)
    return (a*I+1j*sigma(b))/np.sqrt(a*a+b@b)


def matrix_checks():
    rng=np.random.default_rng(10371038)
    error=unitary_error=lipschitz_ratio=derivative_error=0.
    for k in range(30):
        p=rng.normal(size=3);n=rng.normal(size=3);n/=np.linalg.norm(n)
        mass=.7+.05*k;xi=-1.6+.1*k
        w=polar_wigner(p,n,xi,mass)
        error=max(error,float(np.linalg.norm(w-closed_wigner(p,n,xi,mass),2)))
        unitary_error=max(unitary_error,float(np.linalg.norm(w.conj().T@w-I,2)))
        dp=.02*rng.normal(size=3);w2=polar_wigner(p+dp,n,xi,mass)
        t=abs(np.tanh(xi/2));L=t/(2*mass*(1-t))
        if L>0:
            ratio=float(np.linalg.norm(w2-w,2)/(L*np.linalg.norm(dp)))
            lipschitz_ratio=max(lipschitz_ratio,ratio);assert ratio<=1+1e-10
        # Check sign and factor in the algebra's norm derivative.
        axis=k%3;direction=np.eye(3)[axis];h=2e-5
        wp=polar_wigner(p,direction,h,mass);wm=polar_wigner(p,direction,-h,mass)
        d=(wp.conj().T@P[axis]@wp-wm.conj().T@P[axis]@wm)/(2*h)
        v=p/(np.sqrt(mass*mass+p@p)+mass)
        derivative_error=max(derivative_error,float(np.linalg.norm(d-sigma(v-v[axis]*direction),2)))
    assert error<1e-12 and unitary_error<1e-12 and derivative_error<1e-8
    cplus=polar_wigner([0,0,4/3],[1,0,0],math.log(3))
    cminus=polar_wigner([0,0,-4/3],[1,0,0],math.log(3))
    assert np.linalg.norm(cplus-(4*I-1j*P[1])/np.sqrt(17))<1e-13
    assert np.linalg.norm(cminus-(4*I+1j*P[1])/np.sqrt(17))<1e-13
    return dict(cases=30,polar_vs_closed_residual=error,unitarity_residual=unitary_error,
                largest_lipschitz_ratio=lipschitz_ratio,derivative_residual=derivative_error,
                wigner_center_signs_verified=True)


def packet_check():
    xr,wr=np.polynomial.legendre.leggauss(14)
    xc,wc=np.polynomial.legendre.leggauss(8)
    width=.01;normalization=0.;states=[np.zeros((2,2),complex) for _ in range(2)]
    up=np.diag([1.,0.]);down=np.diag([0.,1.]);axis=np.array([1.,0.,0.])
    for r,dr in zip((xr+1)/2,wr/2):
        for z,dz in zip(xc,wc):
            for j in range(16):
                phi=2*np.pi*j/16
                p=np.array([0.,0.,4/3])+width*r*np.array([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
                weight=width**3*dr*dz*(2*np.pi/16)*r*r*np.exp(-2/(1-r*r))/(2*np.sqrt(1+p@p))
                wp=polar_wigner(p,axis,math.log(3));wm=polar_wigner(-p,axis,math.log(3))
                states[0]+=weight*(wp@up@wp.conj().T+wm@down@wm.conj().T)/2
                states[1]+=weight*(wp@down@wp.conj().T+wm@up@wm.conj().T)/2
                normalization+=weight
    states=[s/normalization for s in states]
    probs=[float(np.trace(s@(I+P[0])/2).real) for s in states]
    lower=F(8,17)-2*F(1,2)*F(1,100)
    assert lower==F(783,1700) and lower/2==F(783,3400)
    assert probs[0]-probs[1]>float(lower)
    saved=json.loads((HERE/'moving_direction_transport_results.json').read_text('utf8'))
    savedgap=saved['normal_packet']['fine']['probability_gap']
    assert abs(probs[0]-probs[1]-savedgap)<1e-8
    return dict(points_per_packet=14*8*16, invariant_measure_bump_norm_squared=normalization,
                probabilities=probs,gap=probs[0]-probs[1],
                difference_from_author_gap=abs(probs[0]-probs[1]-savedgap),
                rigorous_gap=str(lower), predictor_error_lower=str(lower/2),
                quadrature_is_not_an_interval_certificate=True)


def trace_spin(a,ref=3):
    return np.trace(a.reshape(2,ref,2,ref),axis1=0,axis2=2)


def reference_tail_check():
    rng=np.random.default_rng(1038037);ref=3
    centers=np.array([[.2,-.1,.7],[-.3,.2,-.8]])
    momenta=np.array([centers[0]+[.004,0,.001],centers[0]+[-.003,.002,0],
                      centers[1]+[0,-.002,.003],centers[1]+[.002,.001,-.004],
                      [3.,-2.,4.]])
    eta=.03;weights=np.array([(1-eta)/4]*4+[eta])
    conditional=rng.normal(size=(5,2*ref))+1j*rng.normal(size=(5,2*ref))
    conditional/=np.linalg.norm(conditional,axis=1)[:,None]
    psi=(np.sqrt(weights)[:,None]*conditional).ravel()
    rho=np.outer(psi,psi.conj());blocks=rho.reshape(5,6,5,6)
    old_reference=sum(trace_spin(blocks[k,:,k,:]) for k in range(5))
    radius=max(np.linalg.norm(momenta[k]-centers[k//2]) for k in range(4))
    rows=[]
    for n,xi,mass in [(np.array([1.,0.,0.]),.9,.7),
                      (np.array([1.,2.,-1.])/np.sqrt(6),-1.3,1.2)]:
        U=np.zeros((30,30),complex)
        for k,p in enumerate(momenta):
            U[6*k:6*k+6,6*k:6*k+6]=np.kron(polar_wigner(p,n,xi,mass),np.eye(ref))
        out=U@rho@U.conj().T
        exact=sum(out.reshape(5,6,5,6)[k,:,k,:] for k in range(5))
        approx=np.zeros((6,6),complex)
        for k,p in enumerate(centers):
            f=blocks[2*k,:,2*k,:]+blocks[2*k+1,:,2*k+1,:]
            u=np.kron(polar_wigner(p,n,xi,mass),np.eye(ref))
            approx+=u@f@u.conj().T
        tail=blocks[4,:,4,:]
        approx+=np.kron(I/2,trace_spin(tail))
        error=float(np.abs(np.linalg.eigvalsh(exact-approx)).sum()/2)
        t=abs(np.tanh(xi/2));budget=t/(2*mass*(1-t))*radius+eta
        referr=float(np.linalg.norm(trace_spin(approx)-old_reference,2))
        assert error<=budget and referr<1e-13
        assert abs(np.trace(approx)-1)<1e-13 and np.linalg.eigvalsh(approx).min()>-1e-13
        rows.append(dict(xi=xi,mass=mass,trace_distance=error,analytic_budget=float(budget),
                         tail_weight=eta,reference_marginal_residual=referr))
    return dict(reference_dimension=3,total_dimension=30,input_is_pure=True,
                momentum_coherences_retained_in_input=True,bin_radius=float(radius),cases=rows)


def run():
    before={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in FROZEN}
    result=dict(round_reviewed=1038,reviewer='1037 research agent',date='2026-10-08',
                representation=matrix_checks(),normal_packet=packet_check(),
                spectator_reference_with_nonzero_tail=reference_tail_check(),
                frozen_author_files_sha256=before,added_scientific_calibration_groups=0,
                physical_detector_or_full_parent_certified=False,passed=True)
    after={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in FROZEN}
    assert before==after
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:
        previous=json.loads(OUT.read_text('utf8'))
        assert previous==result
    print(json.dumps({k:v for k,v in result.items() if k!='frozen_author_files_sha256'},ensure_ascii=False))
