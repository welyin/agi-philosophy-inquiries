"""All-degree double-degeneracy bound and physical scalar response of two color paths."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ARCH=HERE.parent


def read(rel):
    return json.loads((ARCH/rel).read_bytes())


def absolute_bound():
    data=read('1397/inputs.json')['high_scale']
    def high(row):
        y,err=map(lambda x:Fraction(str(x)),row)
        return (y+2*err)**2
    a,b=(high(data[k][2]) for k in ('up','down'))
    da,db=(high(data[k][1]) for k in ('up','down'))
    ratio=4*max(a,b)/3
    assert ratio<1
    harmonic=sum(Fraction(1,j) for j in range(1,7))
    bound=18*da*db*a*a*b*b*(harmonic+ratio/(7*(1-ratio)))
    return dict(a=float(a),b=float(b),delta_a_upper=float(da),delta_b_upper=float(db),
                tail_ratio=float(ratio),strict_rational_upper=str(bound),
                outward_upper_rad=math.ceil(bound*10**26)/10**26)


def word(A,B,w):
    v=np.eye(3,dtype=complex)
    for index in w:
        v=v@(A if index=='A' else B)
    return v


def verify_trace_identities(A,B):
    eig,U=np.linalg.eigh(A)
    spec,V=np.linalg.eigh(B)
    CKM=U.conj().T@V
    J=(CKM[0,0]*CKM[1,1]*CKM[0,1].conjugate()*CKM[1,0].conjugate()).imag
    def vandermonde(x):
        return (x[1]-x[0])*(x[2]-x[0])*(x[2]-x[1])
    theory=abs(J*vandermonde(eig)*vandermonde(spec))
    # Compute Im trace in the basis splitting each matrix's light pair,
    # using the commutator determinant. This avoids the huge real trace term.
    comm=A@B-B@A
    measured=abs(np.linalg.det(comm).imag)/2
    np.testing.assert_allclose(measured/theory,1,rtol=3.e-10)
    # The mixed finite difference identity is independent of number of trace
    # factors. Test on nonhierarchical examples as well, to resolve the scale.
    sampler=np.random.default_rng(1686)
    V=np.linalg.qr(sampler.normal(size=(3,3))+1j*sampler.normal(size=(3,3)))[0]
    A=np.diag([.02,.08,.4])
    B=(V*np.array([.01,.07,.25]))@V.conj().T
    A0=np.diag([.02,.02,.4])
    B0=(V*np.array([.01,.01,.25]))@V.conj().T
    def trace_product(A,B):
        return np.trace(word(A,B,'AABABB'))*np.trace(word(A,B,'BAA'))/3
    v=[trace_product(x,y).imag for x,y in ((A,B),(A0,B),(A,B0),(A0,B0))]
    assert abs(v[0])>1.e-12 and max(abs(z) for z in v[1:])<1.e-20
    return dict(J_at_1681_point=J,abs_first_degree_six_CP_trace=measured,
                factorized_value=theory,mixed_difference_example=v)


def solve():
    data=read('1397/inputs.json')['high_scale']
    yb_center,yb_sigma=data['down'][2]
    yb=yb_center+yb_sigma*read('1678/independent_results.json')['joint_witness']['quark_pulls'][2]
    real=yb*np.array(read('1681/results.json')['quark_path_amplitudes'])
    powers=np.array(read('1653/results.json')['powers_A_B_C_E_F_R_D'])
    phi=-math.pi/6
    yu=np.array(data['up'])[:,0]*np.exp(1j*np.array([2,3,3])*phi)
    def down(t):
        a,b,c,e,f,z1,z2=real*np.exp(1j*powers*t)
        return np.array([[0,a,0],[b,c,z1+z2],[0,e,f]]),z1,z2
    y,z1,z2=down(phi)
    if (np.prod(yu)*np.linalg.det(y)).real<0:
        yu[0]*=-1
    assert abs(np.angle(np.prod(yu)*np.linalg.det(y)))<1.e-14
    weak=y@y.conj().T
    bound=absolute_bound()
    cpinv=verify_trace_identities(np.diag(abs(yu)**2),weak)
    shift=2*(z1.conjugate()*z2).imag
    shift_prime=2*(powers[6]-powers[5])*(z1.conjugate()*z2).real
    n=powers[6]-powers[5]
    def states(t):
        y,z,w=down(t)
        cols=real*np.exp(1j*powers*t)*1j*powers
        a,b,c,e,f,r,d=cols
        dy=np.array([[0,a,0],[b,c,r+d],[0,e,f]])
        U=np.diag([1.,1.,np.exp(2j*(z.conjugate()*w).imag)])
        new=y@U
        np.testing.assert_allclose(new@new.conj().T,y@y.conj().T,atol=1.e-18)
        l,s,rt=np.linalg.svd(y)
        diag_vertex=np.diag(l.conj().T@dy@rt.conj().T)
        frob_prime=2*(s*diag_vertex.real).sum()
        return new,frob_prime
    new,vertex=states(phi)
    dtheta=np.angle(np.linalg.det(new)/np.linalg.det(y))
    np.testing.assert_allclose(dtheta,shift,atol=4.e-15)
    np.testing.assert_allclose(-vertex/n,shift,atol=1.e-17)
    step=1.e-4
    def angles(t):
        m,_,_=down(t)
        mt,_=states(t)
        return np.angle(np.linalg.det(mt)/np.linalg.det(m))
    fd=(angles(phi+step)-angles(phi-step))/(2*step)
    np.testing.assert_allclose(shift_prime,fd,atol=4.e-12)
    edm=read('1389/results.json')['conditional_edm']['theta_limit_rad']
    return dict(all_order=bound,word_check=cpinv,power_difference=int(n),tree_bar_theta_rad=0.,
                resolved_unit_kappa_shift_rad=shift,angular_slope_delta_bar_theta=shift_prime,
                down_squared_singular_value_sum_angular_slope=vertex,
                kappa_limit_on_small_branch_single_EDM_source=edm/abs(shift),
                relative_left_gram_error=float(np.linalg.norm(new@new.conj().T-weak)/np.linalg.norm(weak)))


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    out=solve()
    path=HERE/'independent_results.json'
    if args.check:
        def check(x,y):
            if isinstance(x,dict):
                for k in x: check(x[k],y[k])
            elif isinstance(x,str): assert x==y
            elif np.isscalar(x) and 0<abs(x)<1.e-23:
                np.testing.assert_allclose(x,y,rtol=1.e-8,atol=1.e-36)
            else: np.testing.assert_allclose(x,y,rtol=1.e-8,atol=1.e-14)
        check(out,json.loads(path.read_bytes()))
        print("PASS: mixed-gap finite bound, first CP trace and physical two-path threshold")
    else:
        with path.open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
