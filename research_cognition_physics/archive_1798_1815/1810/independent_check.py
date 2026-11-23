"""Inverse-spectrum one-singlet test and a conditioned CP/neutral-current bound."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
EDGES=[(0,1),(1,2),(2,0)]


def edge(a):
    return np.array([a[i,j] for i,j in EDGES])


def cp(a):
    return float(np.prod(edge(a)).imag)


def adj(a):
    return np.linalg.det(a)*np.linalg.inv(a)


def inverse_poly(light,m2):
    z=np.eye(3)-light@light.conj().T
    k=(light/m2)@light.conj().T
    p=np.array([1.+0j])
    for a,b in zip(edge(k),edge(z)):
        p=np.polynomial.polynomial.polymul(p,[a,b])
    assert abs(p[-1].imag)<1e-15
    return p.imag[:3]


def complete(m2,V):
    K=(V/m2)@V.conj().T
    assert abs(cp(K)) < 1e-12
    # Build a real gauge along a spanning forest, including zero-edge cases.
    phase=np.zeros(3,complex)
    for start in range(3):
        if phase[start]!=0:
            continue
        phase[start]=1
        pending=[start]
        while pending:
            i=pending.pop()
            for j in range(3):
                if phase[j]==0 and abs(K[i,j])>1e-12:
                    phase[j]=phase[i]*np.exp(-1j*np.angle(K[i,j]))
                    pending.append(j)
    rows=phase.conj()[:,None]*V
    nu=np.linalg.svd(rows,full_matrices=True)[2][-1]
    U=np.vstack([rows,nu])
    np.testing.assert_allclose(U@U.conj().T,np.eye(4),atol=2e-15)
    H=(U*m2)@U.conj().T
    m0=math.sqrt(float(H[3,3].real))
    b=H[:3,3:4]/m0
    s=H[:3,:3]-b@b.conj().T
    assert np.max(abs(s.imag))<1e-10
    lam,o=np.linalg.eigh(s.real)
    assert min(lam)>0
    m=(o*np.sqrt(lam))@o.T
    M=np.block([[m,b],[np.zeros((1,3)),np.array([[m0]])]])
    return np.linalg.norm(M@M.conj().T-H),float(np.linalg.det(M).real)


def example():
    m=np.array([[1,0,0],[1,2,0],[1,1,3.]])
    b=np.array([[1],[1j],[1+1j]])
    M=np.block([[m,b],[np.zeros((1,3)),np.array([[4.]])]])
    h=M@M.conj().T
    m2,U=np.linalg.eigh(h)
    V=U[:3,:]
    coefficients=inverse_poly(V[:,:3],m2[:3])
    roots=np.sort(np.roots(coefficients[::-1]))
    assert min(abs(roots-1/m2[-1]))<2e-13
    positive=np.array(roots[(roots>0)&(roots<1/m2[2])])
    # V keeps all coherent weak data while the heavy mass changes.
    bad=m2.copy();bad[-1]*=2
    ibad=cp((V/bad)@V.conj().T)
    assert abs(ibad)>1e-8
    err,det=complete(m2,np.exp(1j*np.array([.2,1.,-.7]))[:,None]*V)
    assert err<1e-12
    assert abs(cp(h[:3,:3]))>1e-3
    return dict(heavy_mass_squared=m2[-1],light_squared=m2[:3].tolist(),
                inverse_quadratic_ascending=coefficients.tolist(),roots_in_inverse_squared_mass=roots.tolist(),
                allowed_heavy_squared_mass=(1/positive).tolist(),
                squared_mass_times_two_CP=ibad,lower_triangle_CP_on_same_physical_matrix=cp(h[:3,:3]),
                weak_quartet=(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag,
                det_positive_real=det,reconstruction_error_below_1e12_inverse=bool(err<1e-12))


def threshold(function,goal):
    lo,hi=0.,.5
    assert function(hi)>goal
    for _ in range(75):
        mid=(lo+hi)/2
        if function(mid)<goal:
            lo=mid
        else:
            hi=mid
    return lo


def physical():
    api=runpy.run_path(str(HERE.parent/'1799/check.py'))
    _,down,_,v,_=api['backend']()
    eig=(down[0]/down)**2
    K=(v*eig)@v.conj().T
    A=(v*np.array([eig[1]*eig[2],eig[0]*eig[2],eig[0]*eig[1]]))@v.conj().T
    np.testing.assert_allclose(A,adj(K),atol=1e-16)
    inv=abs(cp(K))
    z=edge(K);za=edge(A)
    for i,j in EDGES:
        assert abs((K[i,j]*A[j,i]).imag-cp(K))<1e-15
    J=abs((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    assert abs(inv-J*math.prod(abs(eig[i]-eig[j]) for i,j in EDGES))<1e-16
    results={}
    def edge_norms(mat):
        lengths=np.linalg.norm(mat,axis=0)
        return np.array([(lengths[i]+lengths[j]+2*abs(mat[i,j]))/2 for i,j in EDGES])
    linear_k=edge_norms(K)
    linear_a=edge_norms(A)
    for label,h in [('only_heavier_than_b',eig[2]),('infinite_mass_limit',0.)]:
        def errors(d):
            a=d/(1+math.sqrt(1-d))
            dk=2*a+h*d
            da=d/(1-d)*eig[1]+2*h*d/(1-d)**2
            return dk,da
        def cubic(d):
            dk,_=errors(d)
            return np.prod(abs(z)+dk)-np.prod(abs(z))
        def joint(d):
            dk,da=errors(d)
            return min(abs(z)*da+abs(za)*dk+dk*da)
        def refined(d):
            a=d/(1+math.sqrt(1-d))
            b=a/math.sqrt(1-d)
            dk=a*linear_k+(a*a+h*d)/2
            da=b*linear_a+b*b*eig[1]/2+2*h*d/(1-d)**2
            return min(abs(z)*da+abs(za)*dk+dk*da)
        db=threshold(cubic,inv)
        dc=threshold(joint,inv)
        de=threshold(refined,inv)
        amp=math.sqrt(de*(1-de))
        bound=math.sqrt(4*math.pi)*246.22/math.sqrt(2)/amp
        # Conservative common-scale upper envelope m_b <= 5 GeV.
        mhmax=(bound+math.hypot(bound,10.))/2
        results[label]=dict(triple_edge_lower_d=db,cofactor_operator_lower_d=dc,cofactor_edge_lower_d=de,
                            light_heavy_Z_total_squared_lower_if_d_le_half=amp**2,
                            scaled_Z_amplitude_sqrt_lower_if_d_le_half=amp,
                            max_heavy_mass_GeV_if_EW_y_norm_cap_sqrt4pi_and_mb_le_5=mhmax)
    return dict(inverse_eigenvalues_normalized=eig.tolist(),absolute_J=J,absolute_CP=inv,
                edge_absolute=abs(z).tolist(),cofactor_edge_absolute=abs(za).tolist(),
                finite_tests=results)


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    if isinstance(x,float):
        return 0. if abs(x)<1e-13 else float(format(x,'.10g'))
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=HERE/'independent_results.json'
    raw=json.dumps(stable(dict(constructed=example(),quark_inputs=physical())),indent=2,allow_nan=False)+'\n'
    if args.check:
        assert json.loads(raw)==json.loads(out.read_bytes())
    else:
        out.write_bytes(raw.encode())
    print('1810 independent',hashlib.sha256(raw.encode()).hexdigest()[:12])
