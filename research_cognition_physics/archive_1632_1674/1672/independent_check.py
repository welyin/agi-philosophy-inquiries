"""Verify adjoint compression, leakage, and a mass-ratio-independent null-mode obstruction."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def common_exchange(p):
    q=abs(p[2])
    theta=np.pi/(1+q)
    # Independent full Hamiltonian exponentiation, ordered by input/output port.
    h=np.zeros((6,6),dtype=complex)
    h[3:,:3]=np.diag(p)
    h[:3,3:]=h[3:,:3].conj().T
    vals,vecs=np.linalg.eigh(h)
    unitary=(vecs*np.exp(-1j*theta*vals))@vecs.conj().T
    success=unitary[3:,:3]
    failure=unitary[:3,:3]
    desired=np.diag(p/abs(p))
    s=np.sin(theta)
    f=np.cos(q*theta)
    np.testing.assert_allclose(success,-1j*s*desired,atol=1e-14)
    np.testing.assert_allclose(failure,f*np.diag([-1.,-1.,1.]),atol=1e-14)
    np.testing.assert_allclose(success.conj().T@success,s*s*np.eye(3),atol=1e-14)
    np.testing.assert_allclose(failure.conj().T@failure,f*f*np.eye(3),atol=1e-14)
    a,b=np.array([0,1,1]),np.array([1,3,1])
    mass=np.outer(a,a)+.2j*np.outer(b,b)
    np.testing.assert_allclose(success.T@mass@success,-s**2*(desired.T@mass@desired),atol=1e-14)
    # Each constant pulse changes only the common physical phase of the same six-port coupling.
    beta=float(np.arccos((1/np.tan(theta)**2-1)/2))
    middle=h.copy()
    middle[3:,:3]*=np.exp(1j*beta)
    middle[:3,3:]=middle[3:,:3].conj().T
    vals,vecs=np.linalg.eigh(middle)
    rotated=(vecs*np.exp(-1j*theta*vals))@vecs.conj().T
    full=unitary@rotated@unitary
    amp=-1j*s*(f*f*(2+np.exp(1j*beta))-s*s*np.exp(-1j*beta))
    np.testing.assert_allclose(abs(amp),1.,atol=1e-14)
    np.testing.assert_allclose(full[:3,:3],np.zeros((3,3)),atol=1e-14)
    np.testing.assert_allclose(full[3:,:3],amp*desired,atol=1e-14)
    np.testing.assert_allclose(full[3:,:3].T@mass@full[3:,:3],amp**2*(desired.T@mass@desired),atol=1e-14)
    return {'common_pulse_gt':round(float(theta),11),'success_probability_per_attempt':round(float(s*s),11),
            'mean_attempts_if_reset_correction_available':round(float(1/(s*s)),11),
            'failure_probability':round(float(f*f),11),
            'three_pulses_common_phases_rad':[0,round(beta,11),0],
            'three_pulse_common_output_phase_rad':round(float(np.angle(amp)),11),
            'coherent_transfer_probability':1,'minimal_number_of_equatorial_pulses':3}


def solve():
    phi=-np.pi/6
    omega=np.exp(4j*phi)
    adj=np.diag([2.,2.,2.,-3.,-3.])/np.sqrt(30)
    singlet=np.eye(5)/np.sqrt(5)
    link=np.diag([omega]*3+[1.]*2)
    image=link@adj
    z=np.vdot(adj,image)
    loss=np.vdot(singlet,image)
    np.testing.assert_allclose(np.linalg.det(link),1,atol=1e-14)
    np.testing.assert_allclose(link@adj@link.conj().T,adj,atol=1e-14)
    np.testing.assert_allclose(image,z*adj+loss*singlet,atol=1e-14)
    np.testing.assert_allclose(z,(3+2*omega)/5,atol=1e-14)
    np.testing.assert_allclose([abs(z)**2,abs(loss)**2],[7/25,18/25],atol=1e-14)

    charge=runpy.run_path(str(HERE.parent/'1659/independent_check.py'))['charged']
    ue,_=charge(read('../1657/independent_results.json')['extended_box_GJ_test'])
    e=ue[:,0]
    a,b=np.array([0,1,1]),np.array([1,3,1])
    n=np.cross(a,b)
    p=np.array([1.,np.exp(2j*phi),np.exp(3j*phi)*z])
    v=n/p
    v/=np.linalg.norm(v)
    actual=abs(np.vdot(e,v))**2
    tested=[]
    for r in (1.e-4,.07,.5,1.,12.,200.):
        for phase in (0.,np.pi/7,2*np.pi/3,np.pi):
            m=p[:,None]*(np.outer(a,a)+r*np.exp(1j*phase)*np.outer(b,b))*p
            _,sv,vh=np.linalg.svd(m)
            assert sv[1]>1e-8 and sv[2]<1e-11
            assert np.linalg.norm(m@v)<2.e-12
            tested.append(abs(abs(np.vdot(e,vh[-1].conj()))**2-actual))
    assert max(tested)<1.e-9

    bounds=read('../1659/results.json')['joint_windows']
    (ilo,ihi),(slo,shi)=bounds[:2]
    low_s13=np.array([(1-np.sqrt(1-t))/2 for t in (ilo,ihi)])
    high_s13=1-low_s13[::-1]
    def electron_cut(s13):
        return np.array([(1-s13[1])*(1-shi),(1-s13[0])*(1-slo)])
    low,high=electron_cut(low_s13),electron_cut(high_s13)
    c=abs(e[0])
    s=np.linalg.norm(e[1:])
    x=abs(v[0])
    y=np.linalg.norm(v[1:])
    blind=np.square([max(0.,c*x-s*y),min(1.,c*x+s*y)])
    assert high[1]<blind[0]<actual<blind[1]<low[0]
    # Invert the angular-distance triangle bound on the physical small-theta13 branch.
    lower_axis=c*np.sqrt(low[0])-s*np.sqrt(1-low[0])
    min_kappa=1/np.sqrt(4/lower_axis**2-5)

    # For the fixed published phases, electron-null weight as a function
    # of arbitrary third-leg magnitude k is |k*A+B|^2/(5*k^2+1).
    fixed=np.array([1.,np.exp(2j*phi),np.exp(3j*phi)*z/abs(z)])
    A=e[:2].conj()@(n[:2]/fixed[:2])
    B=e[2].conj()*n[2]/fixed[2]
    def weight(k):
        return abs(k*A+B)**2/(5*k*k+1)
    def intersections(target):
        roots=np.roots([abs(A)**2-5*target,2*(A*B.conjugate()).real,abs(B)**2-target])
        return [v.real for v in roots if abs(v.imag)<1e-11 and 0<v.real<1]
    cuts=sorted([0.,1.,*intersections(low[0]),*intersections(low[1])])
    allowed=[[left,right] for left,right in zip(cuts,cuts[1:]) if low[0]<weight((left+right)/2)<low[1]]
    assert len(allowed)==1 and allowed[0][0]>min_kappa and allowed[0][1]==1
    assert abs(weight(abs(z))-actual)<1e-12 and low[0]<weight(1)<low[1]
    def rounded(xs):
        return np.asarray(xs).round(11).tolist()
    return {'adjoint_projection_re_im':rounded([z.real,z.imag]),'adjoint_retention':7/25,
            'singlet_leakage':18/25,'unrotated_electron_null_weight':7/15,
            'direct_physical_electron_null_weight':rounded(actual),
            'phase_independent_physical_electron_null_weight_bounds':rounded(blind),
            'required_small_theta13_electron_null_weight':rounded(low),
            'alternative_large_theta13_electron_null_weight':rounded(high),
            'phase_independent_min_third_leg_magnitude_small_theta13':rounded(min_kappa),
            'fixed_phase_magnitude_intervals_small_theta13_contractions':rounded(allowed),
            'phase_only_physical_electron_null_weight':rounded(weight(1.)),
            'compressed_hamiltonian_common_exchange':common_exchange(p),
            'independent_svd_configurations':len(tested)}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1672 independent PASS: adjoint-singlet decomposition, both reactor branches, all mass ratios')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as output:
            output.write(json.dumps(result,indent=2)+'\n')
