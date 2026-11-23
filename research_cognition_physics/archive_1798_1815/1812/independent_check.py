"""Two-projection lower bounds on hidden heavy-current overlap, Yukawa strength and hard CP."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def ckm(sines):
    s12,s13,s23,delta=sines
    c12,c13,c23=[math.sqrt(1-x*x) for x in sines[:3]]
    z=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/z],
                     [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                     [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])


def quantify(point):
    v=ckm(point['high_ckm_sines_delta'])
    spectra=np.asarray(point['high_yukawas'])
    y=spectra[3:6];yb=y[-1];ytau=spectra[-1]
    h=(v*y**2)@v.conj().T
    z=np.array([h[0,1],h[1,2],h[2,0]])
    inv=abs(np.prod(z).imag)
    J=abs((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    gaps=np.array([y[1]**2-y[0]**2,y[2]**2-y[1]**2,y[2]**2-y[0]**2])
    np.testing.assert_allclose(inv,J*np.prod(gaps),rtol=1e-11,atol=0)
    pairs=np.abs(z*np.roll(z,1))
    low,high=inv/np.linalg.norm(pairs),inv/max(pairs)
    D,T=point['D_T_poles_GeV']
    r=D/T
    frame=yb**2+ytau**2
    cap=math.sqrt(4*math.pi)
    def estimates(r):
        ymax=math.sqrt(frame)
        required=low/(2*yb*r)
        remaining=(ymax/required)**2
        # Conditional when remaining>=1 the possible overlap may be zero.
        xi=max(0,1-remaining)
        min_ai_over_T=max(0,low/(2*yb*cap)-r)/5
        return dict(r=r,Y_required_even_without_screening=required,
                    one_minus_xi_upper=remaining,xi_lower=xi,
                    heavy_to_heavy_inclusive_norm_squared_lower=xi**2,
                    light_to_light_trace_and_lower_bound=1+xi**2,
                    min_ImA_over_lightest_T_if_Y_cap=min_ai_over_T)
    # When both actual color poles are equal b is precisely sigma_min(M_T).
    return dict(yb=yb,ytau=ytau,abs_J=J,abs_CP_cycle=inv,
                rephase_optimized_ImHd_norm_bracket=[low,high],upper_over_lower=high/low,
                central_poles_GeV=[D,T],at_recorded_poles=estimates(r),
                critical_r_if_Y_norm_cap=low/(2*yb*cap),
                ImA_GeV_lower_given_cap_and_central_T=estimates(r)['min_ImA_over_lightest_T_if_Y_cap']*T)


def projection_test():
    # Real A, soft CP B=i*epsilon; deliberately nonorthogonal T,D row spaces.
    rng=np.random.default_rng(1812)
    worst=0
    for index in range(36):
        imag=.2*rng.normal(size=(2,5))
        t=rng.normal(size=(2,5))+1j*imag
        d=rng.normal(size=(2,5))*.3+1j*imag
        y=rng.normal(size=(3,5))
        project=[]
        for x in [t,d]:
            q,_=np.linalg.qr(x.conj().T)
            project.append(q@q.conj().T)
        pt,pd=project
        qt,qd=np.eye(5)-pt,np.eye(5)-pd
        xi=np.linalg.norm(pt@pd,2)
        hd=y@qt@y.T
        he=y@qd@y.T
        lhs=np.linalg.norm(hd.imag,2)
        sw=np.linalg.norm(pt-pt.conj(),2)
        sigma=np.linalg.svd(t,compute_uv=False)[-1]
        check=np.sqrt(np.linalg.norm(hd,2))*np.linalg.norm(y,2)*sw
        assert lhs<=check+2e-14
        assert sw<=2*np.linalg.norm(d,2)/sigma+1e-12
        assert np.linalg.eigvalsh(hd+he-(1-xi)*(y@y.T))[0]>-5e-13
        # Nontrivial principal overlaps of complementary subspaces coincide.
        up=np.linalg.eigh(qt)[1][:,-3:];dp=np.linalg.eigh(qd)[1][:,-3:]
        np.testing.assert_allclose(np.linalg.svd(up.conj().T@dp,compute_uv=False)[:2],[1,xi],atol=1e-12)
        worst=max(worst,lhs/check)
    return dict(samples=36,largest_ratio_to_analytic_projection_bound=worst)


def compute():
    points=read('../1811/results.json')
    windows=read('../1809/results.json')
    rows={}
    for key in ('common_MZ','independent_lattice'):
        q=quantify(points[key])
        gap=windows[key]['original_window']['B_log_min_T_over_D'][0]
        rw=math.exp(-gap)
        q['mass_window_only_r_upper']=rw
        lo=q['rephase_optimized_ImHd_norm_bracket'][0]
        q['min_Y_with_that_r_and_same_high_flavor']=lo/(2*q['yb']*rw)
        q['same_high_flavor_1_minus_xi_upper_for_mass_window']=(
            (q['yb']**2+q['ytau']**2)/q['min_Y_with_that_r_and_same_high_flavor']**2)
        rows[key]=q
    return dict(actual_two_points=rows,projection_examples=projection_test())


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    if isinstance(x,float):
        return float(format(x,'.12g'))
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    f=HERE/'independent_results.json'
    raw=json.dumps(stable(compute()),indent=2,allow_nan=False)+'\n'
    if args.check:
        assert json.loads(raw)==json.loads(f.read_bytes())
    else:
        f.write_bytes(raw.encode())
    print('1812 independent',hashlib.sha256(raw.encode()).hexdigest()[:12])
