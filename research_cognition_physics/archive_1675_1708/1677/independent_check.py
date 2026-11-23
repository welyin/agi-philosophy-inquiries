"""Signed one-loop quartet flow and a renormalizable finite matching counterexample."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
OLD = runpy.run_path(str(ARCH/"1675/independent_check.py"))


def read(path):
    return json.loads((ARCH/path).read_bytes())


def generator(spectrum, source):
    a=spectrum
    t=np.zeros((3,3),complex)
    for i in range(3):
        for k in range(3):
            if i!=k:
                t[i,k]=-1.5*(a[i]+a[k])/(a[k]-a[i])*source[i,k]
    np.testing.assert_allclose(t,-t.conj().T,atol=1.e-15)
    return t


def compare(a,b,v):
    tu=generator(a,(v*b)@v.conj().T)
    td=generator(b,(v.conj().T*a)@v)
    direct=(-tu@v+v@td)/v
    J=(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    out=np.zeros((3,3))
    for i,j in np.ndindex(3,3):
        k,l=(i+1)%3,(i+2)%3
        m,n=(j+2)%3,(j+1)%3
        A=(a[k]-a[i])*(a[l]-a[i])
        B=(b[m]-b[j])*(b[n]-b[j])
        out[i,j]=-3*J*(a[k]-a[l])*(b[m]-b[n])*(a[i]/A+b[j]/B)/abs(v[i,j])**2
    np.testing.assert_allclose(out,direct.imag,atol=4.e-14,rtol=2.e-10)
    signed=out[0,0]-out[0,2]-out[1,0]+out[1,2]
    assert J*signed<=0
    return float(signed),float(np.max(abs(out-direct.imag)))


def sign_bounds(tube):
    eps_u,eps_c=tube["adjacent_up_squared_ratio_bounds"]
    eps_d,eps_s=tube["adjacent_down_squared_ratio_bounds"]
    v=tube['running_CKM_magnitude_bounds']
    # |d log(bs/bb)/d tau| <= 3(at+bb), tau=log(mu)/(16 pi^2).
    data=np.array(read('1397/inputs.json')['high_scale']['down'])
    bs_bb_0=((data[1,0]-2*data[1,1])/(data[2,0]+2*data[2,1]))**2
    integrals=tube['int_largest_up_squared_upper']+tube['int_largest_down_squared_upper']
    bs_bb_floor=bs_bb_0*math.exp(-3*integrals)
    # Round towards a weaker decimal input before evaluating the two tests.
    bs_bb_floor=math.floor(bs_bb_floor*1.e9)/1.e9
    up_adverse=(eps_u/v['ud_min']**2
                + eps_s/((1-eps_s)*v['cb'][0]**2))
    ratio=eps_d/(bs_bb_floor*(1-eps_d)**2)
    down_adverse=ratio*v['ub'][1]**2*(1/v['ud_min']**2
                                      +1/(v['cd_min']**2*(1-eps_c)))
    assert up_adverse<.4 and down_adverse<.05
    return dict(bs_bb_floor=bs_bb_floor,up_adverse_fraction_upper=math.ceil(up_adverse*1.e7)/1.e7,
                down_adverse_fraction_upper=math.ceil(down_adverse*1.e7)/1.e7)


def finite_match():
    point=read('1675/results.json')['two_width_overlapping_direct_phase']
    rows=point['quark_yukawas_and_angles']
    masses=np.array(rows[:3])
    theta12,theta13,theta23=rows[3:]
    center=math.radians(63.2)
    k=math.cos(theta12)*math.tan(theta23)*math.sin(theta13)/math.sin(theta12)
    new_delta=center+math.asin(k*math.sin(center))
    high=OLD['texture'](rows,math.acos(point['cos_delta']))[1]
    low=OLD['texture'](rows,new_delta)[1]
    initial,target=high*masses,low*masses
    change=target-initial
    _,sv,vh=np.linalg.svd(change)
    assert sv[1]>1.e-9 and sv[2]<1.e-16
    # Complete the two-row right-hand mixing R at the threshold. Its rows
    # span the entire correction, including canonical right-field normalization.
    kappa=.01
    R=kappa*vh[:2]
    P=R.conj().T@R/kappa**2
    C=np.eye(3)+(math.sqrt(1+kappa*kappa)-1)*P
    L=(initial-target@C)@R.conj().T/kappa**2
    after=np.linalg.solve(C.T,(initial-L@R).T).T
    np.testing.assert_allclose(after,target,atol=1.e-16)
    gamma0=OLD['gamma'](high)
    gamma1=OLD['gamma'](low)
    p0=abs(high)**2
    p1=abs(low)**2
    np.testing.assert_allclose(p1[[0,0,1],[1,2,2]],p0[[0,0,1],[1,2,2]],atol=1.e-15)
    check=[]
    for heavy_over_vev in (100.,1000.):
        full=np.block([[initial,L],[heavy_over_vev*R,heavy_over_vev*np.eye(2)]])
        u,s,_=np.linalg.svd(full)
        light=u[:3,::-1][:,:3]
        light_y=s[::-1][:3]
        np.testing.assert_allclose(light_y,masses,rtol=2.e-7)
        assert abs(OLD['gamma'](light)-gamma1)<1.e-5
        check.append([heavy_over_vev,OLD['gamma'](light),max(abs(light_y/masses-1))])
    return dict(high_gamma=gamma0,matched_gamma=gamma1,canonical_shift_deg=gamma1-gamma0,
                relative_correction_norm_over_yb=sv[0]/masses[2],two_nonzero_singular_values=sv[:2].tolist(),
                right_mixing_norm=kappa,left_yukawa_norm=np.linalg.norm(L,2),
                full_mass_svd_M_over_vev_gamma_relative_mass_error=check)


def main():
    data=read('1397/inputs.json')['high_scale']
    u=np.array(data['up'])[:,0]**2
    rows=np.array(data['down']+data['angles_theta12_theta13_theta23_delta'][:3])[:,0]
    b=rows[:3]**2
    ds=[-2.7,-1.24,-.1,0.,.1,1.24,2.7]
    diffs=[]
    rates=[]
    for delta in ds:
        v=OLD['texture'](rows,delta)[1]
        rate,error=compare(u,b,v)
        rates.append(rate)
        diffs.append(error)
    return dict(normalized_high_log_slope_vs_delta=rates,tests_delta_radians=ds,
                max_absolute_quartet_generator_error=max(diffs),
                full_tube_signed_proof=sign_bounds(read('1677/results.json')),finite_matching=finite_match())


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=main()
    if args.check:
        def check(x,y):
            if isinstance(x,dict):
                for key in x: check(x[key],y[key])
            else:
                np.testing.assert_allclose(x,y,rtol=1.e-8,atol=5.e-12)
        check(result,read('1677/independent_results.json'))
        print("PASS: signed quartet evolution, all-history sign and two-heavy finite matching")
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as f:
            f.write(json.dumps(result,indent=2)+'\n')
