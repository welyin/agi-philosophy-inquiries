"""A shared, non-unit neutral transfer tested on the full charged and neutral mass pair."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ARCH=HERE.parent
BASE=runpy.run_path(str(ARCH/'1679/independent_check.py'))
OLD=BASE['OLD']
PHI=-np.pi/6
ZL=1/math.sqrt(3)
ZN=math.sqrt(7)/5
AA=np.outer(BASE['PIV']['A'],BASE['PIV']['A'])
BB=np.outer(BASE['PIV']['B'],BASE['PIV']['B'])
ETA=np.exp(2j*np.pi/3)
UNIT=np.exp(1j*np.r_[0.,2*PHI,3*PHI+np.angle(3+2*np.exp(4j*PHI))])


def stat(r,x):
    return r*(x+1)/(x+r*r)


def neutral_ratio(s):
    x=ZL*(1-s*ZL)/(s-ZL)
    assert np.all(x>0)
    return x,stat(ZN,x)


def roots(h,R):
    aa,bb,ab=1+h*h,10+h*h,3+h*h
    co=(aa*bb-ab*ab)*(R+1/R)+ab*ab
    disc=co**2-4*(aa*bb)**2
    assert np.all(disc>0)
    return 2*aa**2/(co+np.sqrt(disc)),(co+np.sqrt(disc))/(2*bb**2)


def modes(h,r):
    phase=np.broadcast_to(UNIT,(len(h),3)).copy()
    phase[:,2]*=h
    n=phase[:,:,None]*(AA+ETA*r[:,None,None]*BB)*phase[:,None,:]
    return np.linalg.eigh(n.conj().transpose(0,2,1)@n)[1]


def point(arr):
    old=BASE['angles_frame'](arr,0,1,63.2)
    R=old['R']
    s=abs(old['t'])
    x,h=neutral_ratio(s)
    ref=modes(np.ones_like(R),BASE['PIV']['roots'](R)[0])
    ud=old['pmns']@ref.conj().transpose(0,2,1)
    r=roots(h,R)[0]
    u=ud@modes(h,r)
    p=abs(u)**2
    pe=p[:,0,2]
    obs=np.array([4*pe*(1-pe),p[:,0,1]/(1-pe),p[:,1,2]/(1-pe)]).T
    return old,obs,h,x


def cost(x):
    old,obs,_,_=point(x)
    a=BASE['ANGLES']
    return np.maximum(np.abs(old['pulls']).max(axis=1)/2,
                      (abs(obs-a.mean(axis=1))/(a[:,1]-a.mean(axis=1))).max(axis=1))


def search():
    rg=np.random.default_rng(1682)
    hi=np.r_[[2.]*6,1.8,1.]
    lo=-hi
    n=300
    a=rg.uniform(lo,hi,(n,8))
    v=cost(a)
    for _ in range(350):
        b,c,d=a[rg.integers(0,n,size=(3,n))]
        y=np.clip(b+rg.uniform(.45,.9)*(c-d),lo,hi)
        y=np.where(rg.random((n,8))<.9,y,a)
        vy=cost(y)
        use=vy<v
        a[use],v[use]=y[use],vy[use]
    j=np.argmin(v)
    print(v[j],a[j].tolist())
    return a[j]


def evaluated(label,data,width,frequency_width):
    xarr=np.array(data)[None,:]
    old,obs,hs,xs=point(xarr)
    s,h,x,R=[float(v[0]) for v in (abs(old['t']),hs,xs,old['R'])]
    assert np.max(abs(xarr[0,:6]))<width and np.max(abs(old['pulls']))<width
    limits=BASE['ANGLES']
    assert np.all((limits[:,0]<obs[0])&(obs[0]<limits[:,1]))
    raw=old['raw'][0]
    sh,long,v,delta=OLD['get_matrices'](raw,63.2)
    E,y,_=np.linalg.svd((sh+s*np.exp(1j*PHI)*long).T*old['scale'][0])
    E,y=E[:,::-1],y[::-1]
    np.testing.assert_allclose(y,BASE['L'][:,0]+old['pulls'][0]*BASE['L'][:,1],rtol=1.e-10)
    cfg=OLD['NEUTRINO']['INPUT']
    ds_lo,ds_hi=OLD['NEUTRINO']['window'](cfg['juno_NO_dm21'],frequency_width)
    da_lo,da_hi=OLD['NEUTRINO']['window'](cfg['juno_NO_dm31'],frequency_width)
    lo,hi=max(ds_lo,da_lo/R**2),min(ds_hi,da_hi/R**2)
    assert lo<hi
    ds=(lo+hi)/2
    r=roots(h,R)[0]
    legs=UNIT.copy()
    legs[2]*=h
    m=legs[:,None]*(AA+ETA*r*BB)*legs
    vals,w=OLD['NEUTRINO']['OLD']['takagi'](m)
    vals=vals*math.sqrt(ds)/vals[1]
    u=E.conj().T@w
    check=OLD['NEUTRINO']['OLD']['angles'](u)
    np.testing.assert_allclose(check,obs[0],atol=2.e-11)
    np.testing.assert_allclose(vals[-1]**2,R**2*ds,rtol=1.e-12)
    ice=OLD['ICE']['interpolate'](*OLD['ICE']['grid'](),check[2],R**2*ds-ds)
    for raw_amp in (1,ZL,ZN):
        block=np.array([[np.sqrt(x),raw_amp],[-raw_amp,np.sqrt(x)]])
        trans=-np.linalg.inv(block)[0,1]
        np.testing.assert_allclose(trans*(1+x),stat(raw_amp,x))
    slopes=np.array([4/math.sqrt(3),12*math.sqrt(3)/7])
    rr=np.array([ZL,ZN])
    log_phi=(x-rr**2)/(x+rr**2)*slopes
    log_x=1/(1+x)-1/(x+rr**2)
    # Change even the common mass parameter with the flavon. This combination
    # then still isolates the return-amplitude response for the specified paths.
    def combined(tau,dx):
        omega=np.exp(4j*(PHI+tau))
        mod=np.array([abs(2+omega)/3,abs(3+2*omega)/5])
        return np.log(stat(mod,x+dx*tau)) @ np.array([log_x[1],-log_x[0]])
    free=float(log_x[1]*log_phi[0]-log_x[0]*log_phi[1])
    for dx in (0.,-.12,.28):
        numerical=(combined(1.e-5,dx)-combined(-1.e-5,dx))/(2.e-5)
        np.testing.assert_allclose(numerical,free,atol=3.e-9)
    return dict(title=label,quark_marginal_pulls=data[:6],lepton_marginal_pulls=old['pulls'][0].tolist(),
                gamma_matching=OLD['CKM']['gamma'](v),mag_t=s,x=x,neutral_modulus=h,
                charged_normalization=float(old['scale'][0]),R=R,source_r=r,
                d21_eV2=ds,d31_eV2=ds*R**2,three_angles=check.tolist(),
                neutrino_delta_over_pi=OLD['NEUTRINO']['phase_delta'](u),
                IceCube_frozen_other_coordinates_q=float(ice),
                common_x_response=dict(fixed_x_log_slope_charged_neutral=log_phi.tolist(),
                                       coefficients_for_dx=log_x.tolist(),x_variation_cancelled_slope=free))


def compute():
    unit_gain=(1+ZN)/(2*np.sqrt(ZN))
    rcrit=math.sqrt(ZN)
    np.testing.assert_allclose(stat(rcrit,ZN),unit_gain)
    prev=OLD['read']('../1678/independent_results.json')['joint_witness']
    xold,oldh=neutral_ratio(prev['common_abs_t'])
    # These points were found after the mass windows were known. Numbers here
    # are rounded before direct verification through both full mass matrices.
    first=evaluated('one width charged; two widths JUNO',
        [-.997,-.997,.99668,-.997,-.75902,.23730,-.3345202,-.98333333],1,2)
    second=evaluated('one and a half widths charged; one width JUNO',
        [-1.499,-1.499,1.499,-1.499,1.499,.33046,-.503534,-.48333333],1.5,1)
    assert first['mag_t']>unit_gain and second['mag_t']>unit_gain
    return dict(static_max_for_exact_neutral_unit=unit_gain,charged_return_at_max=rcrit,
                original_1678_x=xold,original_1678_new_neutral_modulus=oldh,witnesses=[first,second])


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    cli.add_argument('--search',action='store_true')
    args=cli.parse_args()
    if args.search:
        search()
    else:
        out=compute()
        target=HERE/'independent_results.json'
        if args.check:
            def close(x,y):
                if isinstance(x,dict):
                    for k in x: close(x[k],y[k])
                elif isinstance(x,list):
                    for a,b in zip(x,y): close(a,b)
                elif isinstance(x,str): assert x==y
                else: np.testing.assert_allclose(x,y,rtol=1.e-10,atol=1.e-12)
            close(out,json.loads(target.read_bytes()))
            print('PASS: common Dirac kernel, charged singular values, Takagi, both reactor gaps and angles')
        else:
            with target.open('x',encoding='utf8') as f:
                f.write(json.dumps(out,indent=2)+'\n')
