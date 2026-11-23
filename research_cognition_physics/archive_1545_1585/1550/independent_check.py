"""Exact neutral spectrum tests for the second Dirac singular direction."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def mixing(s12,s13,s23,phase):
    c12,c13,c23=np.sqrt(1-np.array([s12,s13,s23])**2)
    p=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/p],[-s12*c23-c12*s23*s13*p,
            c12*c23-s12*s23*s13*p,s23*c13],[s12*s23-c12*c23*s13*p,
            -c12*s23-s12*c23*s13*p,c23*c13]])


def sigma(matrix):
    return np.linalg.svd(matrix,compute_uv=False)


def inverse_root(matrix):
    w,vec=np.linalg.eigh(matrix)
    return (vec/np.sqrt(w))@vec.conj().T


def neutral_example():
    light_ev=np.sqrt(.006**2+np.array([0.,.0000749,.002513]))
    light=light_ev*1e-9
    heavy=np.array([3000.,5000.,8000.])
    U=mixing(*np.sqrt([.308,.02215,.55]),math.radians(200))
    def build(y):
        c,s=math.cosh(y),math.sinh(y)
        O=np.array([[c,1j*s,0],[-1j*s,c,0],[0,0,1]])
        X=1j*np.sqrt(light)[:,None]*O/np.sqrt(heavy)[None,:]
        n=inverse_root(np.eye(3)+X@X.conj().T)
        b=inverse_root(np.eye(3)+X.conj().T@X)
        W=np.block([[U@n,U@n@X],[-b@X.conj().T,b]])
        masses=np.r_[light,heavy]
        M=(W.conj()*masses)@W.conj().T
        return W,M
    lo,hi=0.,19.
    for _ in range(60):
        mid=(lo+hi)/2
        w,mat=build(mid)
        if sigma(mat[:3,3:])[0]<150.:
            lo=mid
        else:
            hi=mid
    w,mat=build((lo+hi)/2)
    np.testing.assert_allclose(w.conj().T@w,np.eye(6),atol=2e-14)
    assert sigma(mat[:3,:3])[0]<8e-14
    N,S=w[:3,:3],w[:3,3:]
    np.testing.assert_allclose(inverse_root(N@N.conj().T)@N,U,atol=2e-13)
    spec=sigma(mat[:3,3:]); x=S@S.conj().T
    nn=N@N.conj().T
    nwidth=float(np.trace(nn@nn).real)
    assert spec[1]**2<=light[-1]*(heavy[-1]+light[-1])
    assert 2.9963-2*.0074<=nwidth<=2.9963+2*.0074
    return {'light_eV':light_ev.tolist(),'heavy_GeV':heavy.tolist(),'complex_boost':(lo+hi)/2,
            'sigma_D_GeV':spec.tolist(),'heavy_active_row_losses':x.diagonal().real.tolist(),
            'Nnu_independent_normalization':nwidth,'polar_pmns_max_nonunitary_deviation':float(np.max(abs(N-U)))}


def rank_one_charged():
    table=json.loads((HERE.parent/'1547'/'inputs.json').read_text(encoding='utf-8'))
    vals=np.array(table['light_scales']['1000'])[:,0]
    down,charged=vals[:3],vals[3:]
    V=mixing(*np.sin([.22704,.00377,.04275]),1.139)
    F=V.conj().T@np.diag([6.15e-6,.00311,.8616])-np.diag(down)
    a=np.sort(sigma(F)); e=charged
    h=[0.]
    for k in (0,1):
        left,right=max(a[k],e[k]), min(a[k+1],e[k+1])
        assert 0<left<right
        h.append(left*right)
    h=np.array(h)
    def template(ascending):
        squares=[-np.prod(q-ascending**2)/np.prod(np.delete(q-h,i)) for i,q in enumerate(h)]
        assert min(squares)>0
        return np.array([[0,np.sqrt(h[1]),0],[0,0,np.sqrt(h[2])],np.sqrt(squares)])
    f0,e0=template(a),template(e)
    l,_,rr=np.linalg.svd(F); l0,_,rr0=np.linalg.svd(f0)
    E=-l@l0.conj().T@e0@rr0.conj().T@rr
    np.testing.assert_allclose(np.sort(sigma(E)),e,rtol=1e-9)
    spec=sigma(F+E)
    assert spec[1]<1e-13 and spec[0]>.8
    return {'charged_measured_yukawas':e.tolist(),'singulars_up_minus_down':a.tolist(),
            'rank_one_dirac_yukawa':spec.tolist(),'intermediate_squared_norms':h.tolist()}


def calculate():
    tab=json.loads((HERE.parent/'1547'/'inputs.json').read_text(encoding='utf-8'))
    rows=np.array(tab['light_scales']['1000']); maximum=rows[:,0]+2*rows[:,1]
    delta=maximum[:3]+maximum[3:]
    yc,yt=.00311-2*.00005,.8616-2*.0043
    cb=math.sin(.04275+2*.00042)
    qc=np.linalg.norm(delta*np.array([1,1,cb])); qt=max(delta)
    alpha=math.hypot(qc/yc,qt/yt)
    y2=yc*(1-alpha)
    # Common, stated lower bound on v(1 TeV) and matched threshold; not the pole VEV.
    v_min=230.
    D2=y2*v_min/math.sqrt(2); m=.5e-9
    return {'absolute_light_mass_cap_eV':m*1e9,'if_heavy_cap_10TeV_D2_upper_MeV':1e3*math.sqrt(m*(1e4+m)),
            'if_D2_at_least_point1_GeV_heaviest_lower_GeV':.1**2/m-m,
            'paired_both_chiralities':{'common_v_lower_GeV':v_min,'two_row_loss_norm':alpha,
            'second_dirac_yukawa_lower':y2,'second_dirac_mass_lower_GeV':D2,'heavy_mass_lower_GeV':D2**2/m-m},
            'charged_rank_one_without_pairing':rank_one_charged(),'six_state_finite_example':neutral_example()}


def compare(now,saved):
    if isinstance(now,dict):
        assert now.keys()==saved.keys()
        for key in now:
            compare(now[key],saved[key])
    else:
        np.testing.assert_allclose(now,saved,rtol=2e-8,atol=2e-12)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check',action='store_true')
    mode.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate(); dest=HERE/'independent_results.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False); stream.write('\n')
    else:
        compare(result,json.loads(dest.read_text(encoding='utf-8')))
        print('1550 independent: second singular direction, charged completion and exact six-state spectrum verified')
