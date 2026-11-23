"""Audit two frozen branches against T2K 2023 Tables 9 and 14."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def projectors(branch, rvalues):
    charged = runpy.run_path(str(HERE.parent/'1659/independent_check.py'))['charged']
    charged, _ = charged(read('../1657/independent_results.json')['extended_box_GJ_test'])
    a, b = np.array([0,1,1]), np.array([1,3,1])
    null = np.cross(a,b) / np.linalg.norm(np.cross(a,b))
    rvalues = np.atleast_1d(rvalues)
    d = 6*(rvalues+1/rvalues)+16
    r = 8/(d+np.sqrt(d*d-44**2))
    angle = -np.pi/6
    if branch == 'third_leg_3_to_2':
        phases = [0., 2*angle, 3*angle+np.angle(3+2*np.exp(4j*angle))]
        sign = 1
    else:
        phases = [0., 10*angle+np.deg2rad(11.5), 7*angle]
        sign = -1
    basis = charged.conj().T * np.exp(-1j*np.array(phases))
    p1 = np.outer(basis@null, (basis@null).conj())
    mass = np.outer(a,a) + (r*np.exp(sign*2j*np.pi/3))[:,None,None]*np.outer(b,b)
    gram = mass.conj().transpose(0,2,1) @ mass
    lo = 6*r/rvalues
    hi = 6*r*rvalues
    p3 = (basis@gram@basis.conj().T-lo[:,None,None]*(np.eye(3)-p1))/(hi-lo)[:,None,None]
    np.testing.assert_allclose(p3@p3,p3,atol=1e-13)
    return p1, p3


def coordinates(branch, rs, h=0):
    p1,p3=projectors(branch,rs)
    z=p3[:,0,0].real
    mu=p3[:,1,1].real
    j=-(p1[0,1]*p3[:,1,0]).imag
    # ||P3'|| < 2 and P1 constant on this interval (1668, proof).
    eps=2*h+2.e-11
    zlo,zhi=z-eps,z+eps
    tlo,thi=(mu-eps)/(1-zlo),(mu+eps)/(1-zhi)
    slo,shi=1-p1[0,0].real/(1-zhi),1-p1[0,0].real/(1-zlo)
    # atan2(y,x) reconstructs the four-quadrant phase, x and y
    # proportional to cos(delta) and sin(delta) with the same positive factor.
    xlo=p1[1,1].real-shi*(1-tlo)-(1-slo)*thi*zhi
    xhi=p1[1,1].real-slo*(1-thi)-(1-shi)*tlo*zlo
    assert np.all(xlo>0)
    ys=[2*(j+q*eps*abs(p1[0,1]))/(1-zb) for q in (-1,1) for zb in (zlo,zhi)]
    angles=np.stack([np.arctan2(y,x) for y in ys for x in (xlo,xhi)])
    return (tlo,thi), (angles.min(axis=0),angles.max(axis=0))


def matrix(phases):
    return np.stack([np.ones_like(phases),np.cos(phases),np.sin(phases),np.cos(2*phases)],axis=-1)


def harmonic_check():
    # Use palindromic but nonconstant density layers at two energies
    # with independent incoherent positive weights on all final flavours.
    rng=np.random.default_rng(1669)
    weights=rng.random((2,2,3))
    h0=np.array([[2.,.8,.2],[.8,.3,1.],[.2,1.,-.2]])
    def count(delta,sign):
        s,c=.75,math.sqrt(1-.75**2)
        rot=np.array([[1,0,0],[0,c,s],[0,-s,c]])@np.diag([1,1,np.exp(1j*sign*delta)])
        total=0.
        for i,energy in enumerate((.325,.647)):
            evol=np.eye(3,dtype=complex)
            for potential in (.07,.16,.07):
                h=(rot@h0@rot.conj().T)/energy + sign*np.diag([potential,0,0])
                vals,vecs=np.linalg.eigh(h)
                evol=((vecs*np.exp(-.2j*vals))@vecs.conj().T)@evol
            total+=np.sum(weights[i]*(abs(evol[:,:2].T)**2))
        return total
    nodes=np.pi*np.array([-.5,0,.5,1.])
    probes=np.linspace(-np.pi,np.pi,37)
    error=0.
    for sign in (-1,1):
        coefs=np.linalg.solve(matrix(nodes),[count(x,sign) for x in nodes])
        error=max(error,max(abs(matrix(probes)@coefs-[count(x,sign) for x in probes])))
    assert error < 3e-13
    return error


def solve():
    old=read('../1668/results.json')['branches']
    nodes=np.pi*np.array([-.5,0,.5,1.])
    counts=np.array([[96.55,81.59,66.89,81.85],[16.56,18.81,20.75,18.49],[9.30,8.10,6.59,7.79]])
    response=np.linalg.solve(matrix(nodes),counts.T)
    result={'verified_harmonic_residual_below':1e-12,'observed':[94,16,14],
            'table9_coefficients':response.T.round(6).tolist(),'branches':{}}
    harmonic_check()
    for branch,oldrow in old.items():
        a,b=oldrow['with_JUNO_two_widths']
        edges=np.linspace(a-2e-9,b+2e-9,1+math.ceil((b-a+4e-9)/5e-6))
        ts,ds=coordinates(branch,.5*(edges[1:]+edges[:-1]),max(np.diff(edges))/2)
        witness=oldrow['witness']
        t0,d0=coordinates(branch,[witness['R']])
        phase=witness['delta_over_pi']*np.pi
        assert d0[0][0] < phase < d0[1][0] and t0[0][0]<witness['three_angles'][2]<t0[1][0]
        delta=np.array([ds[0].min(),ds[1].max()])
        theta=np.array([ts[0].min(),ts[1].max()])
        samples=matrix(phase)@response
        assert np.all(samples>0)
        if branch.startswith('second'):
            assert theta[0]>.5985 and delta[0]>0 and delta[1]<3.095
        else:
            assert theta[0]>.4365 and theta[1]<.5975
            assert -.285<phase<-.275
        result['branches'][branch]={'sin2theta23_enclosure':theta.round(10).tolist(),
            'delta_rad_enclosure':delta.round(10).tolist(),'witness_delta_rad':round(float(phase),9),
            'fixed_nuisance_fixed_other_osc_parameters_event_means':samples.round(6).tolist(),
            'table_rounding_error_per_mean':round(.005*np.sum(abs(matrix(phase)@np.linalg.inv(matrix(nodes)))),8)}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1669 independent PASS: projector enclosures, table harmonics, both T2K projections')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
