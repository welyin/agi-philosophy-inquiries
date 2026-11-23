"""Flavor-resolved scalar decay certificate with free right-handed rotations."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def sharp(d, e):
    b, mu = d[2], e[1]
    others = (d[1]+e[2])**2 + (d[0]+e[0])**2
    assert others + mu**2 >= b*mu and b>mu
    return (b-mu)**2/((b-mu)**2+others)


def chiral(d, e, left, right):
    down, charged = np.diag(d),np.diag(e)
    g = down@right-left@charged
    h = down@left-right@charged
    np.testing.assert_allclose(abs(h)**2-abs(g)**2,
        (d[:,None]**2-e[None,:]**2)*(abs(left)**2-abs(right)**2),atol=1e-17)
    return g, h


def branch(h):
    return float(abs(h[2,1])**2/np.linalg.norm(h)**2)


def unitary(rng, n):
    m = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,r = np.linalg.qr(m)
    return q@np.diag(np.diag(r)/abs(np.diag(r)))


def physical(beta, phase):
    return beta*phase/(1-beta+beta*phase)


def calculate():
    cfg=json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    ds=np.array(json.loads((HERE/cfg['down_anchor']).read_text(encoding='utf-8'))['down_yukawa_dr_1TeV'])
    es=np.array(cfg['charged_leptons_DR_1TeV'])
    d, e = ds[:,0], es[:,0]
    dlo,dhi=d-2*ds[:,1],d+2*ds[:,1]
    elo,ehi=e-2*es[:,1],e+2*es[:,1]
    assert (dlo[1]+elo[2])**2 + (dlo[0]+elo[0])**2 + elo[1]**2 >= dhi[2]*ehi[1]
    lower = sharp(np.array([dhi[0],dhi[1],dlo[2]]),ehi)
    upper_angles=[]
    for signs in itertools.product([-2,2],repeat=6):
        dv,ev=d+np.array(signs[:3])*ds[:,1],e+np.array(signs[3:])*es[:,1]
        upper_angles.append(sharp(dv,ev))
    assert abs(lower-min(upper_angles))<1e-12
    swap=np.eye(3)[[0,2,1]]
    extremizer=np.diag([-1,-1,1])@swap
    g,h=chiral(d,e,swap,extremizer)
    assert abs(branch(h)-sharp(d,e))<1e-12
    rng=np.random.default_rng(1548)
    for _ in range(100):
        left=swap.astype(complex)
        left[:2,:]=unitary(rng,2)@left[:2,:]
        right=unitary(rng,3)
        _,h=chiral(d,e,left,right)
        assert branch(h)>=sharp(d,e)-1e-12
    # Right-handed direction can delete the OTHER pole's b mu vertex exactly.
    angle=math.acos(e[1]/d[2]); c,s=math.cos(angle),math.sin(angle)
    rot=np.array([[1,0,0],[0,c,-s],[0,s,c]])
    g,h=chiral(d,e,swap,rot@swap)
    assert abs(g[2,1])<1e-15
    zeros=chiral(d,e,np.eye(3),np.eye(3))
    assert all(a[2,1]==0 for a in zeros)
    m=cfg['mass_GeV']
    x,y=(5/m)**2,(2/m)**2
    phase=(1-x-y)*math.sqrt((1-x-y)**2-4*x*y)
    pref=np.linalg.norm(dhi)/(dlo[2]-ehi[1])
    def guaranteed(theta):
        perturb=2*math.sin(theta/2)*pref
        new_angle=math.acos(math.sqrt(lower))+math.asin(min(1.,perturb))
        return physical(math.cos(min(new_angle,math.pi/2))**2,phase)
    limits={str(deg):guaranteed(math.radians(deg)) for deg in cfg['misalignment_degrees']}
    btarget=cfg['cms']['branching_fraction']
    threshold=2*math.asin(math.sin(math.acos(math.sqrt(btarget/(phase+btarget*(1-phase))))-
                                 math.acos(math.sqrt(lower)))/(2*pref))
    return {'central_sharp_branch':sharp(d,e),'two_width_sharp_lower':lower,
            'on_shell_min_phase_factor':phase,'physical_bounds_at_angle_degrees':limits,
            'to_guarantee_half_max_angle_degrees':math.degrees(threshold),
            'min_left_probability_for_angle_certificate':math.cos(threshold)**2,
            'extra_width_relative_to_charged_max_for_half_at_zero_angle':limits['0']/btarget-1,
            'cms_min_visible_pair_response_at_zero_angle':(btarget/limits['0'])**2,
            'right_cancels_other_pole_coupling_bmu':float(abs(g[2,1])),
            'named_pole_branch_for_this_cancellation':branch(h)}


def compare(now,saved):
    if isinstance(now,dict):
        assert now.keys()==saved.keys()
        for k in now:
            compare(now[k],saved[k])
    else:
        np.testing.assert_allclose(now,saved,rtol=2e-10,atol=2e-12)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=calculate(); path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write('\n')
    else:
        compare(out,json.loads(path.read_text(encoding='utf-8')))
        print('1548: sharp flavor branch, unknown right direction, phase space and CMS transfer verified')
