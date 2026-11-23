"""Test both flavour mechanisms on the two-dimensional IceCube normal-order map."""
import argparse
import hashlib
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = runpy.run_path(str(HERE.parent/'1668/check.py'))
PHASE = MODEL['PHASE']
OLD = MODEL['OLD']


def read(path):
    return json.loads((HERE/path).read_bytes())


def grid():
    data=read('inputs.json')
    old=read(data['map_reference'])['icecube']
    x=np.array(old['grid_x']+data['additional_x'])
    m=np.array(old['mass_squared_eV2'])
    q=np.array(old['q_rows']+data['additional_q_rows'])
    assert q.shape==tuple(old['source_shape'])==(len(x),len(m))
    return x,m,q


def verify_public_map(path):
    saved = read('../1555/inputs.json')['icecube']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == saved['sha256']
    raw = np.load(path, allow_pickle=True).item()
    xx, mm, qq = (raw[k] for k in ('sin2_theta23','deltam232','deltallh'))
    x, m, q = grid()
    np.testing.assert_array_equal(x, np.unique(xx))
    np.testing.assert_array_equal(m, np.unique(mm))
    np.testing.assert_array_equal(q, qq[np.lexsort((mm,xx))].reshape(q.shape))


def interpolate(x,m,table,xx,mm):
    xx,mm=np.asarray(xx),np.asarray(mm)
    assert xx.min()>=x[0] and xx.max()<=x[-1] and mm.min()>=m[0] and mm.max()<=m[-1]
    i=np.clip(x.searchsorted(xx)-1,0,len(x)-2)
    j=np.clip(m.searchsorted(mm)-1,0,len(m)-2)
    u=(xx-x[i])/(x[i+1]-x[i])
    v=(mm-m[j])/(m[j+1]-m[j])
    return (1-u)*((1-v)*table[i,j]+v*table[i,j+1])+u*((1-v)*table[i+1,j]+v*table[i+1,j+1])


def rectangle_min(axes,xlo,xhi,mlo,mhi):
    x,m,q=axes
    # Cross all grid lines clipped to each rectangle, keeping array memory
    # bounded by the few grid intervals intersected by this family of boxes.
    xs=[xlo,xhi]+[np.clip(v,xlo,xhi) for v in x if xlo.min()<v<xhi.max()]
    ms=[mlo,mhi]+[np.clip(v,mlo,mhi) for v in m if mlo.min()<v<mhi.max()]
    best=np.full(len(xlo),np.inf)
    for xx,mm in itertools.product(xs,ms):
        best=np.minimum(best,interpolate(x,m,q,xx,mm))
    return best


def upper_x_at_level(axes,level):
    x,m,q=axes
    candidates=[v for v,row in zip(x,q) if min(row)<=level]
    for i in range(len(x)-1):
        for a,b in zip(q[i],q[i+1]):
            if (a-level)*(b-level)<0:
                candidates.append(x[i]+(level-a)*(x[i+1]-x[i])/(b-a))
    return max(candidates)


def phase_legs(name):
    phi=-math.pi/6
    if name=='third_leg_3_to_2':
        return np.array([2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))])
    return np.array([10*phi+math.radians(11.5),7*phi])


def branch(name, sign, bounds, charged, angles, axes, ds, da, count):
    grid_r=np.linspace(*bounds,count+1)
    R=(grid_r[1:]+grid_r[:-1])/2
    h=(bounds[1]-bounds[0])/(2*count)
    coeff=PHASE['coefficients'](charged,phase_legs(name))
    r=PHASE['root'](R)
    eta=np.exp(sign*2j*math.pi/3)
    pe=PHASE['projector'](coeff,eta,R,r,0,0).real
    pm=PHASE['projector'](coeff,eta,R,r,1,1).real
    p1=abs(coeff[2][0])**2
    x=pm/(1-pe)
    dslo=np.maximum(ds[0],da[0]/(R+h)**2)
    dshi=np.minimum(ds[1],da[1]/(R-h)**2)
    pp=2*h+2e-8
    xlo=(pm-pp)/(1-pe+pp)
    xhi=(pm+pp)/(1-pe-pp)
    mlo=((R-h)**2-1)*dslo
    mhi=((R+h)**2-1)*dshi
    assert np.all(mlo<=mhi)
    vals=rectangle_min(axes,xlo,xhi,mlo,mhi)

    angle_values=np.stack([4*pe*(1-pe),1-p1/(1-pe),x])
    valid=np.all((angle_values>angles[:,0,None]+1e-10) & (angle_values<angles[:,1,None]-1e-10),axis=0)
    # Construct actual admissible points using the same R, d21 and d31.
    lo=np.maximum(ds[0],da[0]/R**2)*(R**2-1)
    hi=np.minimum(ds[1],da[1]/R**2)*(R**2-1)
    valid &= lo<=hi
    best=(math.inf,None)
    for knot in axes[1]:
        mm=np.clip(knot,lo,hi)
        value=np.where(valid,interpolate(*axes,x,mm),np.inf)
        j=int(value.argmin())
        if value[j]<best[0]:
            best=(value[j],{'R':R[j],'d21':mm[j]/(R[j]**2-1),'d32':mm[j],'sin2theta23':x[j]})
    assert np.isfinite(best[0])
    current=read('../1668/results.json')['branches'][name]['witness']
    old_readout=interpolate(*axes,current['three_angles'][2],current['d32'])
    return {'allowed_x_outer':[min(xlo),max(xhi)],'allowed_d32_outer':[min(mlo),max(mhi)],
            'infimum_bracket':[min(vals)-1e-8,best[0]+1e-8],
            'feasible_near_minimum':best[1],'round1668_witness_q':old_readout}


def compute():
    inputs=read('inputs.json')
    cfg=read(inputs['solar_and_atmospheric'])
    old=read(inputs['previous_fit'])
    ds=MODEL['window'](cfg['juno_NO_dm21'],2)
    da=MODEL['window'](cfg['juno_NO_dm31'],2)
    axes=grid()
    charged,_=OLD['charged_basis'](read(cfg['charged']))
    angles=np.array(read(cfg['angle_windows'])['joint_windows'])
    x_edge=upper_x_at_level(axes,inputs['working_q'])
    results={}
    for name,sign,_,_ in cfg['branches']:
        angle_range=MODEL['allowed_r'](charged,phase_legs(name),sign,angles)
        juno=MODEL['gap_box'](2)
        bounds=[max(angle_range[0],juno[0]),min(angle_range[1],juno[1])]
        np.testing.assert_allclose(bounds,old['branches'][name]['with_JUNO_two_widths'],rtol=1e-10)
        results[name]=branch(name,sign,bounds,charged,angles,axes,ds,da,inputs['curve_cells'])
    third,second=results.values()
    assert third['infimum_bracket'][1]<1 and second['infimum_bracket'][0]>4
    return {'branches':results, 'upper_x_at_q4_profiling_all_masses':x_edge,
            'necessary_second_branch_downward_x_shift':second['allowed_x_outer'][0]-x_edge,
            'necessary_second_branch_q_reduction':second['infimum_bracket'][0]-inputs['working_q'],
            'separation_in_best_q_at_least':second['infimum_bracket'][0]-third['infimum_bracket'][1]}


def tidy(data):
    if isinstance(data,dict):
        return {k:tidy(v) for k,v in data.items()}
    if isinstance(data,(tuple,list)):
        return [tidy(v) for v in data]
    return float(f'{data:.11g}') if isinstance(data,(float,np.floating)) else data


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--verify-public-map',type=Path)
    args=parser.parse_args()
    if args.verify_public_map:
        verify_public_map(args.verify_public_map)
    out=tidy(compute())
    path=HERE/'results.json'
    if args.check:
        assert out==read('results.json')
        print('PASS: IceCube 2D grid, continuous source curves, extremum brackets and shift bounds')
    else:
        with path.open('x',encoding='utf8') as stream:
            stream.write(json.dumps(out,indent=2)+'\n')
