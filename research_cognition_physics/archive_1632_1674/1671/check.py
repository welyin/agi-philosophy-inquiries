"""Pull back Super-K's four-dimensional atmospheric table to two flavour mechanisms."""
import argparse
import hashlib
import io
import itertools
import json
import runpy
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PREV = runpy.run_path(str(HERE.parent / '1670/check.py'))
MODEL = PREV['MODEL']
PHASE = PREV['PHASE']


def read(path):
    return json.loads((HERE / path).read_bytes())


def section(raw):
    return [np.array(nodes) for nodes in raw['axes']], np.array(raw['q'])


def verify_release(path, data):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == data['archive_sha256']
    with tarfile.open(path) as archive:
        content = archive.extractfile(data['member']).read()
    assert hashlib.sha256(content).hexdigest() == data['member_sha256']
    original = np.loadtxt(io.BytesIO(content))
    assert original.shape == (141120, 5) and original[:, 4].min() == 0
    for raw in data['sections'].values():
        axes, values = section(raw)
        mask = np.logical_and.reduce([np.isin(original[:, i], axis) for i, axis in enumerate(axes)])
        rows = original[mask]
        assert len(rows) == values.size
        keys = tuple(axis.searchsorted(rows[:, i]) for i, axis in enumerate(axes))
        np.testing.assert_array_equal(values[keys], rows[:, 4])


def interp(axes, q, coordinates):
    loc, frac = [], []
    for axis, coord in zip(axes, coordinates):
        assert np.all(coord >= axis[0]) and np.all(coord <= axis[-1])
        i = np.clip(np.searchsorted(axis, coord) - 1, 0, len(axis)-2)
        loc.append(i)
        frac.append((coord-axis[i])/(axis[i+1]-axis[i]))
    result = 0.
    for bits in itertools.product((0, 1), repeat=len(axes)):
        key = tuple(i + bit for i, bit in zip(loc, bits))
        weight = np.prod(np.broadcast_arrays(*[u if bit else 1-u for u, bit in zip(frac, bits)]), axis=0)
        result = result + weight*q[key]
    return result


def clipped_knots(axis, lo, hi):
    return [lo, hi] + [np.clip(k, lo, hi) for k in axis if lo.min() < k < hi.max()]


def evaluate(name, sign, section_data, count, charged, windows, cfg):
    axes, q = section(section_data)
    angles = PREV['phase_legs'](name)
    rrange = MODEL['allowed_r'](charged, angles, sign, windows)
    jb = MODEL['gap_box'](2)
    edges = np.linspace(max(jb[0],rrange[0]), min(jb[1],rrange[1]), count+1)
    rlo, rhi = edges[:-1], edges[1:]
    R = (rlo+rhi)/2
    h = np.diff(edges).max()/2
    coeff = PHASE['coefficients'](charged, angles)
    eta = np.exp(sign * 2j*np.pi/3)
    roots = PHASE['root'](R)
    def p3(i, j):
        return PHASE['projector'](coeff, eta, R, roots, i, j)
    v = coeff[2]
    z = p3(0,0).real
    mu = p3(1,1).real
    x = mu/(1-z)
    s = 1-abs(v[0])**2/(1-z)
    J = -(v[0]*v[1].conjugate()*p3(1,0)).imag
    real = abs(v[1])**2 - s*(1-x) - (1-s)*x*z
    imag = 2*J/(1-z)
    delta = np.angle(real+1j*imag) % (2*np.pi)
    assert max(z)+2*h < .04 and min(real)-16*h > .059
    direction_pad = 400*h + 3e-6
    assert min(delta)>direction_pad and max(delta)<2*np.pi-direction_pad
    ds = MODEL['window'](cfg['juno_NO_dm21'], 2)
    da = MODEL['window'](cfg['juno_NO_dm31'], 2)
    mlo = (rlo**2-1)*np.maximum(ds[0], da[0]/rhi**2)
    mhi = (rhi**2-1)*np.minimum(ds[1], da[1]/rlo**2)
    pad = 2*h+2e-8
    ranges = [(mlo,mhi), ((mu-pad)/(1-z+pad),(mu+pad)/(1-z-pad)),
              (delta-direction_pad,delta+direction_pad), (z-pad,z+pad)]
    # Every multilinear cell takes its minimum on a corner, also after clipping.
    vertices = [clipped_knots(axis, lo, hi) for axis, (lo,hi) in zip(axes, ranges)]
    lower = min(float(interp(axes,q,coords).min()) for coords in itertools.product(*vertices))
    angle_vector = np.stack([4*z*(1-z), s, x])
    feasible = np.all((angle_vector>windows[:,0,None]+1e-10) & (angle_vector<windows[:,1,None]-1e-10),axis=0)
    dlo = np.maximum(ds[0], da[0]/R**2)
    dhi = np.minimum(ds[1], da[1]/R**2)
    feasible &= dlo<=dhi
    best = float('inf')
    for k in axes[0]:
        m = np.clip(k, (R**2-1)*dlo, (R**2-1)*dhi)
        values = np.where(feasible,interp(axes,q,(m,x,delta,z)),np.inf)
        j = int(values.argmin())
        if values[j] < best:
            best = float(values[j])
            witness = {'R':R[j], 'd21':m[j]/(R[j]**2-1),'m':m[j],'x':x[j],'delta':delta[j],'z':z[j]}
    assert np.isfinite(best)
    # Free just the CP coordinate at that same physical point. The profile can
    # be found exactly at a published phase node for piecewise linear interpolation.
    w = witness
    scores = [interp(axes,q,(w['m'],w['x'],angle,w['z'])) for angle in axes[2]]
    i = int(np.argmin(scores))
    return {'infimum_bracket':[lower-1e-7,best+1e-7],'near_minimum_witness':w,
            'free_phase_at_same_witness':{'q':scores[i],'delta':axes[2][i],'q_reduction':best-scores[i]},
            'outer_coordinates':[[lo.min(),hi.max()] for lo,hi in ranges]}


def compute(data):
    cfg = read('../1668/inputs.json')
    E, _ = PREV['OLD']['charged_basis'](PREV['read'](cfg['charged']))
    windows = np.array(PREV['read'](cfg['angle_windows'])['joint_windows'])
    out = {name:evaluate(name,sign,data['sections'][name],data['cells'],E,windows,cfg)
           for name,sign,*_ in cfg['branches']}
    first, second = [result['infimum_bracket'] for result in out.values()]
    # Freeze the counterexample registered by the IceCube review, before
    # consulting its SK readout. Both source ports give legitimate x,z coordinates.
    phi = -np.pi/6
    R, d = 5.8, 7.35e-5
    angle = np.angle(2-np.exp(4j*phi))
    u = MODEL['neutrino_mix'](E,[10*phi+angle,7*phi],-1,R)
    delta = MODEL['phase_delta'](u)*np.pi % (2*np.pi)
    z = abs(u[0,2])**2
    x = abs(u[1,2])**2/(1-z)
    m = (R**2-1)*d
    axes,q = section(data['sections']['third_leg_3_to_2'])
    check_point = {'R':R,'d21':d,'m':m,'x':x,'delta':delta,'z':z,
                   'q_SK':interp(axes,q,(m,x,delta,z))}
    return {'branches':out,'second_minus_third_infimum_interval':[second[0]-first[1],second[1]-first[0]],
            'second_leg_weights_2_and_minus1_fixed_witness':check_point}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--verify-release',type=Path)
    args = parser.parse_args()
    data = read('inputs.json')
    if args.verify_release:
        verify_release(args.verify_release,data)
    out = PREV['tidy'](compute(data))
    if args.check:
        assert out==read('results.json')
        print('PASS: SK four-coordinate source curves and conditional phase profiles')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
