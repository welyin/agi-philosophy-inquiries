"""Profile both real source strengths and the full positive Schur-source phase at fixed ports."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
PORT = runpy.run_path(str(ARCH/'1682/independent_check.py'))
CFG = json.loads((HERE/'inputs.json').read_bytes())


def roots(poly):
    return [float(x.real) for x in poly.roots() if abs(x.imag)<1e-7 and x.real>0]


def charged_sources(w):
    base = PORT['BASE']
    q = base['Q']
    raw = q[:,0]+q[:,1]*np.array(w['quark_marginal_pulls'])
    sh, long, _, _ = PORT['OLD']['get_matrices'](raw, CFG['quark_phase_degrees'])
    charged = (sh+w['mag_t']*np.exp(1j*PORT['PHI'])*long).T*w['charged_normalization']
    e = np.linalg.svd(charged)[0][:, ::-1]
    legs = PORT['UNIT'].copy()
    legs[2] *= w['neutral_modulus']
    return e.T @ (legs*base['PIV']['A']), e.T @ (legs*base['PIV']['B'])


def profile(w):
    a, b = charged_sources(w)
    aa, bb = np.vdot(a,a).real, np.vdot(b,b).real
    pair = np.vdot(a,b)
    assert abs(pair.imag)<2e-14 and pair.real>0
    g = pair.real
    gram = aa*bb-g*g
    R = w['R']
    t, gap = R+1/R, R-1/R
    span = np.column_stack([a.conj(),b.conj()])
    active = (span @ np.linalg.pinv(span)).diagonal().real
    axis = Poly([0.,1.])
    norm = 2*g*g*axis
    num = Poly([-aa*aa,gram*t,-bb*bb])
    discr = norm**2-num**2
    windows = PORT['BASE']['ANGLES']
    small = (1-np.sqrt(1-windows[0]))/2
    assert active[0]<1-small[1]  # Large-electron solution of the reactor observable excluded.
    # p_alpha3 = (P(r)+k sqrt(discr(r))) / r.
    polynomials, radical = [], []
    for i in range(3):
        overlap = np.conj(a[i])*b[i]
        polynomials.append((Poly([aa*abs(a[i])**2,0,bb*abs(b[i])**2])
                            +overlap.real/g*num - gram*active[i]/R*axis)/(gram*gap))
        radical.append(-overlap.imag/(g*gram*gap))
    e, mu = np.eye(3)[:2]
    # All cuts are in the form v . diag(P3) >= limit.
    sl, su = windows[1]
    al, au = windows[2]
    cuts = [(e,small[0]),(-e,-small[1]), (-(1-sl)*e,sl-active[0]),
            ((1-su)*e,active[0]-su), (mu+al*e,al),(-mu-au*e,-au)]

    def cos_phase(r):
        return float(num(r)/norm(r))

    def probabilities(r):
        rad = math.sqrt(max(0.,discr(r)))
        return np.array([(p(r)+k*rad)/r for p,k in zip(polynomials,radical)])

    def inside(r):
        if abs(cos_phase(r))>1+1e-12:
            return False
        p = probabilities(r)
        return all(v@p>=limit-2e-13 for v,limit in cuts)

    def eta_deg(r):
        return math.degrees(math.acos(np.clip(cos_phase(r),-1,1)))

    def delta(r):
        eta = math.radians(eta_deg(r))
        return 2+2*math.sqrt(3)/math.tan(eta)

    breaks = roots(num-norm)+roots(num+norm)
    for v, limit in cuts:
        p = sum((coef*f for coef,f in zip(v,polynomials)),Poly([0]))-limit*axis
        k = np.dot(v,radical)
        breaks.extend(roots(p*p-k*k*discr))
    breaks.sort()
    edges = []
    for r in breaks:
        if not edges or abs(r-edges[-1])>1e-8:
            edges.append(r)
    intervals = []
    for left,right in zip(edges,edges[1:]):
        if inside((left+right)/2):
            if intervals and abs(intervals[-1][1]-left)<1e-8:
                intervals[-1][1] = right
            else:
                intervals.append([left,right])
    for r in edges:
        if inside(r) and not any(l-1e-8<=r<=u+1e-8 for l,u in intervals):
            intervals.append([r,r])
    intervals.sort()
    ranges = []
    for l,u in intervals:
        points = [l,u]+([aa/bb] if l<aa/bb<u else [])
        assert all(-1<cos_phase(r)<1 for r in points)
        ranges.append(dict(r=np.round([l,u],10).tolist(),
                           phase_degrees=np.round([min(map(eta_deg,points)),max(map(eta_deg,points))],7).tolist(),
                           schur_over_v=np.round([min(map(delta,points)),max(map(delta,points))],8).tolist()))

    def at_shift(d):
        eta = np.angle((d-2)+2j*math.sqrt(3))
        ans = []
        for r in sorted(roots(num-math.cos(eta)*norm)):
            mass = np.outer(a,a)+r*np.exp(1j*eta)*np.outer(b,b)
            vals, U = np.linalg.eigh(mass.conj().T@mass)
            p = abs(U)**2
            np.testing.assert_allclose(math.sqrt(vals[-1]/vals[-2]), R, rtol=2e-12)
            np.testing.assert_allclose(probabilities(r),p[:,2],atol=4e-13)
            pe = p[0,2]
            obs = np.array([4*pe*(1-pe),p[0,1]/(1-pe),p[1,2]/(1-pe)])
            accept = bool(np.all((windows[:,0]<=obs)&(obs<=windows[:,1])))
            assert accept==inside(r)
            ans.append(dict(r=round(r,10),mixing=np.round(obs,9).tolist(),inside=accept))
        return ans

    tests = {str(d):at_shift(d) for d in CFG['schur_examples']}
    for i, value in enumerate(w['three_angles']):
        np.testing.assert_allclose(tests['0'][0]['mixing'][i],value,atol=1e-8)
    assert tests['0'][0]['inside']
    return dict(frozen_R=R, allowed_regions=ranges, profiles=tests)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    arg = p.parse_args()
    points = json.loads((HERE/CFG['charged_and_neutral_witnesses']).read_bytes())['witnesses']
    result = [profile(w) for w in points]
    dest = HERE/'results.json'
    if arg.check:
        assert result == json.loads(dest.read_bytes())
        print('PASS: polynomial boundary elimination, both source roots, direct three-flavor diagonalization')
    else:
        with dest.open('x',encoding='utf8') as fp:
            fp.write(json.dumps(result,indent=2)+'\n')
