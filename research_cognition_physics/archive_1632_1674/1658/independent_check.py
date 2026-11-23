"""Constructive three-mass witnesses: fixed quark pulls and no electron calibration."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(rel):
    return json.loads((HERE/rel).read_bytes())


def reconstruct(vec):
    yd, ys, yb, theta12, theta13, theta23, delta = vec
    cp, cu, ct = np.cos([theta12, theta13, theta23])
    sp, su, st = np.sin([theta12, theta13, theta23])
    ph = np.exp(1j*delta)
    V = np.array([[cp*cu, sp*cu, su/ph],
                  [-sp*ct-cp*st*su*ph, cp*ct-sp*st*su*ph, st*cu],
                  [sp*st-cp*ct*su*ph, -cp*st-sp*ct*su*ph, ct*cu]])
    diag = np.array([yd/yb, ys/yb, 1])
    H = (V*diag**2) @ V.conj().T
    A = math.sqrt(H[0,0].real)
    C, E = np.abs(H[0,1:])/A
    # A^2 F^2 is a principal minor, given by a positive sum using unitarity.
    minor = sum(diag[i]**2*diag[j]**2*abs(V[1,k])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))
    F = math.sqrt(minor)/A
    B = -diag[0]*diag[1]/(A*F)
    Z = (H[1,2]*np.exp(1j*(np.angle(H[0,1])-np.angle(H[0,2])))-C*E)/F
    down = np.array([[0, A, 0], [B, C, Z], [0, E, F]])
    np.testing.assert_allclose(np.linalg.svd(down,compute_uv=False)[::-1],diag,atol=6e-15)
    return down


def calibrated(down, cg, relative_Z, target):
    ca, cb, cc = cg
    obs_e, mu, tau = target
    def shape(t):
        factors = np.array([[0,cb,0],[ca,cc,relative_Z*t],[0,t,t]])
        return np.linalg.svd(down*factors,compute_uv=False)

    def residual(t):
        trial = shape(t)
        return trial[0]/trial[1]-tau/mu

    lo, hi = .3, 4.
    assert residual(lo)<0<residual(hi)
    for _ in range(52):
        mid=(lo+hi)/2
        if residual(mid)>0:
            hi=mid
        else:
            lo=mid
    t = (lo+hi)/2
    shape_values = shape(t)
    gain = tau/shape_values[0]
    prediction = shape_values[::-1]*gain
    np.testing.assert_allclose(prediction[1:],target[1:],rtol=2e-13)
    return {'t':t,'gain_in_yb_units':gain,'predicted_ye':prediction[0],
            'ye_over_reference':prediction[0]/obs_e}


def evaluate():
    source = read('../1397/inputs.json')['high_scale']
    rows = source['down'] + source['angles_theta12_theta13_theta23_delta']
    leptons = read('../1402/inputs.json')['scales']['1e16'][-3:]
    target = np.array([v[0] for v in leptons])
    cgs = read('../1486/inputs.json')['clebsch']
    dirs = np.array([1,-1,-1,-1,1,-1,1])
    pulls = 1.992*dirs
    quark = [center + width*shift for (center,width),shift in zip(rows,pulls)]
    down = reconstruct(quark)
    witnesses = {}
    for name, cg in cgs.items():
        competitor = 1 if name=='GJ' else 3
        result = calibrated(down, cg, competitor, target)
        independent_e_pull = (result['predicted_ye']-target[0])/leptons[0][1]
        assert abs(independent_e_pull)<1.9
        if competitor!=1:
            # Three equal gains on E,F,Z and the explicitly split Z case differ.
            nominal = calibrated(down,cg,1,target)
            assert nominal['ye_over_reference'] < .8
            result['reference_ye_ratio_same_quarks_unsplit_Z']=nominal['ye_over_reference']
        result.update({'electron_pull':independent_e_pull,'relative_Z':competitor})
        witnesses[name]=result

    # A cross-check of both determinant-calibration branches at independent
    # box corners. Main interval proof covers interiors; this audit probes the
    # shape invariant directly using singular values, without that reconstruction.
    e0, mu0, tau0 = target
    tr = 1+(e0/tau0)**2+(mu0/tau0)**2
    volume = e0*mu0/tau0**2/tr**1.5
    observed = ((e0*mu0/tau0**2)**2+(e0/tau0)**2+(mu0/tau0)**2)/tr**2
    mins={}
    for name in ('small_left','large_left'):
        ca,cb,cc=cgs[name]
        gaps = []
        for signs in itertools.product((-1,1),repeat=7):
            mat=reconstruct([x+2*d*g for (x,d),g in zip(rows,signs)])
            block=mat*np.array([[0,cb,0],[ca,cc,0],[0,0,0]])
            rest=mat*np.array([[0,0,0],[0,0,1],[0,1,1]])
            K,L=(np.linalg.norm(m)**2 for m in (block,rest))
            target_product = volume*K*math.sqrt(L)/abs(ca*cb*np.linalg.det(mat))
            roots = np.roots([1.,-2.,1.,-target_product**2])
            valid = [w.real for w in roots if abs(w.imag)<1e-9 and 0<w.real<1]
            assert len(valid)==2
            for w in valid:
                for sign in (-1,1):
                    t=sign*math.sqrt(K*w/(L*(1-w)))
                    singular=np.linalg.svd(block+t*rest,compute_uv=False)
                    l=singular**2/np.dot(singular,singular)
                    np.testing.assert_allclose(np.sqrt(l.prod()),volume,rtol=1e-8)
                    gaps.append(abs((1-l@l)/2-observed))
        mins[name]=min(gaps)
    return {'quark_pulls':pulls.tolist(),'vectors_inside_the_common_box':witnesses,
            'corner_minimum_second_shape_gap_at_central_leptons':mins}


def tidy(o):
    if isinstance(o,dict): return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,list): return [tidy(i) for i in o]
    return float(f'{o:.10g}') if isinstance(o,(float,np.floating)) else o


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=tidy(evaluate())
    if args.check:
        assert read('independent_results.json')==result
        print('PASS: independent in-box witnesses, split-Z competitor, two signed determinant branches')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
