"""Pairing of a common mass-basis left current and its electron detection test."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PERMS = list(itertools.permutations(range(3)))


def quark_boxes(cfg):
    q = json.loads((HERE/cfg['quark_anchor']).read_text(encoding='utf-8'))
    n = cfg['quark_widths']
    a, err = q['ud_beta']
    lo, hi = a - n*err, a + n*err
    f, *ferr = q['F_times_cb']
    c = (f + n*sum(ferr)) / q['F_lower']
    ratio, *rerr = q['ub_over_cb']
    b = c*(ratio+n*sum(rerr))
    t = c*c+b*b
    lower = [[lo,math.sqrt(1-hi*hi-b*b),0], [math.sqrt(1-hi*hi-t),math.sqrt(lo*lo-c*c),0],
             [0,0,math.sqrt(1-t)]]
    upper = [[hi,math.sqrt(1-lo*lo),b], [math.sqrt(1-lo*lo),math.sqrt(hi*hi+b*b),c],
             [math.sqrt(t),math.sqrt(t),1]]
    return np.array(lower), np.array(upper)


def conservative_cp_free():
    # A deliberately rounded outer box for the three angle ranges in NuFit Table 1;
    # all delta are allowed. Product bounds follow the triangle and reverse triangle.
    s12,c12 = [np.sqrt(a) for a in ([.27,.35],[.65,.73])]
    s23,c23 = [np.sqrt(a) for a in ([.38,.65],[.35,.62])]
    s13,c13 = [np.sqrt(a) for a in ([.020,.025],[.975,.980])]
    def two(a,b):
        return (max(0.,a[0]-b[1],b[0]-a[1]),min(1.,a[1]+b[1]))
    pairs = [[c12*c13,s12*c13,s13],
             [two(s12*c23,c12*s23*s13),two(c12*c23,s12*s23*s13),s23*c13],
             [two(s12*s23,c12*c23*s13),two(c12*s23,s12*c23*s13),c23*c13]]
    data = np.array(pairs)
    return data[:,:,0],data[:,:,1]


def bound(vlo, vhi, ulo, uhi):
    va, vb = np.arccos(vhi), np.arccos(vlo)
    ua, ub = np.arccos(uhi), np.arccos(ulo)
    candidates = []
    for pi in PERMS:
        for rho in PERMS:
            # pi(i)=neutrino index matched to up row i, rho(j)=lepton matched to down row j.
            a = ua[np.ix_(rho,pi)].T
            b = ub[np.ix_(rho,pi)].T
            gap = np.maximum(0,np.maximum(va-b,a-vb))
            candidates.append(float(gap.max()))
    return min(candidates)


def pmns(s12,s13,s23,delta):
    a,b,c = np.sqrt(1-np.array([s12,s13,s23])**2)
    p = np.exp(1j*delta)
    return np.array([[a*b,s12*b,s13/p],[-s12*c-a*s23*s13*p,a*c-s12*s23*s13*p,s23*b],
                     [s12*s23-a*c*s13*p,-a*s23-s12*c*s13*p,c*b]])


def leakage(matrix):
    return min(max(1-abs(matrix[i,p[i]])**2 for i in range(3)) for p in PERMS)


def main_results():
    cfg = json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    vlo,vhi = quark_boxes(cfg)
    ubox = cfg['pmns_union_of_IC19_IC24_modulus_projections']
    ulo,uhi = np.array(ubox['low']),np.array(ubox['high'])
    gap = bound(vlo,vhi,ulo,uhi)
    entries = {}
    for key,(low,high) in {'eq4':(ulo,uhi),'all_delta':conservative_cp_free()}.items():
        angle = bound(vlo,vhi,low,high)
        gap12 = low[0,0]**2-high[0,1]**2
        gap23 = low[0,1]**2-high[0,2]**2
        harmonic = 1/(1/gap12 + 1/gap23)
        angles = [0.,.01, math.sin(angle/2)**2]
        data = []
        for d in angles:
            u = math.sin(max(0,angle-math.asin(math.sqrt(d))))**2
            data.append({'max_down_loss':d,'min_up_loss':u,'min_far_electron_distance':u*harmonic})
        entries[key] = {'gap_radians':angle, 'gap_degrees':math.degrees(angle),
                        'at_least_one_side_max_missed':math.sin(angle/2)**2,
                        'electron_target_gaps':[gap12,gap23], 'conditional':data}
    ref=cfg['synthetic']
    s13=ref['Vcb']*ref['Vub_over_Vcb']
    c13=math.sqrt(1-s13**2)
    V=pmns(math.sqrt(1-(ref['Vud']/c13)**2),s13,ref['Vcb']/c13,math.radians(ref['quark_phase_deg']))
    U=pmns(*np.sqrt(ref['lepton_sin2']),math.radians(ref['lepton_phase_deg']))
    assert np.all((abs(V)>=vlo-1e-12) & (abs(V)<=vhi+1e-12))
    assert np.all((abs(U)>=ulo-1e-12) & (abs(U)<=uhi+1e-12))
    T=V.conj().T@U.conj().T
    eig, axes=np.linalg.eig(T)
    half=axes@np.diag(np.exp(.5j*np.angle(eig)))@np.linalg.inv(axes)
    p=abs(U[0])**2
    witnesses={}
    for name,kd in [('down_perfect',np.eye(3)),('split',half),('up_perfect',T)]:
        ku=V@kd@U
        np.testing.assert_allclose(kd.conj().T@kd,np.eye(3),atol=2e-14)
        np.testing.assert_allclose(ku@U.conj().T,V@kd,atol=2e-14)
        x=abs(ku)**2 @ p
        delta=min(float(np.max(np.abs(x-p[list(perm)]))) for perm in PERMS)
        ld,lu=max(0,leakage(kd)),max(0,leakage(ku))
        assert math.asin(math.sqrt(ld))+math.asin(math.sqrt(lu)) >= gap-1e-12
        assert delta + 1e-12 >= lu/(1/(p[0]-p[1])+1/(p[1]-p[2]))
        witnesses[name]={'down_loss':ld,'up_loss':lu,'electron_near':(abs(V@kd)[:,0]**2).tolist(),
                         'electron_far':x.tolist(),'far_distance':delta}
    return {'boxes':entries,'examples':witnesses,'synthetic_p':p.tolist()}


def compare(result, saved):
    if isinstance(result,dict):
        assert result.keys()==saved.keys()
        for key in result:
            compare(result[key],saved[key])
    elif isinstance(result,list) and result and isinstance(result[0],dict):
        assert len(result)==len(saved)
        for new,old in zip(result,saved):
            compare(new,old)
    else:
        np.testing.assert_allclose(result,saved,rtol=1e-9,atol=4e-13)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    opt=parser.parse_args()
    result=main_results()
    out=HERE/'results.json'
    if opt.write:
        with out.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
    else:
        compare(result,json.loads(out.read_text(encoding='utf-8')))
        print('1546: 36 permutations, CP widths, tagged near/far probes verified')
