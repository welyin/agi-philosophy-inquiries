"""Reconstruct loop-phase harmonics from spectral projectors and solve their inverse."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def box(row,n):
    return (row[0]-n*row[1],row[0]+n*row[-1])


def parameters(data,name):
    if name=='q':
        y=np.square([a[0] for a in data['yukawa']['down']])
        row=(y-y[0])/(y[2]-y[0])
        a12,a13,a23=[x[0] for x in data['quark_sines_and_phase'][:3]]
    else:
        vals=data[name]
        d=vals['solar_gap'][0]
        atm=vals['atmospheric_gap'][0]
        row=[0,d/(atm+d),1] if name=='NO' else [0,d/atm,-1]
        a12=math.sqrt(vals['solar_sine_squared'][0])
        a23=math.sqrt(vals['atmospheric_sine_squared'][0])
        a13=math.sin(math.asin(math.sqrt(data['DayaBay_amplitude'][0]))/2)
    return (a12,a13,a23),np.array(row)


def matrix(parameters,phase):
    (a,b,c),eig=parameters
    A,B,C=np.sqrt(1-np.square([a,b,c]))
    z=cmath.exp(1j*phase)
    U=np.array([[A*B,a*B,b/z],[-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])
    return (U*eig)@U.conj().T


def loop(h):
    return h[0,1]*h[1,2]*h[2,0]


def harmonics(param):
    p,q,r=[loop(matrix(param,t)) for t in (0.,math.pi/2,math.pi)]
    a=q.real
    b=(p.real-r.real)/2
    c=(p.real+r.real)/2-a
    d=q.imag
    scale=max(abs(x) for x in (a,b,c,d))
    for t in np.linspace(-math.pi,math.pi,81):
        model=a+b*math.cos(t)+c*math.cos(t)**2+1j*d*math.sin(t)
        assert abs(model-loop(matrix(param,t)))<scale*2e-12
    return np.array([a,b,c,d])/scale


def distance_real(z):
    return math.atan2(abs(z.imag),abs(z.real))


def extrema(coeff,angular_range):
    a,b,c,d=coeff
    xs=[angular_range[0],angular_range[1]]
    for poly in ([-c,0,a+2*c,b],[c,b,a]):
        for v in np.roots(poly):
            if abs(v.imag)<1e-8 and angular_range[0]<v.real<angular_range[1]:
                xs.append(v.real)
    result=[distance_real(complex(a+b*x+c*x*x,d*math.sqrt(max(0,1-x*x)))) for x in xs]
    return min(result),max(result)


def inverse(coeff,phase):
    a,b,c,d=coeff
    s,h=math.sin(phase),math.cos(phase)
    poly=[-s*(a-b+c),2*d*h,-s*(2*a-2*c),2*d*h,-s*(a+b+c)]
    result=[]
    for root in np.roots(poly):
        if abs(root.imag)>1e-7:
            continue
        t=2*math.atan(root.real)
        x=math.cos(t)
        z=complex(a+b*x+c*x*x,d*math.sin(t))
        if (z*cmath.exp(-1j*phase)).real>0:
            result.append(t)
    # delta=pi is missing in the tangent half-angle chart.
    if abs(cmath.phase(complex(a-b+c)*cmath.exp(-1j*phase)))<1e-10:
        result.append(math.pi)
    return result


def evaluate():
    cfg=json.loads((HERE/'../1386/inputs.json').read_bytes())
    p={k:parameters(cfg,k) for k in ['q','NO','IO']}
    h={k:harmonics(pp) for k,pp in p.items()}
    ranges={k:extrema(h[k],[-1,1]) for k in ['NO','IO']}
    # Enlarge the two measured tree-angle branches as in the admission.
    q_ranges=[]
    for a,b in [(56,82),(236,262)]:
        ends=sorted([math.cos(math.radians(a)),math.cos(math.radians(b))])
        q_ranges.append(extrema(h['q'],ends))
    gap=min(lo for lo,hi in q_ranges)-max(hi for lo,hi in ranges.values())
    assert gap>math.radians(3)
    phase=cmath.phase(loop(matrix(p['q'],math.radians(68.7))))
    counts={}
    for branch in ['NO','IO']:
        counts[branch]=[len(inverse(h[branch],t)) for t in [phase,-phase,phase+math.pi,-phase+math.pi]]
        assert max(counts[branch])==0
        probe=.63
        allowed=cmath.phase(loop(matrix(p[branch],probe)))
        assert min(abs(cmath.phase(cmath.exp(1j*(v-probe)))) for v in inverse(h[branch],allowed))<1e-8
    # Sufficiency check for a nonzero cyclic triple, with arbitrary shifts on its diagonal.
    first=matrix(p['q'],.73)
    gain=np.array([.3,1.9,7.])*np.exp(1j*np.array([.2,-.85,.9]))
    second=np.diag(gain)@first@np.diag(gain.conj())+np.diag([.08,.2,-.03])
    edges=[(0,1),(1,2),(2,0)]
    ratios=[second[i,j]/first[i,j] for i,j in edges]
    u,v,w=map(abs,ratios)
    mag=np.sqrt([u*w/v,u*v/w,v*w/u])
    angle=[0.,-cmath.phase(ratios[0]),cmath.phase(ratios[2])]
    recovered=mag*np.exp(1j*np.array(angle))
    residual=second-np.diag(recovered)@first@np.diag(recovered.conj())
    assert np.max(abs(residual-np.diag(np.diag(residual))))<1e-12
    return {'at_central_masses_and_angles':{'quark_branches_deg':np.degrees(q_ranges).round(6).tolist(),
            'neutrino_all_delta_deg':{k:np.degrees(v).round(6).tolist() for k,v in ranges.items()},
            'minimum_gap_degrees': round(math.degrees(gap),6),
            'one_changed_edge_fraction_at_least':round(math.sin(gap),7),
            'three_changed_edges_common_cap_at_least':round(math.sin(gap/3),7)},
            'quartic_solution_counts_at_charged_tree_center':counts}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    target=HERE/'independent_results.json'
    data=evaluate()
    if args.check:
        assert data==json.loads(target.read_bytes())
        print('PASS: exact phase quartic, independent extrema and diagonal-congruence inverse')
    else:
        with target.open('x',encoding='utf-8') as out:
            out.write(json.dumps(data,indent=2)+'\n')
