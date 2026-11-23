"""Reconstruct charged spectra by pivoted Cholesky and evolve the quark triangle."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def texture(data, delta):
    d = np.array(data[:3])/data[2]
    s12,s13,s23 = np.sin(data[3:6])
    c12,c13,c23 = np.cos(data[3:6])
    phase = np.exp(1j*delta)
    v = np.array([[c12*c13,s12*c13,s13*phase.conjugate()],
                  [-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13],
                  [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13]])
    gram = (v*d**2)@v.conj().T
    order = [0,2,1]
    lower = np.linalg.cholesky(gram[np.ix_(order,order)])
    matrix = lower[order,:][:,[2,0,1]]
    np.testing.assert_allclose(gram,matrix@matrix.conj().T,atol=1e-13)
    return matrix,v


def branches(data, leptons, delta):
    matrix,_ = texture(data,delta)
    light = np.zeros_like(matrix)
    light[:2,:2] = matrix[:2,:2]
    heavy = matrix-light
    light[1,1] *= 3
    k0 = (np.prod(leptons)/(data[0]*data[1]/data[2]**2))**(1/3)
    roots = np.roots([np.linalg.norm(heavy)**2,0,-np.dot(leptons,leptons)/k0**2,np.linalg.norm(light)**2])
    positive = sorted(w.real for w in roots if abs(w.imag)<1e-12 and w.real>0)
    assert len(positive)==2
    rows = []
    for w in positive:
        for sign in (1,-1):
            t = sign*w**1.5
            k = k0/np.sqrt(w)
            sv = np.linalg.svd(k*(light+t*heavy),compute_uv=False)[::-1]
            rows.append((t,k,sv))
    return rows


def gamma(v):
    return np.degrees(np.angle(-v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate()))) % 180


def low_energy(data, leptons, delta, count):
    # One-loop SM in the convention that the doublet index is on the left.
    # Input gauge couplings are evolved analytically up from M_Z in the same scheme.
    low = math.log(91.1876)
    high = math.log(1e16)
    b = np.array([41/10,-19/6,-7])
    norm = 16*np.pi**2
    reference = read('../1369/inputs.json')
    g0 = np.array([reference[k][0] for k in ('g1_GUT','g2','g3')])
    up = read('../1397/inputs.json')['high_scale']['up']
    _,v = texture(data,delta)
    y = np.stack([np.diag([row[0] for row in up]),v@np.diag(data[:3]),np.diag(leptons)])
    def rhs(t,y):
        u,d,e = [a@a.conj().T for a in y]
        g1,g2,g3 = 1/(1/g0**2-2*b*(t-low)/norm)
        trace = (3*np.trace(u+d)+np.trace(e)).real
        universal = np.array([trace-17*g1/20-9*g2/4-8*g3,
                              trace-g1/4-9*g2/4-8*g3,
                              trace-9*g1/4-9*g2/4])
        return (1.5*np.stack([u-d,d-u,e])@y+universal[:,None,None]*y)/norm
    h = (low-high)/count
    for n in range(count):
        t = high+n*h
        k1 = rhs(t,y)
        k2 = rhs(t+h/2,y+h*k1/2)
        k3 = rhs(t+h/2,y+h*k2/2)
        k4 = rhs(t+h,y+h*k3)
        y += h*(k1+2*k2+2*k3+k4)/6
    # Direct SVD, with ascending mass columns. Arbitrary eigenvector phases
    # cancel out of the rephasing-invariant angle gamma.
    basis = [np.linalg.svd(a)[0][:,::-1] for a in y[:2]]
    return gamma(basis[0].conj().T@basis[1])


def solve():
    raw = read('../1397/inputs.json')['high_scale']
    bounds = raw['down']+raw['angles_theta12_theta13_theta23_delta'][:3]
    charged = read('../1402/inputs.json')['scales']['1e16'][3:]
    central = np.array(bounds)[:,0]
    l0 = np.array(charged)[:,0]
    a,b = 1.2,1.4
    for _ in range(48):
        mid = (a+b)/2
        if branches(central,l0,mid)[2][2][1]>l0[1]:
            a=mid
        else:
            b=mid
    root = (a+b)/2
    windows = read('results.json')
    configs = {'central':(central,l0,root,2)}
    for name,slot in [('two_width_CP_conserving',3),('two_width_overlapping_direct_phase',2)]:
        row = windows[name]
        configs[name]=(row['quark_yukawas_and_angles'],row['three_charged_yukawas'],np.arccos(row['cos_delta']),slot)
    out = {}
    for name,(data,leptons,delta,slot) in configs.items():
        t,k,spectrum = branches(data,leptons,delta)[slot]
        np.testing.assert_allclose(spectrum,leptons,rtol=2e-9)
        value = gamma(texture(data,delta)[1])
        at_z = low_energy(data,leptons,delta,96)
        fine = low_energy(data,leptons,delta,192)
        difference = lambda x,y: (x-y+90)%180-90
        assert abs(difference(at_z,fine))<1e-6
        cg = 0. if abs(np.sin(delta))<1e-12 else value
        weak = 0. if cg==0 else fine
        out[name]={'cos_delta':float(np.cos(delta)), 't':float(t),'normalization':float(k),
                   'gamma_at_matching_deg':float(cg),'gamma_at_MZ_one_loop_deg':float(weak),
                   'MZ_minus_high_in_deg_mod180':float(difference(weak,cg))}
    assert out['central']['gamma_at_MZ_one_loop_deg']>67.8
    assert 58.8<out['two_width_overlapping_direct_phase']['gamma_at_MZ_one_loop_deg']<67.8
    return {case:{key:round(value,9) for key,value in row.items()} for case,row in out.items()}


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    out = solve()
    if args.check:
        assert out==read('independent_results.json')
        print('PASS: permuted Cholesky, three charged spectra and one-loop SM angle transport')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as stream:
            stream.write(json.dumps(out,indent=2)+'\n')
