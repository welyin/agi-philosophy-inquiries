"""Test neutral and charged phase transfers on identical full three-family matrices."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
CKM = runpy.run_path(str(HERE.parent/'1675/independent_check.py'))
NEUTRINO = runpy.run_path(str(HERE.parent/'1668/check.py'))
ICE = runpy.run_path(str(HERE.parent/'1670/check.py'))
PHI = -np.pi/6


def read(path):
    return json.loads((HERE/path).read_bytes())


def get_matrices(raw, gamma):
    a,b,c = raw[3:]
    ga = np.radians(gamma)
    delta = ga+np.arcsin(np.cos(a)*np.tan(c)*np.sin(b)*np.sin(ga)/np.sin(a))
    m,v = CKM['texture'](raw,delta)
    # Keep the reference weak basis used for the separately specified neutrino source.
    m = np.diag([1,np.exp(-1j*np.angle(m[1,1])),np.exp(-1j*np.angle(m[2,1]))])@m
    m = m@np.diag([-np.exp(-1j*np.angle(m[1,0])),1,np.exp(-1j*np.angle(m[2,2]))])
    short = m.copy()
    short[2,:] = 0
    short[:,2] = 0
    long = m-short
    short[1,1] *= 3
    return short,long,v,delta


def inverse_angles(raw,l,ga):
    short,long,_,delta = get_matrices(raw,ga)
    candidates = CKM['branches'](raw,l,delta)
    target = np.sum([l[i]**2*l[j]**2 for i,j in [(0,1),(0,2),(1,2)]])
    solutions = []
    for size,k,_ in candidates[::2]:
        def e2(theta):
            sv = np.linalg.svd(k*(short+size*np.exp(1j*theta)*long),compute_uv=False)
            return sum(sv[i]**2*sv[j]**2 for i,j in [(0,1),(0,2),(1,2)])
        q0,q90,q180 = [e2(th) for th in [0,np.pi/2,np.pi]]
        even = (q0+q180)/2
        b,c = q0-even,q90-even
        d = (target-even)/np.hypot(b,c)
        if abs(d)<=1:
            for theta in np.arctan2(c,b)+np.array([-1,1])*np.arccos(d):
                sv = np.linalg.svd(k*(short+size*np.exp(1j*theta)*long),compute_uv=False)[::-1]
                np.testing.assert_allclose(sv,l,rtol=1e-8)
                solutions.append(float((np.degrees(theta)+180)%360-180))
    return sorted(solutions)


def compute():
    cfg = read('../1397/inputs.json')['high_scale']
    qc,qe = np.array(cfg['down']+cfg['angles_theta12_theta13_theta23_delta'][:3]).T
    lc,le = np.array(read('../1402/inputs.json')['scales']['1e16'][3:]).T
    angles = {str(g):inverse_angles(qc,lc,g) for g in (63.2,-116.8)}
    pulls = np.array([-.9945598157425,.0483264007185,.6081167801086,.0643854883131,-.4143859681949,.5453269838802])
    raw = qc+qe*pulls
    short,long,v,delta = get_matrices(raw,63.2)
    mag = 1.17-.07*.2375804525632
    k = .0082+.0004*.260621711182
    def yuk(t):
        return k*(short+t*long).T
    monomial = mag*np.exp(1j*PHI)
    neutral = mag*(2+np.exp(4j*PHI))/np.sqrt(3)
    np.testing.assert_allclose(neutral,monomial,rtol=1e-14)
    lepton = yuk(neutral)
    e,sv,_ = np.linalg.svd(lepton)
    e,sv = e[:,::-1],sv[::-1]
    lp = (sv-lc)/le
    assert max(abs(lp))<1 and max(abs(pulls))<1
    R,ds = 5.83,7.388e-5
    phases = ICE['phase_legs']('third_leg_3_to_2')
    m,_,_ = NEUTRINO['OLD']['rank_two'](ds,ds*(R**2-1),2*np.pi/3,'small')
    P = np.diag(np.exp(1j*np.r_[0.,phases]))
    physical = P@m@P
    masses,modes = NEUTRINO['OLD']['takagi'](physical)
    u = e.conj().T@modes
    predicted = NEUTRINO['OLD']['angles'](u)
    windows = np.array(read('../1659/results.json')['joint_windows'])
    assert np.all((predicted>windows[:,0])&(predicted<windows[:,1]))
    d21,d32 = np.diff(masses**2)
    for value,key in [(d21,'juno_NO_dm21'),(d21+d32,'juno_NO_dm31')]:
        lo,hi = NEUTRINO['window'](NEUTRINO['INPUT'][key],1)
        assert lo<value<hi
    ice = ICE['interpolate'](*ICE['grid'](),predicted[2],d32)
    weakgamma = CKM['low_energy'](raw,sv,delta,192)
    # Angular spurion derivatives of the determinant ratio, in the fixed basis.
    def slopes(transfer):
        step=1e-5
        a,b = [np.linalg.det(yuk(transfer(PHI+sign*step))) for sign in (-1,1)]
        return np.log(b/a)/(2*step)
    slope1 = slopes(lambda t:mag*np.exp(1j*t))
    slope2 = slopes(lambda t:mag*(2+np.exp(4j*t))/np.sqrt(3))
    np.testing.assert_allclose(slope1,1j,atol=1e-8)
    np.testing.assert_allclose(slope2,4/np.sqrt(3),atol=1e-8)
    return {'continuous_t_phases_central_deg':angles,'joint_witness':{
        'quark_pulls':pulls.tolist(),'charged_lepton_pulls':lp.tolist(),'quark_gamma_high':CKM['gamma'](v),
        'quark_gamma_MZ_one_loop':weakgamma,'common_abs_t':mag,'k_yb':k,'R':R,'d21':d21,
        'd31':d21+d32,'d32':d32,'three_neutrino_angles':predicted.tolist(),
        'neutrino_delta_over_pi':NEUTRINO['phase_delta'](u),'icecube_conditional_q':float(ice)},
        'log_abs_and_phase_slopes':{'one_power':[slope1.real,slope1.imag],
                                    'neutral_two_paths':[slope2.real,slope2.imag]}}


def stable(x):
    if isinstance(x,dict): return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list): return [stable(a) for a in x]
    return round(float(x),12)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=stable(compute())
    if args.check:
        assert out==read('independent_results.json')
        print('PASS: independent three-family angular fit, neutral spurion and distinct response')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as stream:
            stream.write(json.dumps(out,indent=2)+'\n')
