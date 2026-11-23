"""Three-flavour appearance enclosures allowing an arbitrary bounded Earth profile."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'1659'/'independent_check.py'
HBARC = 1.973269804e-7  # eV m
VE = math.sqrt(2)*1.1663787e-23*6.02214076e23*(100*HBARC)**3
BASE = 295e3/HBARC


def read(name):
    return json.loads((HERE/name).read_bytes())


def physical_port(charged, direction, ratio, gap):
    a = np.array([0,1,1])
    b = np.array([1,3,1])
    d = 6*(ratio+1/ratio)+16
    r = 8/(d+math.sqrt(d*d-44**2))
    phi = -math.pi/6
    if direction=='third':
        phase = [0,2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))]
        sign = 1
    else:
        phase = [0,10*phi+math.radians(11.5),7*phi]
        sign = -1
    mass = np.outer(a,a)+r*np.exp(sign*2j*math.pi/3)*np.outer(b,b)
    phases=np.exp(1j*np.array(phase))
    mass=phases[:,None]*mass*phases
    _,sv,vh=np.linalg.svd(mass)
    mixing=charged.conj().T@vh.conj().T[:,::-1]
    squared=gap*(sv[::-1]/sv[1])**2
    np.testing.assert_allclose(squared,[0,gap,gap*ratio**2],atol=1e-16)
    weights=abs(mixing)**2
    s13=weights[0,2]
    s12,s23=weights[0,1]/(1-s13),weights[1,2]/(1-s13)
    j=(mixing[0,0]*mixing[1,1]*np.conj(mixing[0,1]*mixing[1,0])).imag
    cc=math.sqrt(s12*(1-s12)*s13*s23*(1-s23))
    sine=j/((1-s13)*cc)
    cosine=(weights[1,0]-s12*(1-s23)-(1-s12)*s23*s13)/(2*cc)
    np.testing.assert_allclose(sine**2+cosine**2,1,atol=5e-14)
    delta_pi=math.atan2(sine,cosine)/math.pi
    return mixing, squared, {'R':ratio,'dm21':gap,'dm31':squared[2],
            'dm32':squared[2]-squared[1], 'A13':4*s13*(1-s13),
            's2theta12':s12,'s2theta23':s23,'J':j,
            'delta_pi':delta_pi, 'no_phase_in_three_sigma':
                any(-1.38 <= delta_pi+n <= .30 for n in (-2,0,2))}


def amplitude(mix,squared,q,y,sign,length=BASE):
    u = mix if sign==1 else mix.conj()
    m=(u*squared)@u.conj().T
    k=(length*.5e-9)*q[...,None,None]*m
    k[...,0,0]+=sign*length*VE*y
    vals,vecs=np.linalg.eigh(k)
    return np.sum(vecs[...,0,:]*vecs[...,1,:].conj()*np.exp(-1j*vals),axis=-1)


def screen(u,ms,lo,hi):
    # q is the reciprocal energy in GeV^-1. Centre each Hermitian perturbation.
    q=np.linspace(1/hi,1/lo,513)
    amplitude_rho=BASE*VE*.1/2
    amplitude_q=(q[1]-q[0])/2*BASE*.25e-9*(ms[-1]-ms[0])
    amps=np.array([abs(amplitude(u,ms,q,1.3,s)) for s in [1,-1]])
    def bounds(eps):
        floors=np.maximum(0,amps-eps)**2
        ceilings=np.minimum(1,amps+eps)**2
        return [min(floors[0]-ceilings[1]),max(ceilings[0]-floors[1])]
    # Also enclose all constant densities by partitioning density, retaining
    # only a small amplitude remainder in each cell.
    rows=[]
    for den in np.linspace(1.2,1.4,33):
        amps=np.array([abs(amplitude(u,ms,q,den,s)) for s in [1,-1]])
        rows.append(bounds(amplitude_q+amplitude_rho/32+1e-10))
    amps=np.array([abs(amplitude(u,ms,q,1.3,s)) for s in [1,-1]])
    broad=bounds(amplitude_q+amplitude_rho+1e-10)
    eps=amplitude_q+amplitude_rho+1e-10
    lows=np.maximum(0,amps.min(axis=1)-eps)**2
    highs=np.minimum(1,amps.max(axis=1)+eps)**2
    separate_weights=[lows[0]-highs[1],highs[0]-lows[1]]
    margin=max(0,separate_weights[0],-separate_weights[1])
    const=[min(x[0] for x in rows),max(x[1] for x in rows)]
    qmid=np.array(2/(lo+hi))
    straight=[abs(amplitude(u,ms,qmid,1.3,s))**2 for s in (1,-1)]
    # Verify several noncommuting density profiles using explicit propagation.
    layered=[]
    for profile in ([1.2,1.4,1.2],[1.4,1.4,1.2],[1.3]*3):
        pair=[]
        for sign in [1,-1]:
            mix=u if sign==1 else u.conj()
            states=np.eye(3,dtype=complex)
            for den in profile:
                m=(mix*ms)@mix.conj().T*qmid*.5e-9
                m[0,0]+=sign*VE*den
                d,v=np.linalg.eigh(m)
                states=((v*np.exp(-1j*d*BASE/3))@v.conj().T)@states
            pair.append(abs(states[0,1])**2)
        layered.append(pair[0]-pair[1])
    assert all(broad[0]<v<broad[1] for v in layered)
    assert abs(layered[-1]-(straight[0]-straight[1]))<1e-12
    return {'energy_GeV':[lo,hi], 'constant_profile_difference':const,
            'all_measurable_profiles_difference':broad,'center_probabilities':straight,
            'separately_normalized_positive_energy_weights_difference':separate_weights,
            'sufficient_equal_maximum_uncontrolled_tail_fraction':margin/(1+margin),
            'entire_band_neutrino_probability':[lows[0],highs[0]],
            'entire_band_antineutrino_probability':[lows[1],highs[1]],
            'three_layer_checks':layered,'Duhamel_density_radius':amplitude_rho}


def solve():
    ch = read('../1657/independent_results.json')['extended_box_GJ_test']
    e,_=runpy.run_path(str(OLD))['charged'](ch)
    pts={'third':((5.63,7.54e-5),(5.74,7.53e-5)),
         'second':((5.8804,7.54e-5),(5.8804,7.32e-5))}
    result={}
    old_data=read('../1661/independent_results.json')
    for key,(old,new) in pts.items():
        mix_old,mass_old,props_old=physical_port(e,key,*old)
        prior=old_data['third_leg_3_to_2' if key=='third' else 'second_leg_counterpart']
        assert abs(props_old['J']-prior['J'])<1e-10
        m,mass,props=physical_port(e,key,*new)
        def in_boxes(q):
            return [2.459e-3<=q['dm31']<=2.563e-3,
                    7.232e-5<=q['dm21']<=7.544e-5,
                    .2908<=q['s2theta12']<=.3164,
                    2.37e-3<=q['dm32']<=2.51e-3]
        assert all(in_boxes(props)) and not all(in_boxes(props_old))
        second=screen(m,mass,.195,.205)
        first=screen(m,mass,.59,.61)
        broad=second['all_measurable_profiles_difference']
        assert (broad[0]>0 if key=='third' else broad[1]<0)
        folded=second['separately_normalized_positive_energy_weights_difference']
        assert (folded[0]>0 if key=='third' else folded[1]<0)
        result[key]={'original':{'values':props_old,'box_checks':in_boxes(props_old)},
                     'frequency_adjusted':props,'second_peak':second,'first_peak':first}
    return result


def compact(obj):
    if isinstance(obj,dict):
        return {k:compact(v) for k,v in obj.items()}
    if isinstance(obj,list):
        return [compact(i) for i in obj]
    if isinstance(obj,(np.bool_,bool)):
        return bool(obj)
    return float(f'{obj:.11g}')


if __name__=='__main__':
    cmd=argparse.ArgumentParser()
    cmd.add_argument('--check',action='store_true')
    args=cmd.parse_args()
    result=compact(solve())
    if args.check:
        assert result==read('independent_results.json')
        print('1668 independent PASS: both phases, calibrated frequencies, arbitrary matter profiles, all energies')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
