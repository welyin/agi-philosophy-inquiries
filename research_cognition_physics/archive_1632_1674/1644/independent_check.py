"""Independent neutral Bd CP test using the light polar block and decay interference."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def light_block(point, pole):
    s12,s23,s13,deg=point['U_input']
    c12,c23,c13=np.sqrt(1-np.square([s12,s23,s13]))
    z=cmath.exp(1j*math.radians(deg))
    U=np.array([[c12*c13,s12*c13,s13/z],
                [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])
    rb=math.sqrt(pole['right_doublet_probability'])
    rs=-complex(*point['F_in_input_basis'])/rb
    M=point['M']
    ell=np.array([0.,.1*rs/M,4.8*rb/M])
    singlet_bottom=2*(pole['left_difference']/pole['k']-abs(ell[2])**2)
    gram=np.eye(3)-np.outer(ell,ell.conj())-np.diag([0,0,singlet_bottom])
    vals,V=np.linalg.eigh(gram)
    return U@(V@np.diag(np.sqrt(vals))@V.conj().T)


def equal_box(x):
    return x*(4-11*x+x*x)/(4*(1-x)**2)-3*x**3*math.log(x)/(2*(1-x)**3)


def mixed_box(x,y):
    def f(t):
        return .25+1.5/(1-t)-.75/(1-t)**2
    return x*y*((f(x)*math.log(x)-f(y)*math.log(y))/(x-y)-.75/((1-x)*(1-y)))


def critical_amplitude(theta, weak, target):
    if math.sin(theta)>=target:
        return 0.
    roots=[]
    for sign in (-1,1):
        coefficients=[math.sin(theta-2*weak)-target,
                      2*sign*(math.sin(theta-weak)-target*math.cos(weak)),
                      math.sin(theta)-target]
        for root in np.roots(coefficients):
            if abs(root.imag)<1e-10 and 0<root.real<1:
                roots.append(root.real)
    assert roots
    r=min(roots)
    def observable(cs):
        return ((math.sin(theta)+2*r*cs*math.sin(theta-weak)+r*r*math.sin(theta-2*weak))
                /(1+2*r*cs*math.cos(weak)+r*r))
    assert abs(max(observable(1),observable(-1))-target)<1e-10
    return r


def evaluate():
    data=read('../1643/independent_results.json')
    pole=read('../1633/results.json')
    working=(1.27/80.377)**2,(172.5/80.377)**2
    xc,xt=working
    target=.66-2*(.10+.06)
    output=[]
    for point in [*data['witnesses'],data['positive_branch']]:
        K=light_block(point,pole)
        np.testing.assert_allclose(np.abs([*K[0],K[1,2]]),point['Kud_Kus_Kub_Kcb'],atol=1e-9,rtol=1e-8)
        lu,lc,lt=K[:,0].conj()*K[:,2]
        assert abs(lu+lc+lt)<1e-14
        tt=lt**2*equal_box(xt)
        ct=2*lc*lt*mixed_box(xc,xt)
        cc=lc**2*equal_box(xc)
        def coefficients(matrix):
            decay=matrix[1,2]*matrix[0,0].conjugate()
            charm=matrix[1,1].conjugate()*matrix[0,1]
            dominant=decay*charm
            small=matrix[0,2]*matrix[1,0].conjugate()*charm.conjugate()
            prod=matrix[:,0].conj()*matrix[:,2]
            mix=prod[2]**2*equal_box(xt)+2*prod[1]*prod[2]*mixed_box(xc,xt)+prod[1]**2*equal_box(xc)
            return cmath.phase(mix/dominant**2),cmath.phase(small/dominant),abs(small/dominant)
        theta, phi, ratio=coefficients(K)
        changed=np.diag(np.exp(1j*np.array([.6,-.8,1.2])))@K@np.diag(np.exp(1j*np.array([.4,1.6,-.25])))
        np.testing.assert_allclose(coefficients(changed),(theta,phi,ratio),atol=1e-12)
        r=critical_amplitude(theta,phi,target)
        delta=max(0.,math.asin(target)-theta) if theta<math.asin(target) else 0.
        assert r+1e-12>=math.sin(delta/2)
        assert abs(cmath.phase((tt+ct+cc)/tt))<=math.asin((abs(ct)+abs(cc))/abs(tt))
        output.append({'M_GeV':point['M'],'gamma_d_deg':round(point['gamma_s_gamma_d_deg'][1],4),
                       'dominant_phase_deg':round(math.degrees(theta),5),
                       'minus_eta_S_single_KK':round(math.sin(theta),7),'weak_ratio':round(ratio,8),
                       'needed_common_shift_deg':round(math.degrees(delta),4),
                       'two_free_decay_disks_min_r':round(math.sin(delta/2),6),
                       'KK_two_tree_min_r_at_assigned_target':round(r,6),
                       'KK_min_hadronic_ratio_at_assigned_target':round(r/ratio,5),
                       'box_phase_cap_deg_if_QCD_ratios_le_3':round(math.degrees(math.asin(3*(abs(ct)+abs(cc))/abs(tt))),6)})
    assert all(row['minus_eta_S_single_KK']<0 and row['KK_min_hadronic_ratio_at_assigned_target']>20 for row in output[:3])
    assert output[-1]['minus_eta_S_single_KK']>target
    return {'signed_target_for_common_fit':round(target,4),'working_xc_xt':list(working),'witnesses':output}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=evaluate()
    path=HERE/'independent_results.json'
    if args.check:
        assert result==json.loads(path.read_bytes())
        print('PASS: polar light currents, Bd boxes, rephasings, exact two-amplitude threshold')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
