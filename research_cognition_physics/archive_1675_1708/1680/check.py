"""Spectral shape and physical spin-sensitive vertices of two flavon profiles."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def parts(z):
    return np.stack((np.real(z),np.imag(z)),axis=-1).round(10).tolist()


def shape(matrix,vertex):
    left,values,right_dag=np.linalg.svd(matrix)
    right=right_dag.conj().T
    normalized=np.diag(left.conj().T@vertex@right)/values
    weights=values**2/np.sum(values**2)
    total=np.dot(weights-1/3,normalized)
    invariant=np.vdot(matrix,vertex)/np.vdot(matrix,matrix)-np.trace(np.linalg.solve(matrix,vertex))/3
    np.testing.assert_allclose(total,invariant,atol=3e-10)
    return total,normalized[::-1],values[::-1],weights[::-1]


def solve():
    cfg=read('inputs.json')
    helper=runpy.run_path(str(HERE/cfg['reconstruction']))
    pt=read(cfg['matched_point'])['joint_witness']
    q=read('../1397/inputs.json')['high_scale']
    base,error=np.array(q['down']+q['angles_theta12_theta13_theta23_delta'][:3]).T
    raw=base+error*np.array(pt['quark_pulls'])
    se,hard,_,_=helper['get_matrices'](raw,pt['quark_gamma_high'])
    sd=se.copy()
    sd[1,1]/=cfg['clebsch_C']
    phi=np.pi*cfg['vacuum_phi_over_pi']
    t=pt['common_abs_t']*np.exp(1j*phi)
    lep=(se+t*hard).T
    down=sd+hard
    we=abs(t)**2*np.linalg.norm(hard)**2/np.linalg.norm(lep)**2
    wd=np.linalg.norm(hard)**2/np.linalg.norm(down)**2
    ratio=we/wd
    assert we>1/3
    np.testing.assert_allclose(abs(t)**2*np.linalg.norm(down)**2/np.linalg.norm(lep)**2,ratio)
    slope_dict={'one_power':1j,'two_neutral_paths':4/np.sqrt(3)}
    output={}
    rng=np.random.default_rng(1680)
    def random_unitary():
        return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    rot_l,rot_r=random_unitary(),random_unitary()
    # Shared variation compatible with the two real-path implementation of Z,
    # plus an arbitrary complex E response. Test independent complex normalizations.
    omega=np.exp(4j*phi)
    drift=np.zeros((3,3),dtype=complex)
    drift[1,2]=4j*omega*hard[1,2].imag/omega.imag
    drift[2,1]=-.37+.16j
    x,y=.31-.29j,-.24+.53j
    common_down=shape(down,drift+y*down)[0]
    for name,eta in slope_dict.items():
        opt=shape(lep,(t*eta*hard).T)
        g=(t*(eta*hard+drift)).T+x*lep
        obs=shape(lep,g)
        control=shape(rot_l.conj().T@lep@rot_r,rot_l.conj().T@g@rot_r)
        np.testing.assert_allclose(control[1],obs[1],atol=4e-11)
        corrected=obs[0]-ratio*common_down
        np.testing.assert_allclose(corrected,(we-1/3)*eta,atol=2e-11)
        def finite_shape(h):
            cur=phi+h
            f=(np.exp(1j*cur) if name=='one_power' else (2+np.exp(4j*cur))/np.sqrt(3))
            mt=(se+pt['common_abs_t']*f*(hard+h*drift)).T*np.exp(x*h)
            md=(sd+hard+h*drift)*np.exp(y*h)
            def logshape(m):
                s=np.linalg.svd(m,compute_uv=False)
                return .5*np.log(s@s)-np.mean(np.log(s))
            return logshape(mt)-ratio*logshape(md)
        step=3e-5
        numeric=(finite_shape(step)-finite_shape(-step))/(2*step)
        assert abs(numeric-corrected.real)<4e-7
        output[name]={'relative_scalar_pseudoscalar_diagonals_if_other_edges_fixed':parts(opt[1]),
                      'with_Z_E_and_complex_common_drifts_shape':parts(obs[0]),
                      'combined_shape':parts(corrected),
                      'mass_singular_value_difference_derivative':round(numeric,7)}
        # Explicit additional low-energy benchmark: SM-aligned CP-even part
        # set to 1; remaining vertex matched without rescaling; no edge drift.
        tau=opt[1][-1]
        output[name]['tau_CP_deg_if_common_SM_plus_b_times_frozen_angular_vertex']={
            str(b):round(float(np.degrees(np.angle(1+b*tau))),7) for b in (0.,.3,.6)}
    # An independent complex susceptibility for the C edge can fake this particular test.
    c=sd[1,1]
    response=(9-abs(t)**2)*np.conj(c)/np.linalg.norm(lep)**2
    required=(we-1/3)*(slope_dict['two_neutral_paths']-slope_dict['one_power'])/response
    direction=np.zeros((3,3),dtype=complex)
    direction[1,1]=required
    corr=shape(lep,(t*1j*hard+3*direction).T)[0]-ratio*shape(down,direction)[0]
    np.testing.assert_allclose(corr,(we-1/3)*slope_dict['two_neutral_paths'],atol=2e-9)
    # One measured charged lepton cannot separate two complex common normalization slopes.
    tau_a=shape(lep,(t*1j*hard).T)[1][-1]
    tau_b=shape(lep,(t*slope_dict['two_neutral_paths']*hard).T)[1][-1]
    tweak=tau_a-tau_b
    np.testing.assert_allclose(shape(lep,(t*slope_dict['two_neutral_paths']*hard).T+tweak*lep)[1][-1],tau_a)
    return {'w_lepton':round(float(we),11),'w_down':round(float(wd),11),
            'shape_weight_ratio':round(float(ratio),11),'profiles':output,
            'common_down_shape_with_Z_and_E_drifts':parts(common_down),
            'C_edge_relative_log_drift_to_fake_combined_shape':parts(required/c),
            'independent_common_lepton_drift_to_fake_tau_vertex':parts(tweak)}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('results.json')
        print('1680 PASS: six physical mass vertices, basis rotations, two common drifts, shared edges, finite differences and counterfeit')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
