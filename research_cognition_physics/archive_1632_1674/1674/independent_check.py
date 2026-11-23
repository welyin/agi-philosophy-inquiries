"""Cross-check neutral singlet poles, charged-current residues, and identifiability."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def full_spectrum(electron, ratio, g=1.):
    omega=np.exp(-2j*np.pi/3)
    z=(3+2*omega)/5
    q=abs(z)
    t=g*ratio
    M=g*np.sqrt(q)
    ports=(1.,1.,z)
    mass=np.zeros((9,6),dtype=complex)
    for j,p in enumerate(ports):
        mass[3*j:3*j+3,2*j:2*j+2]=[[t,0],[M,g*p],[-g*np.conjugate(p),M]]

    # Full SVD and physical charged-current amplitudes on all nine states.
    modes,masses,_=np.linalg.svd(mass,full_matrices=True)
    mass2=np.r_[masses**2,np.zeros(3)]
    current=electron.conj().T@modes[::3]
    projected=lambda mask: current[:,mask]@current[:,mask].conj().T
    residue_light=projected(mass2==0)
    v=electron[2,:].conjugate()
    axis=np.outer(v,v.conjugate())
    all_visible=[]
    qs=[1.,q]
    for i,port_abs in enumerate(qs):
        H2=M*M+g*g*port_abs**2
        vis2=H2+t*t
        # The test ratios are chosen away from accidental heavy pole coincidences.
        invisible=projected(np.isclose(mass2,H2,rtol=2e-12,atol=1e-15*g*g))
        measured=projected(np.isclose(mass2,vis2,rtol=2e-12,atol=1e-15*g*g))
        np.testing.assert_allclose(invisible,0,atol=1e-11)
        expected=(np.eye(3)-axis if i==0 else axis)*t*t/vis2
        np.testing.assert_allclose(measured,expected,atol=1e-11)
        all_visible.append((vis2,measured))
    np.testing.assert_allclose(sum((r for h,r in all_visible),residue_light),np.eye(3),atol=1e-11)
    np.testing.assert_allclose(sum((h*r for h,r in all_visible)),t*t*np.eye(3),atol=1e-10*g*g)
    u=t*t/(g*g*(q+1))
    f=(1-q)*u/(q+u)
    predicted=(np.eye(3)-f*axis)/(1+u)
    np.testing.assert_allclose(residue_light,predicted,atol=1e-11)
    rows=np.diag(residue_light).real
    weights=abs(v)**2
    r=rows[0]/rows[1]
    if abs(weights[0]-weights[1])>1e-5:
        inferred=(r-1)/(r*weights[1]-weights[0])
        u_from_ratio=q*inferred/(1-q-inferred)
        np.testing.assert_allclose(u,u_from_ratio,rtol=4e-10)
    appearance=abs(residue_light)**2/np.outer(rows,rows)
    return {'residues_e_mu_tau':rows,'electron_to_muon_width_ratio_factor':r,
            'tau_to_muon_width_ratio_factor':rows[2]/rows[1],
            'normalized_zero_distance_appearance_e_mu':appearance[0,1],
            'normalized_zero_distance_appearance_mu_tau':appearance[1,2],
            'inclusive_neutral_current_trace_in_units_of_true_tree_coupling':
                np.trace(residue_light@residue_light).real,
            'physical_visible_poles':np.sqrt([h for h,r in all_visible]),
            'u':u}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    return np.asarray(x).round(13).tolist()


def solve():
    helper=runpy.run_path(str(HERE.parent/'1659/independent_check.py'))
    e,_=helper['charged'](read('../1657/independent_results.json')['extended_box_GJ_test'])
    record={}
    for ratio in (.03,.35,1.7):
        point=full_spectrum(e,ratio)
        dilated=full_spectrum(e,ratio,10.)
        for name,val in point.items():
            factor=10. if name=='physical_visible_poles' else 1.
            np.testing.assert_allclose(np.asarray(val)*factor,dilated[name],atol=1.e-11)
        record[str(ratio)]=rounded(point)
    fourier=np.exp(2j*np.pi*np.outer(range(3),range(3))/3)/np.sqrt(3)
    blind=full_spectrum(fourier,.35)
    rows=blind['residues_e_mu_tau']
    np.testing.assert_allclose(rows,rows[0],atol=1.e-12)
    assert rows[0]<1 and blind['normalized_zero_distance_appearance_mu_tau']>1.e-6
    return {'frozen_charged_basis_third_row_weights':rounded(abs(e[2,:])**2),
            'frozen_examples_t_over_g':record,'balanced_basis_control_t_over_g_035':rounded(blind),
            'scale_dilation_factor_checked':10}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1674 independent PASS: full active residues, spectral moments, Fourier control and scale dilation')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
