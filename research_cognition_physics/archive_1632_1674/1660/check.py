"""Test equivalence of two-singlet Dirac and Majorana mass mechanisms and a data bound on source overlap."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def overlap_bound(cfg):
    spec=importlib.util.spec_from_file_location('interval1486',HERE.parent/'1486/check.py')
    iv=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(iv)
    iv.SCALE=10**22
    I=iv.I
    osc=read(cfg['frequency_inputs'])
    R=(1+iv.width(osc['T2K_dm32_NO_eV2'])/iv.width(osc['KamLAND_dm21_eV2'])).sqrt()
    A=I.point(cfg['charged_electron_modulus_floor'])
    B=(1-A.square()).sqrt()
    limit=I.point(cfg['rho_in_exclusion'])
    dot=iv.hull([-limit,limit])
    g=iv.hull([(6*(1-limit.square())).sqrt(),I.point(6).sqrt()])
    small_max=0
    large_min=iv.SCALE
    large_max=0
    for k in R.parts(cfg['R_cells']):
        drive=g*(k+1/k)-8*dot
        disc=drive.square()-88
        assert disc.lo>0
        xsmall=4/(drive+disc.sqrt())
        xlarge=(drive+disc.sqrt())/22
        small=(xsmall/g-1/(3*k))/(k-1/k)
        large=(xlarge/g-1/(3*k))/(k-1/k)
        # The physical vector has a^2+b^2=1; bounded separately to allow
        # any charged rephasing, including those not realized by the texture.
        small_hi=I(small.hi,small.hi).sqrt()
        small_max=max(small_max,(small_hi+B).square().hi)
        diff=A*I(large.lo,large.lo).sqrt()-B
        assert diff.lo>0
        large_min=min(large_min,diff.square().lo)
        large_max=max(large_max,(I(large.hi,large.hi).sqrt()+B).square().hi)
    assert large_max<iv.SCALE//2
    small_A=4*I(small_max,small_max)*(1-I(small_max,small_max))
    large_A=4*I(large_min,large_min)*(1-I(large_min,large_min))
    observed=iv.width(osc['DayaBay_sin2_2theta13'])
    assert small_A.hi<observed.lo and large_A.lo>observed.hi
    return {'R_window':R.values(),'excluded_overlap_modulus_through':float(cfg['rho_in_exclusion']),
            'small_branch_reactor_amplitude_upper':small_A.hi/iv.SCALE,
            'large_branch_reactor_amplitude_lower':large_A.lo/iv.SCALE,
            'large_branch_electron_probability_upper':large_max/iv.SCALE,
            'measured_reactor_window':observed.values()}


def compute():
    cfg=read('inputs.json')
    prev=read(cfg['charged_and_majorana'])
    impl=runpy.run_path(str(HERE.parent/'1659/check.py'))
    E,_=impl['charged_basis'](read(cfg['charged_recipe']))
    assert abs(E[0,0])-float(cfg['charged_electron_modulus_floor'])>1.e-5
    W=np.array([[0.,1.],[1.,3.],[1.,1.]])
    G=W.T@W
    evals,V=np.linalg.eigh(G)
    sqrtG=(V*np.sqrt(evals))@V.T
    pinv=np.linalg.solve(G,W.T)
    omega=4/math.sqrt(22)
    out={}
    for name,block in prev['scenarios'].items():
        out[name]={}
        for sign,old in block['small_root_by_eta_deg'].items():
            eta=math.radians(float(sign))
            M,_,r=impl['rank_two'](block['dm21'],block['dm32'],eta,'small')
            C=np.diag([1,r*np.exp(1j*eta)])
            m0=float((pinv@M@pinv.T)[0,0].real)
            P=np.diag(np.exp(1j*np.r_[0.,old['mass_P_alpha_beta_rad']]))
            majorana=P@M@P
            D=m0*P.conj()@W@C.conj()@sqrtG
            H=D@D.conj().T
            np.testing.assert_allclose(H,majorana.conj().T@majorana,atol=2e-18)
            reconstructed=pinv@P@H@P.conj()@pinv.T
            row_ratio=math.sqrt(reconstructed[1,1].real/reconstructed[0,0].real)
            rho=reconstructed[0,1]/math.sqrt(reconstructed[0,0].real*reconstructed[1,1].real)
            np.testing.assert_allclose(rho,omega*np.exp(1j*eta),atol=2e-14)
            np.testing.assert_allclose(row_ratio,r*math.sqrt(11/2))
            squares, modes=np.linalg.eigh(H)
            assert abs(squares[0])<3e-18
            U=E.conj().T@modes
            angles=impl['angles'](U)
            J=(U[0,0]*U[1,1]*U[0,1].conjugate()*U[1,0].conjugate()).imag
            mbeta=math.sqrt((abs(U[0,:])**2)@squares)
            np.testing.assert_allclose(angles,old['A13_solar_atmospheric'],atol=2e-9)
            np.testing.assert_allclose([J,mbeta],[old['J_dirac'],old['m_beta_eV']],atol=2e-10)
            # On a Dirac neutrino the absolute beta mass uses H. Exact lepton
            # number forbids the light Majorana-exchange double-beta channel.
            out[name][sign]={'row_norm_ratio':row_ratio,'right_overlap_modulus':abs(rho),
                             'matched_J':J,'matched_m_beta_eV':mbeta,
                             'Majorana_light_exchange_eV':old['m_beta_beta_eV'],
                             'Dirac_light_exchange_eV':0}
    return {'Gram':G.tolist(),'shared_right_angle_deg':math.degrees(math.acos(omega)),
            'rank':int(np.linalg.matrix_rank(W)),'witnesses':out,'reactor_exclusion':overlap_bound(cfg)}


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)):return float(f'{x:.10g}')
    return x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    data=tidy(compute())
    if args.check:
        assert read('results.json')==data
        print('PASS: common left Hamiltonian, Gram inverse, beta observables, full overlap exclusion')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(data,indent=2)+'\n')
