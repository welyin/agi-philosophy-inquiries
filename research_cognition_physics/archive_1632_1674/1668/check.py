"""Use the two mass frequencies and the two signs of appearance to test frozen flavour links."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
OLD = runpy.run_path(str(ARCH/'1659/check.py'))
PHASE = runpy.run_path(str(ARCH/'1661/check.py'))
INPUT = json.loads((HERE/'inputs.json').read_bytes())


def read(path):
    return json.loads((HERE/path).read_bytes())


def window(row, width):
    return np.array([row[0]-width*row[1], row[0]+width*row[2]])


def gap_box(width):
    solar=window(INPUT['juno_NO_dm21'],width)
    heavy=window(INPUT['juno_NO_dm31'],width)
    return np.sqrt([heavy[0]/solar[1], heavy[1]/solar[0]])


def allowed_r(charged, phases, sign, angles):
    grid=np.linspace(5.4924,6.2245,65537)
    centers=.5*(grid[1:]+grid[:-1])
    h=np.diff(grid).max()/2
    coefs=PHASE['coefficients'](charged,phases)
    r=PHASE['root'](centers)
    phase=np.exp(1j*sign*2*math.pi/3)
    ok,_=PHASE['enclose'](coefs,phase,centers,r,h,angles)
    opposite,_=PHASE['enclose'](coefs,phase.conjugate(),centers,r,h,angles)
    assert ok.any() and not opposite.any()
    return np.array([centers[ok].min()-h, centers[ok].max()+h])


def neutrino_mix(charged, phases, sign, R):
    # Phase orientation is that of the symmetric mass, hence the conjugate
    # transfer multiplies the Takagi modes in the observable mixing matrix.
    core,modes,_=OLD['rank_two'](1.e-4,(R**2-1)*1.e-4,sign*2*np.pi/3,'small')
    gram=core.conj().T@core
    assert np.linalg.norm(gram@modes-modes*(1.e-4*np.array([0,1,R**2])))<1e-14
    return charged.conj().T @ (np.exp(-1j*np.r_[0.,phases])[:,None]*modes)


def phase_delta(v):
    p=abs(v)**2
    z=p[0,2]
    s=p[0,1]/(1-z)
    t=p[1,2]/(1-z)
    j=(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    d=math.sqrt(s*(1-s)*t*(1-t)*z)
    sin_delta=j/((1-z)*d)
    cos_delta=(p[1,0]-s*(1-t)-(1-s)*t*z)/(2*d)
    assert math.isclose(sin_delta**2+cos_delta**2,1,abs_tol=1e-11)
    return math.atan2(sin_delta,cos_delta)/math.pi


def lin_cells(a,b,max_step):
    count=max(1,math.ceil((b-a)/max_step))
    points=np.linspace(a,b,count+1)
    return (points[1:]+points[:-1])/2,(b-a)/(2*count)


def appearance_amplitudes(H):
    eig, modes=np.linalg.eigh(H)
    adj=modes.swapaxes(-1,-2).conj()
    assert np.max(np.linalg.norm(H@modes-modes*eig[...,None,:],axis=(-1,-2)))<1e-11
    assert np.max(np.linalg.norm(modes@adj-np.eye(3),axis=(-1,-2)))<1e-12
    return abs(np.sum(modes[...,0,:]*modes[...,1,:].conjugate()*np.exp(-1j*eig),axis=-1))


def certify_band(charged, phases, sign, Rbox):
    prop=INPUT['propagation']
    e0,e1=prop['energy_GeV']
    conv=prop['km_per_inverse_eV_conversion']*prop['baseline_km']
    g=conv/(2e9)
    v_unit=conv*prop['V_per_rhoYe_eV']
    solar=window(INPUT['juno_NO_dm21'],2)
    heavy=window(INPUT['juno_NO_dm31'],2)
    rcenter,rh=lin_cells(*Rbox,.000075)
    rho=[sum(prop['rhoYe_g_cm3'])/2]
    ph=(prop['rhoYe_g_cm3'][1]-prop['rhoYe_g_cm3'][0])/2
    projector_e=np.diag([1.,0.,0.])
    lower,upper=1.,-1.
    independent_low,independent_high=1.,-1.
    count=0
    for R in rcenter:
        # Enclose every allowed d21 in this R cell before making a grid in d21/E.
        ds_min=max(solar[0],heavy[0]/(R+rh)**2)
        ds_max=min(solar[1],heavy[1]/(R-rh)**2)
        if ds_min>ds_max:
            continue
        u,uh=lin_cells(g*ds_min/e1,g*ds_max/e0,.00004)
        mix=neutrino_mix(charged,phases,sign,R)
        mass=(mix*np.array([0.,1.,R**2]))@mix.conj().T
        vac=u[:,None,None,None]*mass
        potential=v_unit*np.array(rho)[None,:,None,None]*projector_e
        nu=appearance_amplitudes(vac+potential)
        anti=appearance_amplitudes(vac.conjugate()-potential)
        r_cap=R+rh
        # Uniform norm bound on dP3/dR < 2; see proof. Common identity
        # phases do not affect the modulus of a transition amplitude.
        err=((u.max()+uh)*(2*r_cap+2*(r_cap*r_cap-1))*rh
             +.5*r_cap**2*uh+.5*v_unit*ph+1e-6)
        nlo,nhi=np.maximum(0,nu-err)**2,np.minimum(1,nu+err)**2
        alo,ahi=np.maximum(0,anti-err)**2,np.minimum(1,anti+err)**2
        low,high=nlo-ahi,nhi-alo
        lower=min(lower,low.min())
        upper=max(upper,high.max())
        # Beam energy weights may differ; both beams still propagate the
        # same neutrino masses. Keep that shared solar frequency per slice.
        slices=np.linspace(ds_min,ds_max,9)
        for slo,shi in zip(slices[:-1],slices[1:]):
            mask=(u+uh+1e-13 >= g*slo/e1) & (u-uh-1e-13 <= g*shi/e0)
            independent_low=min(independent_low,nlo[mask].min()-ahi[mask].max())
            independent_high=max(independent_high,nhi[mask].max()-alo[mask].min())
        count+=nu.size
    assert count>0
    return {'enclosure':[lower,upper],'independent_weights_enclosure':[independent_low,independent_high],
            'grid_cells':count}


def compute():
    charged,_=OLD['charged_basis'](read(INPUT['charged']))
    saved=read(INPUT['phase_witnesses'])
    angles=np.array(read(INPUT['angle_windows'])['joint_windows'])
    branches={}
    for name,sign,R,ds in INPUT['branches']:
        phi=-math.pi/6
        if name=='third_leg_3_to_2':
            phases=np.array([2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))])
        else:
            phases=np.array([10*phi+math.radians(11.5),7*phi])
        np.testing.assert_allclose(phases,saved[name]['mass_P_rad'],rtol=0,atol=1e-9)
        outer=allowed_r(charged,phases,sign,angles)
        juno=gap_box(2)
        clip=[max(outer[0],juno[0]),min(outer[1],juno[1])]
        assert clip[0]<R<clip[1]
        mix=neutrino_mix(charged,phases,sign,R)
        d31=R**2*ds
        d32=d31-ds
        observed=OLD['angles'](mix)
        assert np.all((observed>=angles[:,0]) & (observed<=angles[:,1]))
        assert window(INPUT['juno_NO_dm21'],2)[0]<=ds<=window(INPUT['juno_NO_dm21'],2)[1]
        assert window(INPUT['juno_NO_dm31'],2)[0]<=d31<=window(INPUT['juno_NO_dm31'],2)[1]
        assert window(INPUT['accelerator_NO_dm32'],1)[0]<=d32<=window(INPUT['accelerator_NO_dm32'],1)[1]
        delta=phase_delta(mix)
        dlo,dhi=INPUT['accelerator_NO_delta_over_pi']
        assert any(dlo<=delta+2*i<=dhi for i in (-1,0,1))
        band=certify_band(charged,phases,sign,clip)
        branches[name]={'angle_allowed_R_outer':outer.tolist(),'with_JUNO_two_widths':clip,
                        'witness':{'R':R,'d21':ds,'d31':d31,'d32':d32,
                                   'delta_over_pi':delta,'three_angles':observed.tolist()},
                        'nu_minus_antinu':band}
    first,second=branches.values()
    assert first['angle_allowed_R_outer'][1]<gap_box(1)[0]
    assert first['nu_minus_antinu']['enclosure'][0]>0>second['nu_minus_antinu']['enclosure'][1]
    assert (first['nu_minus_antinu']['independent_weights_enclosure'][0] > 0
            > second['nu_minus_antinu']['independent_weights_enclosure'][1])
    return {'JUNO_R_one_width':gap_box(1).tolist(),'JUNO_R_two_widths':gap_box(2).tolist(),
            'branches':branches}


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f'{x:.11g}')
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: two reactor-gap comparisons, full CP quadrants, finite energy and density certificate')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
