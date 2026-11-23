"""Construct the finite three-heavy Majorana sector of the registered mixed-source witnesses."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
VU = 246.22  # GeV, explicit reference normalization for the up-type Dirac mass
EV = 1e-9


def sig(a):
    return np.linalg.svd(a, compute_uv=False)[::-1]


def frequency_cover(a,b):
    # Relative change of the type-I coefficient only; every other entry stays fixed.
    data=json.loads((HERE.parent/'1668'/'inputs.json').read_bytes())
    solar, atm = data['juno_NO_dm21'], data['juno_NO_dm31']
    # Enlarge asymmetric errors symmetrically, since only a necessary window is needed.
    win = [(solar[0]-2*max(solar[1:]))/(atm[0]+2*max(atm[1:])),
           (solar[0]+2*max(solar[1:]))/(atm[0]-2*max(atm[1:]))]
    radius=np.linalg.norm(b,2)*(1+5e-13)
    kept=[]
    pending=[(-.02,.02)]
    while pending:
        l,h=pending.pop()
        mid=(l+h)/2
        s=sig(a+(1+mid)*b)
        dev=radius*(h-l)/2+2e-13
        sl=np.maximum(s-dev,0.)**2
        su=(s+dev)**2
        if su[1]-sl[0]<win[0]*(sl[2]-su[0]) or sl[1]-su[0]>win[1]*(su[2]-sl[0]):
            continue
        if h-l>1e-6:
            pending.extend([(mid,h),(l,mid)])
        else:
            if kept and kept[-1][1]==l:
                kept[-1][1]=h
            else:
                kept.append([l,h])
    assert any(l<=0<=h for l,h in kept)
    probes={}
    for x in (-.01,-.001,0.,.001,.01):
        s=sig(a+(1+x)*b)**2
        probes[str(x)]=(s[1]-s[0])/(s[2]-s[0])
    return {'solar_atm_ratio_two_marginal_window':win,
            'allowed_relative_I_variation_outer_in_minus_plus_2_percent':kept,
            'ratio_at_fixed_relative_changes':probes}


def realize(x, data, api):
    u, d, e, _ = data[:4]
    U, D, E, F, _, ue, _, _ = api['matrices'](x,data)
    # Recover the two Yukawa weights directly from the three physical matrices.
    pencil = np.stack([D.ravel(),U.ravel()], axis=1)
    (a,b), *_ = np.linalg.lstsq(pencil,E.ravel(),rcond=None)
    np.testing.assert_allclose(a*D+b*U,E,atol=2e-15)
    rh, rf = (1-a)/b,(-3-a)/b
    H = D-F
    np.testing.assert_allclose(rh*H+rf*F,U,atol=2e-14)
    Y = rh*H-3*rf*F
    I = Y@np.linalg.solve(F,Y.T)
    fn, inn = np.linalg.norm(F),np.linalg.norm(I)
    zeta = x[7]+1j*x[8]
    C = F/fn + zeta*I/inn
    L = ue@data[4]
    low = L.conj().T@C@L.conj()
    measured = api['neutral'](low)[1]
    k = measured['common_light_scale_eV']*EV
    mL = k*F/fn
    mD = VU*Y/math.sqrt(2)
    R = -VU**2*inn/(2*k*zeta)*F
    tilt = np.linalg.solve(R,mD.T)
    schur = mL-mD@tilt
    np.testing.assert_allclose(schur,k*C,rtol=1e-11,atol=1e-23)
    norm = np.linalg.norm(tilt,2)
    smax = math.hypot(1,norm/2)+norm/2
    # Singular values of the block shear are exp(+-asinh(norm/2)).
    relative = math.expm1(2*math.asinh(norm/2))
    light = sig(schur)
    heavy = sig(R)
    assert heavy[0] > 1e7 and heavy[0] > light[-1]*smax**4
    sterile = norm/(1-light[-1]*smax**2/heavy[0])
    # Frobenius conditioning of the two separate light-mass contributions.
    total = np.linalg.norm(schur)
    sensitivity = (np.linalg.norm(mL)+np.linalg.norm(mD@tilt))/total
    return {'rH':[rh.real,rh.imag],'rF':[rf.real,rf.imag],
            'fixed_source_frequency_test':frequency_cover(L.conj().T@(F/fn)@L.conj(),
                                                         zeta*L.conj().T@(I/inn)@L.conj()),
            'light_without_active_heavy_mixing_at_boundary_eV':(light/EV).tolist(),
            'heavy_singular_masses_GeV':heavy.tolist(),
            'heaviest_vs_registered_unification_scale_5p9254e15':heavy[-1]/5.9254e15,
            'lightest_heavy_over_10PeV':heavy[0]/1e7,
            'maximum_relative_six_mass_error':relative,
            'sterile_norm_upper_per_light_vector':sterile,
            'norm_mD_GeV':np.linalg.norm(mD,2),
            'relative_light_source_sensitivity':sensitivity,
            'low_mbeta_eV':measured['mbeta_eV'],'low_mbb_eV':measured['mbb_eV']}


def compute():
    api=runpy.run_path(str(HERE/'check.py'))
    data=api['backend']()
    entries=json.loads((HERE/'inputs.json').read_bytes())['mixed_type_I_II']
    return {'up_Dirac_vev_GeV':VU,'cases':{key:realize(x,data,api) for key,x in entries.items()}}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(compute(),indent=2,allow_nan=False)+'\n'
    path=HERE/'finite_neutral_results.json'
    if args.check:
        assert json.loads(raw)==json.loads(path.read_bytes())
    else:
        path.write_bytes(raw.encode())
    print('1799 finite neutral',hashlib.sha256(raw.encode()).hexdigest()[:12])
