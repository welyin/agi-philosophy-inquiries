"""Two resolved sterile poles: exact seesaw minors and charged-kaon line strengths."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"../1386/inputs.json").read_bytes())
HEAVY=np.array([250.,300.])*1e6  # eV, chosen and resolved, not inferred
ROW_NUISANCE=.01
SMALL_ROW=1e-5


def rotation(s12,s13,s23,delta):
    a,b,c=np.sqrt([s12,s13,s23])
    A,B,C=np.sqrt(1-np.array([s12,s13,s23]))
    z=np.exp(1j*delta)
    return np.array([[A*B,a*B,b/z],[-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                     [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])


def small_branch(amplitude):
    return (1-math.sqrt(1-amplitude))/2


def mass_pair(order,solar,atm):
    return np.sqrt([solar,solar+atm] if order=="NO" else [atm,atm+solar])


def kaon_factor(neutral,charged):
    # Ratio of K -> l N phase-space/helicity factor to K -> l nu_massless.
    x=(charged/493.677)**2
    y=(neutral/493.677)**2
    return (x+y-(x-y)**2)*math.sqrt((1-x-y)**2-4*x*y)/(x*(1-x)**2)


def construction(order):
    d=CFG[order]
    weights=d["atmospheric_sine_squared"]
    m=mass_pair(order,d["solar_gap"][0],d["atmospheric_gap"][0])
    mixing=rotation(d["solar_sine_squared"][0],small_branch(CFG["DayaBay_amplitude"][0]),weights[0],math.pi/2)
    active_cols=[1,2] if order=="NO" else [0,1]
    zero=3-sum(active_cols)
    u=np.zeros((5,5),dtype=complex)
    u[:3,zero]=mixing[:,zero]
    masses=np.zeros(5)
    for j,i in enumerate(active_cols):
        p=math.sqrt(HEAVY[j]/(HEAVY[j]+m[j]))
        q=math.sqrt(m[j]/(HEAVY[j]+m[j]))
        u[:3,i]=1j*p*mixing[:,i]
        u[3+j,i]=-1j*q
        u[:3,3+j]=q*mixing[:,i]
        u[3+j,3+j]=p
        masses[i],masses[3+j]=m[j],HEAVY[j]
    mass=(u.conj()*masses)@u.conj().T
    assert np.linalg.norm(u.conj().T@u-np.eye(5))<1e-14
    assert np.max(abs(mass[:3,:3]))<1e-15
    assert np.max(abs((mass@u-u.conj()*masses)[:,:3]))<1e-8
    c=u[:2,active_cols]
    b=u[:2,3:]
    minor=abs(np.linalg.det(c))**2
    lhs=HEAVY.prod()*abs(np.linalg.det(b))**2
    rhs=m.prod()*minor
    assert math.isclose(lhs,rhs,rel_tol=1e-12)
    x=abs(b)**2
    norm=1-x.sum(axis=1)
    charged=[.51099895,105.6583755]
    rates=[[kaon_factor(HEAVY[j]/1e6,charged[i])*x[i,j]/norm[i] for j in range(2)] for i in range(2)]
    shape=abs(u[:2,:3])**2/norm[:,None]
    return dict(massive_light_eV=m.tolist(), conditional_light_shape_distortion=float(np.max(abs(shape-abs(mixing[:2,:])**2))),
                zero_active_block_residual_eV=float(np.max(abs(mass[:3,:3]))),
                electron_muon_heavy_line_mixings=x.tolist(),K_heavy_over_light_rates=rates,
                measured_minor_light_squared=minor,mass_weighted_minor_eV2=lhs)


def accessible_bound(order):
    d=CFG[order]
    a=CFG["DayaBay_amplitude"]
    s13hi=small_branch(a[0]+2*a[1])
    s12lo=d["solar_sine_squared"][0]-2*d["solar_sine_squared"][1]
    # Work box uses the larger asymmetric one-sided atmospheric half-width on both sides.
    a23=d["atmospheric_sine_squared"]
    s23lo,s23hi=a23[0]-2*max(a23[1:]),a23[0]+2*max(a23[1:])
    if order=="NO":
        amp=max(0,math.sqrt(s12lo*s23lo)-math.sqrt((1-s12lo)*(1-s23lo)*s13hi))
        tau_zero_min=amp**2
    else:
        tau_zero_min=(1-s23hi)*(1-s13hi)
    gap=d["solar_gap"]
    atm=d["atmospheric_gap"]
    masses=mass_pair(order,gap[0]-2*gap[1],atm[0]-2*atm[1])
    assert tau_zero_min>2*ROW_NUISANCE
    minor=(1-SMALL_ROW)**2*(tau_zero_min-2*ROW_NUISANCE)-2*SMALL_ROW-SMALL_ROW**2
    root=math.sqrt(masses.prod()*minor/HEAVY.prod())
    max_line=root/2
    assert 0<max_line<SMALL_ROW/2
    return dict(two_width_mass_lower_eV=masses.tolist(),unitary_reference_tau_zero_floor=tau_zero_min,
                physical_light_minor_squared_lower=minor,heavy_minor_magnitude_lower_small_rows=root,
                guaranteed_at_least_one_line_above=max_line,uniform_upper_limit_product_needed=max_line**2)


def run():
    return dict(heavy_poles_MeV=(HEAVY/1e6).tolist(),shape_transfer_per_row=ROW_NUISANCE,
                small_heavy_e_or_mu_row_max=SMALL_ROW,
                mass_orderings={key:{"tree_example":construction(key),"area_bound":accessible_bound(key)}
                                for key in ("NO","IO")})


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    payload=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if check:
        assert json.loads(dest.read_bytes())==json.loads(payload)
    else:
        dest.write_text(payload,encoding="utf-8")
    print("1782 exact two-pole charged-kaon minor",hashlib.sha256(payload.encode()).hexdigest()[:12])
