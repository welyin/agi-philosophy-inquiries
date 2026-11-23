"""Propagate the common charged zero into a light LNV band, with heavy pole controls."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def old(n):
    return json.loads((HERE/f"../{n}/results.json").read_bytes())


def calculate():
    osc=old(1782)
    common=old(1784)
    spot=common["two_source_witness_02"]
    # Include a guard for the 10-significant-digit archived direction.
    d=common["class_bounds"]["first_axis_tilt_sine"]+1.e-10
    box=osc["experimental_rectangles"]
    m2lo,m2hi=map(math.sqrt,box["solar_frequency_eV2"])
    m3lo,m3hi=map(math.sqrt,osc["assignments"]["NO"]["absolute_atmospheric_frequency_eV2"])
    xmin,xmax=box["isolated_weight_roots"][0]
    smin,smax=box["solar_share"]
    amin=math.sqrt(m2lo*(1-xmax)*smin)
    bmin,bmax=math.sqrt(m3lo*xmin),math.sqrt(m3hi*xmax)
    distance=d*math.sqrt(m2hi+m3hi)
    lower=amin**2-bmax**2
    upper=distance*(2*bmax+distance)
    assert amin>bmax and lower>0
    max_defect=4*math.asin(math.sqrt((distance**2-(amin-bmax)**2)/(4*amin*bmin)))
    witness_m2=math.sqrt(spot["delta21_eV2"])
    witness_m3=math.sqrt(spot["delta31_eV2"])
    x=(1-math.sqrt(1-spot["reactor_fast_amplitude"]))/2
    s=spot["calibration_s12_s23_mass_ratio"][0]
    p=witness_m2*s*(1-x)
    q=witness_m3*x
    mbb=spot["light_exchange_mbb_meV"]/1000
    relative_cos=(mbb**2-p**2-q**2)/(2*p*q)
    phase_deviation=math.acos(-relative_cos)
    assert lower < mbb < upper < p+q and phase_deviation < max_defect
    # Demonstrate the geometrical bound on explicit complex isotropic lines.
    rng=np.random.default_rng(1786)
    for _ in range(40):
        m2=rng.uniform(m2lo,m2hi)
        m3=rng.uniform(m3lo,m3hi)
        b=np.array([0,math.sqrt(m3),1j*math.sqrt(m2)])/math.sqrt(m2+m3)
        w=rng.normal(size=3)+1j*rng.normal(size=3)
        w-=b*np.vdot(b,w)
        w/=np.linalg.norm(w)
        e=math.sqrt(1-d*d)*b+d*w
        val=abs(m2*e[1]**2+m3*e[2]**2)
        tag=d*math.sqrt(m2+m3)
        assert val<=tag*(2*math.sqrt(m3)*abs(e[2])+tag)+1e-14
    # L-conserving spectator masses/couplings; only reuse the spectrum norm of 1698.
    spectator=old(1698)["shared_witness"]
    vev=246.22
    qmax=.2  # GeV, a fixed space-like propagator test, not an isotope integration
    gap=5000.  # conservative GeV lower sterile Dirac mass
    y2=sum(spectator["neutral_pole_Yukawa_squared"])
    spectator_D2=vev**2*y2/2
    leg=spectator_D2/gap**2
    light_meV=m3hi*1000
    # PeV singlet bound: complete active electron mass row, including spectators.
    high_GeV=1.e7
    pole_meV=1e12*qmax**2*(4*math.pi*vev**2/2+spectator_D2)/high_GeV**3
    # The 1784 point has a much smaller two-source Dirac norm.
    coeff=np.array(spot["raw_real_D_coefficients_before_common_scale"])
    dimless=coeff*.5**np.array([[2,0,1],[4,2,3]])
    high_m0=2.e7
    point_D2=spot["common_mass_scale_eV"]*1e-9*high_m0*np.vdot(dimless,dimless).real
    point_GeV=1.23e7  # below 0.618*20 PeV, allowing bounded tiny mixing
    point_meV=1e12*qmax**2*(point_D2+spectator_D2)/point_GeV**3
    half_life=(3.8e26*(28/(1000*upper))**2,3.8e26*(122/(1000*lower))**2)
    return dict(light_band_meV=[1000*lower,1000*upper],
                relative_Majorana_contribution_distance_from_pi_deg_max=math.degrees(max_defect),
                witness_meV=mbb*1000,witness_Majorana_phase_deviation_deg=math.degrees(phase_deviation),
                aligned_phase_same_oscillation_meV=(p+q)*1000,
                generic_oscillation_upper_meV=(m2hi*smax*(1-xmax)+m3hi*xmax)*1000,
                spacelike_qmax_GeV=qmax,heavy_electron_row_yukawa_squared_budget=4*math.pi,
                all_PeV_exchange_bound_meV=pole_meV,
                witness_fixed_high_D_Frobenius_squared_GeV2=point_D2,
                witness_PeV_exchange_bound_with_spectators_meV=point_meV,
                spectator_squared_Yukawa_sum=y2,largest_static_leg_deviation=leg,
                unmatched_static_spectator_shift_meV=2*leg*light_meV,
                residual_after_static_matching_meV=2*leg*(qmax/gap)**2*light_meV,
                xenon_light_only_conditional_years=half_life)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    payload=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if check:
        assert json.loads(dest.read_bytes())==json.loads(payload)
    else:
        dest.write_text(payload,encoding="utf-8")
    print("1786 charged texture Majorana interval",hashlib.sha256(payload.encode()).hexdigest()[:12])
