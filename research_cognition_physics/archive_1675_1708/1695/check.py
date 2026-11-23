"""Use all four poles of the 1689 chain in a specifically declared non-SUSY Higgs matching."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def port_data(c, h, d):
    mass=np.block([[np.zeros_like(c),c],[c.T,h]])
    port=np.column_stack([d,np.zeros_like(d)])
    return mass,port


def wilson(mass, port):
    f=np.linalg.solve(mass.T,port.T).T
    return -f@port.T, f@f.conj().T


def calc():
    cfg=read("inputs.json")
    chain=read(cfg["double_seesaw"])
    old=read(cfg["double_seesaw_reference"])
    reg=read(cfg["registered_matching"])
    aset=read(cfg["higgs_pole"])["atlas"]["mass_GeV"]
    v=reg["v246_GeV"]
    u=v/math.sqrt(2)
    mH,err=aset
    mlo,mhi=mH-2*err,mH+2*err
    wl,wh=reg["positive_Higgs_threshold_transport"]
    unit=chain["sterile_unit_GeV"]
    flavor=read(chain["flavor_fields"])
    s=flavor["scalar_magnitude_over_cutoff"]*np.exp(1j*np.pi*flavor["scalar_angle_over_pi"])
    x=flavor["real_X_over_cutoff"]
    k=chain["coupling_to_s2"]/chain["m_N0H0"]
    c=unit*np.array([[chain["m_N0H0"],0],[chain["coupling_to_s2"]*s**2,chain["coupling_to_x"]*x]])
    h=unit*np.diag(chain["bare_H0_H6_masses_in_sterile_units"])
    nold=read(chain["old_model"])
    r=read(chain["old_result"])["rank_two_witness"]["small_root_relative_amplitude"]
    nu=read(chain["neutrino_calibration"])["joint_witness"]
    a,b=np.array(nold["source_a"]),np.array(nold["source_b"])
    rows=(s/abs(s))**np.array(nold["su5_generation_charges_5bar"])
    d=rows[:,None]*np.column_stack([math.sqrt(r)*b,s**2*(a*x/abs(s)**2+k*math.sqrt(r)*b)])
    _,port0=port_data(c,h,d)
    mass=port_data(c,h,d)[0]
    strength=np.sqrt(math.sqrt(nu["d21"])*1e-9/np.linalg.svd(wilson(mass,port0)[0],compute_uv=False)[1])
    d=d*strength
    np.testing.assert_allclose(strength/unit,old["active_D_prefactor_in_sterile_units"],rtol=1.e-9,atol=0)
    mass,port=port_data(c,h,d)
    light,w=wilson(mass,port)
    light_spec=np.linalg.svd(light,compute_uv=False)
    np.testing.assert_allclose(light_spec[:2]*1e9, np.sqrt(nu["d21"])*np.array([nu["R"],1]),rtol=5e-10)
    denom=8*np.pi**2*u**2
    def threshold_and_poles(m,y):
        vals, modes=np.linalg.eigh(m.conj().T@m)
        rows=np.sum(abs(y@modes)**2,axis=0)
        shares=vals*rows/denom
        np.testing.assert_allclose(sum(shares),np.linalg.norm(y@m.conj().T)**2/denom,rtol=2e-13)
        return sum(shares),np.sqrt(vals),shares
    total,heavy,shares=threshold_and_poles(mass,port)
    yukawa_norm2=shares*8*np.pi**2/heavy**2
    soft_mass2=1e6
    soft_curvature=sum(yukawa_norm2*((heavy**2+soft_mass2)*np.log1p(soft_mass2/heavy**2)
                                     -soft_mass2))/(8*np.pi**2)
    assert soft_curvature>0
    np.testing.assert_allclose(heavy[::-1],old["physical_sterile_singular_masses_EW_off_GeV"])
    np.testing.assert_allclose(total,np.linalg.norm(d@c.conj())**2/denom,rtol=2e-13)
    k_mass=-c@np.linalg.solve(h,c.T)
    baseline=wilson(k_mass,d)[1]
    excess=w-baseline
    q=d@d.conj().T
    flavor_moment=(d@c.conj())@(d@c.conj()).conj().T
    np.testing.assert_allclose(q@np.linalg.pinv(excess,rcond=1e-12)@q,flavor_moment,rtol=1e-12)
    assert np.trace(excess).real*denom*total>=np.trace(q).real**2
    mass2,port2=port_data(2*c,4*h,d)
    high2=threshold_and_poles(mass2,port2)[0]
    m2,w2=wilson(mass2,port2)
    np.testing.assert_allclose(m2,light,rtol=1e-12,atol=1e-25)
    np.testing.assert_allclose(w2-baseline,excess/4,rtol=2e-12,atol=1e-29)
    np.testing.assert_allclose(high2,4*total)
    # The unique minimal two-singlet completion of both low Wilson matrices.
    val,evec=np.linalg.eigh(w)
    f=evec[:,-2:]*np.sqrt(val[-2:])
    inv=np.linalg.pinv(f)
    m_fake=-inv@light@inv.T
    d_fake=f@m_fake
    lm,wm=wilson(m_fake,d_fake)
    np.testing.assert_allclose(lm,light,rtol=2e-12,atol=1e-25)
    np.testing.assert_allclose(wm,w,rtol=2e-12,atol=1e-29)
    fake_thr,fake_poles,_=threshold_and_poles(m_fake,d_fake)
    # Change whole neutral unit and Y so the light mass is identical. Never identify this with a measured scale.
    bcentral=(mH*mH/total)**(1/3)
    blo=(mlo*mlo/(total*wh))**(1/3)
    bhi=(mhi*mhi/(total*wl))**(1/3)
    mres,wr=wilson(bcentral*mass,math.sqrt(bcentral)*port)
    np.testing.assert_allclose(mres,light,rtol=1e-12,atol=1e-25)
    np.testing.assert_allclose(wr,w/bcentral,rtol=1e-12,atol=1e-30)
    fixed_cutoff=chain["flavor_cutoff_over_sterile_unit"]*unit
    bcut=fixed_cutoff/max(heavy)
    t=1/64
    bfixed=bcentral/t**(1/3)
    c_fixed=bfixed*math.sqrt(t)*c
    h_fixed=bfixed*t*h
    d_fixed=math.sqrt(bfixed)*d
    mass_fixed,port_fixed=port_data(c_fixed,h_fixed,d_fixed)
    new_light,new_w=wilson(mass_fixed,port_fixed)
    new_threshold,new_poles,new_shares=threshold_and_poles(mass_fixed,port_fixed)
    np.testing.assert_allclose(new_light,light,rtol=2e-12,atol=1e-25)
    np.testing.assert_allclose(new_w,(baseline+excess/t)/bfixed,rtol=2e-12,atol=1e-29)
    np.testing.assert_allclose(new_threshold,mH**2,rtol=3e-14)
    assert max(new_poles)<fixed_cutoff
    active_ymax=np.linalg.norm(d_fixed,2)/u
    a2_vev=abs(s)**2*fixed_cutoff
    x_vev=x*fixed_cutoff
    s_yuk=abs(c_fixed[1,0])/a2_vev
    x_yuk=abs(c_fixed[1,1])/x_vev
    b_vev=fixed_cutoff/2
    h_yuk=np.abs(np.diag(h_fixed))/b_vev
    assert max(active_ymax,s_yuk,x_yuk,*h_yuk)<1
    return {
        "per_pole_in_ascending_mass_order": [
            {"mass_GeV": round(m,6),"threshold_GeV2":float(f"{sd:.10g}"),
             "fraction":round(sd/total,8)} for m,sd in zip(heavy,shares)],
        "benchmark_sum_GeV2": float(f"{total:.12g}"),
        "soft_split_common_1TeV_scalar_CW_abs_H_squared_coefficient_GeV2":float(f"{soft_curvature:.12g}"),
        "benchmark_max_transport_over_lower_ATLAS_square":float(f"{wh*total/mlo**2:.10g}"),
        "scaled_unit_central_and_window_PeV":[round(unit*z/1e6,8) for z in [bcentral,blo,bhi]],
        "scaled_poles_central_PeV":[round(m*bcentral/1e6,8) for m in heavy],
        "scaled_largest_Yukawa_singular_value":float(f"{np.linalg.norm(d,2)*math.sqrt(bcentral)/u:.10g}"),
        "central_total_weak_weight":float(f"{np.trace(wr).real:.12g}"),
        "full_to_minimal_two_state_threshold_ratio":round(total/fake_thr,9),
        "minimal_two_state_fake_poles_GeV":[round(z,7) for z in fake_poles],
        "two_state_at_full_scale_relative_to_target":round(fake_thr/total,10),
        "allowed_positive_transport_bands_disjoint_for_both_models":bool(total*wl>fake_thr*wh),
        "fixed_previous_flavor_cutoff_PeV":fixed_cutoff/1e6,
        "fixed_cutoff_max_neutral_scale_ratio":round(bcut,7),
        "fixed_cutoff_max_fraction_higgs_lower_target":round(wh*total*bcut**3/mlo**2,9),
        "fixed_cutoff_changed_heavy_shape_witness": {
            "t":t,"b":round(bfixed,7),
            "poles_PeV":[round(z/1e6,8) for z in new_poles],
            "per_pole_threshold_GeV2":[round(z,6) for z in new_shares],
            "trace_weak_current_weight":float(f"{np.trace(new_w).real:.10g}"),
            "largest_active_Yukawa":float(f"{active_ymax:.10g}"),
            "A2_X_neutral_Yukawas":[round(s_yuk,8),round(x_yuk,8)],
            "B_vev_PeV_and_two_Yukawas":[b_vev/1e6,*[round(y,8) for y in h_yuk]],
        },
        "two_step_shape_deformation_C_times_2_H_times_4_fixed_D":
            {"Higgs_threshold_ratio":round(high2/total,10),
             "excess_current_weight_ratio":round(np.trace(w2-baseline).real/np.trace(excess).real,10)}
    }


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=calc()
    if args.check:
        assert result==read("results.json")
        print("1695 checks passed")
    else:
        print(json.dumps(result,indent=2,allow_nan=False))
