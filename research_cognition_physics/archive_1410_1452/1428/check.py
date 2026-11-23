"""Number of resolvable eigenphases of a finite coherent transfer measured by reciprocal rank-one ports."""
import argparse
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def load_transport():
    spec=spec_from_file_location("profile1423",HERE.parent/"1423"/"check.py")
    module=module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def two_phase_witness(s):
    # Orthogonal e, v; physical second port has overlap s with e.
    e=np.array([1,0],complex)
    v=np.array([0,1],complex)
    other=s*e+np.sqrt(1-s*s)*v
    line=(e+1j*v)/np.sqrt(2)
    unitary=np.eye(2)+(1j-1)*np.outer(line,line.conj())
    delta=abs(np.vdot(other,unitary@e))**2-abs(np.vdot(e,unitary@other))**2
    assert np.isclose(abs(delta),2*s*np.sqrt(1-s*s))
    return float(abs(delta))


def run():
    read=lambda p: json.loads((HERE/p).read_text(encoding="utf-8"))
    inputs=read("inputs.json")
    d=read(inputs["profiles_and_masses"])
    s=np.sqrt(read(inputs["neutral_port_bound"])["precision_comparison"]["Z_reanalysis"]["zero_distance_upper"])
    ceiling=two_phase_witness(s)
    old=load_transport()
    positions=np.asarray(d["position_km"])
    lengths=np.diff(np.r_[positions[0],(positions[:-1]+positions[1:])/2,positions[-1]])
    phase=1e-13/(2*d["hbarc_GeV_cm"])
    matter=2e18*np.sqrt(2)*d["GF"]*d["NA"]*d["hbarc_GeV_cm"]**3*d["electron_fraction"]
    s13=np.sqrt((1-np.sqrt(1-d["theta13_amplitude"]))/2)
    alpha,alo,ahi=d["alpha_rad"]
    sin_min=min(np.sin(alo),np.sin(ahi))
    report={}
    for lo,hi in inputs["finite_windows_GeV"]:
        energies=np.linspace(lo,hi,inputs["intervals"]+1)
        idx=len(energies)//2
        for profile,densities in d["density_g_cm3"].items():
            for order,(x12,x23,ds,da) in d["mass_and_angles"].items():
                r=old.rotation(0,2,s13)@old.rotation(0,1,np.sqrt(x12))
                m=(r*[0,ds,da])@r.T
                fw=old.transport(m,energies,lengths,densities,phase,matter)
                bw=old.transport(m,energies,lengths[::-1],densities[::-1],phase,matter)
                assert np.allclose(bw,fw.transpose(0,2,1))
                odd=(np.imag(fw[:,1,0].conj()*fw[:,2,0])+np.imag(bw[:,1,0].conj()*bw[:,2,0]))/2
                lip=phase*lengths.sum()*np.ptp([0,ds,da])/lo**2
                margin=lip*(hi-lo)/(2*inputs["intervals"])
                assert min(odd)*max(odd)>0 and min(abs(odd))>margin
                lower=4*np.sqrt(x23*(1-x23))*sin_min*(min(abs(odd))-margin)
                # Across the three-flavor target's two matter directions the diagonal
                # electron survival is identical. The test itself uses measured appearance
                # sums, so using this forecast does not grant row closure to competitors.
                appearance_sum_upper=2*(max(1-abs(fw[:,0,0])**2)+margin)
                rate_dependent_ceiling=2*s*np.sqrt(2*appearance_sum_upper)
                lo_error,hi_error=0.,lower/2
                for _ in range(65):
                    error=(lo_error+hi_error)/2
                    if lower-2*error>2*s*np.sqrt(2*(appearance_sum_upper+2*error)):
                        lo_error=error
                    else:
                        hi_error=error
                for delta in (-alpha,alpha,np.pi-alpha,np.pi+alpha):
                    u=old.rotation(1,2,np.sqrt(x23))@np.diag([1,1,np.exp(1j*delta)])
                    f=u@fw[idx]@u.conj().T
                    b=u@bw[idx]@u.conj().T
                    direct=(abs(f[1,0])**2-abs(f[0,1])**2+abs(b[1,0])**2-abs(b[0,1])**2)/2
                    assert np.isclose(direct,-4*np.sqrt(x23*(1-x23))*odd[idx]*np.sin(delta))
                advantage=max(0.,lower-ceiling)
                report[f"{lo}-{hi}_{profile}_{order}"]={
                    "absolute_four_rate_lower":float(lower),"strict_gap_to_two_phases":float(advantage),
                    "per_rate_common_absolute_error_work_budget":float(advantage/2),
                    "target_mean_pair_rate_sum_upper":float(appearance_sum_upper),
                    "rate_dependent_ceiling_at_target":float(rate_dependent_ceiling),
                    "rate_dependent_per_probability_error_budget":float(lo_error),
                    "some_transfer_operator_distance_floor":float(advantage/2),
                    "some_three_dim_transfer_min_pair_angle_floor_rad":float(4*np.arcsin(advantage/4))}
    return {"round":1428,"max_port_overlap_assumed":float(s),"two_phase_absolute_ceiling":ceiling,"windows":report}


def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
    return bool(np.allclose(a,b,rtol=2e-9,atol=1e-13))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=run()
    if args.check:
        assert close(result,json.loads((HERE/"results.json").read_text(encoding="utf-8")))
        print("1428 two-phase saturation, reciprocal four probabilities and finite spectra verified")
    else:
        print(json.dumps(result,indent=2))
