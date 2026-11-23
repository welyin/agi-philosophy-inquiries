"""Independent PMNS-determinant test at the electron-fraction crossing."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8"))


def rotation(i,j,angle,phase=0.):
    r=np.eye(3,dtype=complex)
    r[i,i]=r[j,j]=math.cos(angle)
    r[i,j]=math.sin(angle)*np.exp(-1j*phase)
    r[j,i]=-r[i,j].conjugate()
    return r


def audit():
    cfg=read("inputs.json")
    s=cfg["solar_sin2"][0]
    amp=read(cfg["reactor_2theta"])["angles"]["reactor_amplitude"][0]
    z=.5*(1-math.sqrt(1-amp))
    mass=read(cfg["splittings"])["juno_207_days"]
    offset=np.array([0.,mass["solar_21"][0],mass["normal_atmospheric_31"][0]])
    t12,t13=np.arcsin(np.sqrt([s,z]))
    aa=math.cos(2*t12)*math.cos(2*t13)
    bb=(1-3*z)*math.sin(2*t12)*math.sqrt(z)/2
    half=.5*bb/math.hypot(aa,bb)
    def branch(t,delta):
        t23=math.asin(math.sqrt(t))
        unitary=rotation(1,2,t23)@rotation(0,2,t13,delta)@rotation(0,1,t12)
        mod=abs(unitary)**2
        e=mod[0]
        singlet=(1-mod[:,0])/2
        lo,hi=0.,.01
        for _ in range(90):
            x=(lo+hi)/2
            tested=(mod@(offset+x))/(sum(offset)+3*x)
            if tested[0]<singlet[0]:
                lo=x
            else:
                hi=x
        x=(lo+hi)/2
        difference=mod@(offset+x)/(sum(offset)+3*x)-singlet
        det=float(np.linalg.det(mod))
        direct=aa*math.cos(2*t23)-bb*math.sin(2*t23)*math.cos(delta)
        factor=x/(sum(offset)+3*x)/(e[2]-e[1])
        np.testing.assert_allclose([difference[0],difference[1],difference[2],det],
                                   [0.,factor*det,-factor*det,direct],atol=5e-15)
        return {"sin2_theta23":t,"delta_degrees":math.degrees(delta),
                "det_mod_squared":det,"triplet_minus_singlet_mu":float(difference[1]),
                "sum_eV":float(np.sqrt(x+offset).sum()),"lightest_eV":math.sqrt(x),
                "beta_eV":math.sqrt(x+e@offset)}
    given=.55
    t23=math.asin(math.sqrt(given))
    cosdelta=aa/(bb*math.tan(2*t23))
    blind=branch(given,math.acos(cosdelta))
    visible=branch(given,0.)
    assert abs(blind["triplet_minus_singlet_mu"])<1e-14
    assert abs(visible["triplet_minus_singlet_mu"])>.005
    limits=cfg["cosmology"]["sum_limits_eV"]
    assert limits[0]<blind["sum_eV"]<limits[1]
    assert blind["beta_eV"]<cfg["beta"]["m_beta_limit_eV"]
    return {"central_octant_window_allowing_full_degeneracy":[.5-half,.5+half],
            "non_mu_tau_symmetric_exact_overlap":blind,"different_phase_same_e":visible,
            "gap_after_20pct_symmetric_mu_tau_mistag":.6*abs(visible["triplet_minus_singlet_mu"])}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
