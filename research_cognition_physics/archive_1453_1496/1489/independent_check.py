"""Test the loop phase directly and restore a cube-root route with one real up tilt."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def rotation(i,j,angle,phase=0.):
    c,s=math.cos(angle),math.sin(angle)
    w=np.eye(3,dtype=complex)
    w[i,i]=w[j,j]=c
    w[i,j]=s*np.exp(-1j*phase)
    w[j,i]=-w[i,j].conj()
    return w


def verify():
    data=json.loads((HERE.parent/"1386/inputs.json").read_text(encoding="utf-8"))
    s13raw=data["quark_sines_and_phase"]
    a,b,c,delta=[x[0] for x in s13raw]
    V=rotation(1,2,math.asin(c)) @ rotation(0,2,math.asin(b),delta) @ rotation(0,1,math.asin(a))
    masses=np.array([i[0] for i in data["yukawa"]["down"]])
    masses/=masses[-1]
    original=(V*masses**2)@V.conj().T
    # Unobservable up mass-eigenstate phases are chosen BEFORE introducing the physical tilt.
    convention=np.diag([1.,1.,np.exp(1j*np.angle(original[0,2]))])
    def material(angle):
        w=rotation(0,2,math.radians(angle))
        vk=w@convention@V
        hd=(vk*masses**2)@vk.conj().T
        P=np.diag([1.,np.exp(1j*np.angle(hd[0,1])),1.])
        vd=P@vk
        hd=(vd*masses**2)@vd.conj().T
        ih=(vd*masses**-2)@vd.conj().T
        inv=-hd[0,1]*hd[2,0]*ih[1,2]
        assert abs(hd[0,1].imag)<1e-12 and abs(hd[0,2].imag)<1e-12
        expected=np.prod([masses[i]**2-masses[j]**2 for i,j in [(1,0),(2,0),(2,1)]])
        jcp=(vd[0,0]*vd[1,1]*vd[0,1].conj()*vd[1,0].conj()).imag
        assert np.isclose(inv.imag,hd[0,0].real*jcp*expected/np.prod(masses**2))
        return w,P,vd,hd,inv
    lo,hi=-.2,-.175
    phase=lambda t: np.angle(material(t)[-1],deg=True)
    assert phase(lo)>60>phase(hi)
    for _ in range(60):
        mid=(lo+hi)/2
        if phase(mid)>60:
            lo=mid
        else:
            hi=mid
    angle=(lo+hi)/2
    w,P,vd,hd,inv=material(angle)
    A=math.sqrt(hd[0,0].real); C=hd[0,1].real/A; E=hd[2,0].real/A
    F=math.sqrt(hd[2,2].real-E*E)
    Z=(hd[1,2]-C*E)/F
    # The sign of the real route is free. phase 60 deg is -120 deg modulo pi.
    assert abs((Z*np.exp(2j*math.pi/3)).imag) < 3e-14
    B=-math.sqrt(hd[1,1].real-C*C-abs(Z)**2)
    yd=np.array([[0,A,0],[B,C,Z],[0,E,F]])
    up=np.array([i[0] for i in data["yukawa"]["up"]])
    yu=w@np.diag(up)@w.T
    assert np.allclose(yd@yd.conj().T,hd,atol=1e-13)
    assert np.allclose(yu.imag,0) and np.linalg.det(yu).real>0 and np.linalg.det(yd).real>0
    # P acts only on row 2, which does not mix in the up block.
    vu=P@w@convention
    assert np.allclose(vu.conj().T@vd,V)
    assert np.allclose(vu@np.diag(up**2)@vu.conj().T,yu@yu.conj().T)
    raw_closed=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    bound=raw_closed["cos_delta_continuous_certificates"]["2"]["upper_loop_angle_deg"]
    central=abs(phase(0.))
    return {"original_central_phase_deg":central,"real_up_tilt_deg":angle,"new_loop_deg":phase(angle),
            "down_normalized_edges":[A,B,C,E,F],"normalized_Z":[Z.real,Z.imag],
            "maximum_CKM_matrix_error":float(np.max(abs(vu.conj().T@vd-V))),
            "all_four_rays_min_relative_any_complex_error_conditional_on_author_outer":
                math.sin(math.radians(60-bound)/4),
            "all_four_central_sharp_distance":math.sin(math.radians(60-central)/4)}


def comparable(a,b):
    return a.keys()==b.keys() and all(np.allclose(a[k],b[k],rtol=1e-10,atol=1e-14) for k in a)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    o=p.add_mutually_exclusive_group()
    o.add_argument("--write",action="store_true");o.add_argument("--check",action="store_true")
    result=verify(); target=HERE/"independent_results.json"
    if p.parse_args().write:
        with target.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False);out.write("\n")
    else:
        assert comparable(result,json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps(result,indent=2))
