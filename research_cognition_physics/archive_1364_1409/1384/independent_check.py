"""Exact folded-ray gap and an explicitly attaining complex-unitary readout change."""
import argparse
import json
from math import asin, acos, cos, degrees, pi, sin, sqrt
from pathlib import Path

import numpy as np


def matrix(a,b,d,phase):
    c13=sqrt(1-b*b)
    s23=d/c13
    c23=sqrt(1-s23*s23)
    c12=a/c13
    s12=sqrt(1-c12*c12)
    e=np.exp(1j*phase)
    return np.array([[a,s12*c13,b/e],
        [-s12*c23-c12*s23*b*e,c12*c23-s12*s23*b*e,d],
        [s12*s23-c12*c23*b*e,-c12*s23-s12*c23*b*e,c23*c13]])


def rotate_ray(v,idx,target):
    # v is a unit vector, e the flavor readout. Tilt e in span{e,v}
    # to the smallest projective distance realising |<e',v>| = target.
    e=np.eye(3,dtype=complex)[:,idx]
    v=v*np.exp(-1j*np.angle(v[idx]))
    overlap=np.vdot(e,v).real
    w=(v-e*overlap)/sqrt(1-overlap**2)
    h=acos(target)-acos(overlap)
    rot=np.eye(3,dtype=complex)+(cos(h)-1)*(np.outer(e,e.conj())+np.outer(w,w.conj()))
    rot+=sin(h)*(np.outer(e,w.conj())-np.outer(w,e.conj()))
    assert np.allclose(rot.conj().T@rot,np.eye(3))
    assert abs(abs(np.vdot(rot@e,v))-target)<1e-14
    return h,rot


def compute_window(n):
    a=.97373+n*.00031
    d=(.0351+n*(.0019+.0018+.0008))/.80
    b=(.083+n*.008)*d
    mb=matrix(a,b,d,pi)
    assert np.allclose(mb@mb.conj().T,np.eye(3))
    u=sqrt(1-a*a-b*b)
    z=sqrt(b*b+d*d)
    c=(u*sqrt(1-z*z)-a*b*d)/(1-b*b)
    assert abs(c-abs(mb[1,0]))<1e-14
    assert abs(sqrt(1-abs(mb[2,2])**2)-z)<1e-14
    gaps=[4*z*z*(1-z*z),4*c*c*(1-c*c)]
    rrange=[.0851-n*.0024,.0851+n*.0024]
    rlo,rhi=(.5*asin(sqrt(r)) for r in rrange)
    # theta=asin sqrt(min(p,1-p)); these are the exact two neighboring distances.
    left=rlo-asin(z)
    right=asin(c)-rhi
    assert 0<right<left and gaps[0]<rrange[0]<rrange[1]<gaps[1]
    target=sin(rhi)
    shift,rot=rotate_ray(mb[:,0],1,target)
    out=rot.conj().T@mb
    assert abs(shift-right)<3e-14
    assert abs(4*abs(out[1,0])**2*(1-abs(out[1,0])**2)-rrange[1])<3e-14
    # Additional all-row test using only ud and cs and the low reactor branch.
    cs=.9589-n*(.0023+.0040+.0015+.0005+.0095)
    a_low=.97373-n*.00031
    row_floor=a_low+cs-1
    solar_low=.3092-n*.0087
    e_max=sqrt((1-sin(rlo)**2)*(1-solar_low))
    assert row_floor>e_max
    # Saturation of |U33| >= a+cs-1 by a real rotation in the plane
    # of e3 and (sqrt(1-a),sqrt(1-cs),0).
    axis=np.array([sqrt(1-a_low),sqrt(1-cs)])
    norm=np.linalg.norm(axis)
    direction=axis/norm
    sn=sqrt(1-row_floor**2)
    att=np.eye(3)
    att[:2,:2]-=(1-row_floor)*np.outer(direction,direction)
    att[:2,2]=sn*direction
    att[2,:2]=-sn*direction
    att[2,2]=row_floor
    assert np.allclose(att@att.T,np.eye(3))
    assert np.allclose(np.diag(att),[a_low,cs,row_floor])
    return {"quark_folded_modulus_band_edges":[z,c],
        "quark_reactor_amplitude_band_edges":gaps,"reactor_amplitude_box":rrange,
        "folded_angle_gaps_deg":[degrees(left),degrees(right)],
        "min_attaining_single_basis_change_deg":degrees(shift),
        "min_attaining_single_projector_trace_distance":sin(shift),
        "source_quark_rephased_cell":float(abs(mb[1,0])),
        "new_lepton_cell_absolute":float(abs(out[1,0])),
        "solar_row_anchor_floor":row_floor,"solar_electron_row_upper_max":e_max,
        "additional_row_ray_angle_lower_deg":degrees(acos(e_max)-acos(row_floor))}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--check",action="store_true")
    out=Path(__file__).with_name("independent_results.json")
    result=rounded({str(n):compute_window(n) for n in (1,2)})
    if ap.parse_args().write:
        with out.open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2,ensure_ascii=False)
            fp.write("\n")
    else:
        assert result==json.loads(out.read_text(encoding="utf-8"))
    print("1384 exact both-branch gaps, attaining ray rotations and solar row certificates verified.")
