"""Neutral color-singlet scalar flavor change forced by a right weak partner."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
baseline=json.loads((HERE/cfg["anchor"]).read_text(encoding="utf-8"))


def ckm(sines, phase):
    s12,s13,s23=sines
    c12,c13,c23=np.sqrt(1-np.square(sines))
    e=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13*e.conjugate()],
        [-s12*c23-c12*s23*s13*e,c12*c23-s12*s23*s13*e,s23*c13],
        [s12*s23-c12*c23*s13*e,-c12*s23-s12*c23*s13*e,c23*c13]])


def cross(mat):
    return float(np.linalg.norm(mat[:2,2])**2+np.linalg.norm(mat[2,:2])**2)


def unitary_from_column(v):
    # Choose two orthogonal vectors in the complement; no phase matching is required there.
    _,_,vh=np.linalg.svd(v.conj()[None,:],full_matrices=True)
    return np.column_stack((vh.conj().T[:,1:],v))


def evaluate():
    matrix=baseline["yukawa"]["up"]
    y=np.array([x[0] for x in matrix])
    angles=baseline["quark_sines_and_phase"]
    V=ckm([x[0] for x in angles[:3]],angles[-1][0])
    A=V.conj().T@np.diag(y)
    w=A.conj().T[:,2]
    form=np.diag(y*y)-2*np.outer(w,w.conj())
    eig,vec=np.linalg.eigh(form)
    sharp=float(np.vdot(w,w).real)
    sharp+=float(eig[0])
    R=unitary_from_column(vec[:,0])
    np.testing.assert_allclose(R.conj().T@R,np.eye(3),atol=1e-13)
    N=A@R
    np.testing.assert_allclose(cross(N),sharp,rtol=2e-11)
    # A down-aligned direction can contain all diagonal down masses with no bq decay.
    down=np.diag([x[0] for x in baseline["yukawa"]["down"]])
    pole=[down/math.sqrt(2),1j*down/math.sqrt(2),N/math.sqrt(2),1j*N/math.sqrt(2)]
    assert cross(pole[0])==0
    np.testing.assert_allclose(sum(map(cross,pole)),sharp,rtol=3e-11)
    rotation=np.linalg.qr(np.array([[1,2,-3,1],[2,1,4,-2],[4,2,1,-3],[1,3,2,5]],float))[0]
    rotated=np.einsum("st,tij->sij",rotation,pole)
    np.testing.assert_allclose(sum(map(cross,rotated)),sharp,rtol=3e-11)
    mult=cfg["widths"]
    top_lo,top_hi=y[2]-mult*matrix[2][1],y[2]+mult*matrix[2][1]
    charm_hi=y[1]+mult*matrix[1][1]
    def q(sign):
        a=angles[1][0]+sign*mult*angles[1][1]
        b=angles[2][0]+sign*mult*angles[2][1]
        return a*a+(1-a*a)*b*b
    qlo,qhi=q(-1),q(1)
    gap=top_lo**2*(1-2*qhi)-charm_hi**2
    assert gap>0
    correction=4*top_hi**2*charm_hi**2*qhi/gap
    box_floor=top_lo**2*qlo-correction
    # Test the conservative analytic lower bound at the spectrum/angle corners.
    from itertools import product
    for st,sc,s13,s23 in product((-1,1),repeat=4):
        masses=y.copy()
        masses[2]+=st*mult*matrix[2][1]
        masses[1]+=sc*mult*matrix[1][1]
        z=ckm([angles[0][0],angles[1][0]+s13*mult*angles[1][1],
               angles[2][0]+s23*mult*angles[2][1]],angles[-1][0])[:,2]*masses
        k=np.diag(masses*masses)-2*np.outer(z,z.conj())
        assert np.vdot(z,z).real+np.linalg.eigvalsh(k)[0]>=box_floor-1e-12
    rate_factor=3/(8*math.pi)
    phase_floor=(1-(cfg["kinematical_bottom_plus_light_mass_GeV_upper"]/
                     cfg["all_poles_mass_GeV_lower"])**2)**2
    physical=rate_factor*box_floor*phase_floor*cfg["partial_width_matching_floor"]
    h_max=cfg["ATLAS_total_width_upper_95_MeV"]/(cfg["h_GeV"]*1000)
    return {"central_sharp_Yukawa_squared":sharp,"Vcb_Vub_mix_squared_box":[qlo,qhi],
            "two_width_analytic_floor":box_floor,"schur_loss_upper":correction,
            "central_width_per_mass_floor":sharp*rate_factor,
            "two_width_massless_floor":rate_factor*box_floor,"finite_phase_floor":phase_floor,
            "matched_width_per_mass_floor":physical,"h125_share_upper_ATLAS":h_max/physical,
            "other_poles_width_per_mass_floor_ATLAS":physical-h_max,
            "h_aligned_zero_witness":cross(pole[0])}


def compare(a,b):
    assert a.keys()==b.keys()
    for key in a:
        np.testing.assert_allclose(a[key],b[key],rtol=2e-9,atol=1e-13)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=evaluate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        compare(out,json.loads(path.read_text(encoding="utf-8")))
    print("1549: sharp right-handed minimization, finite width floor and 125-GeV coverage verified")
