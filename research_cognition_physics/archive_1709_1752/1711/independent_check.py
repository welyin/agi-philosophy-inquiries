"""Holomorphic flavour orbits, the unique two-leg form, and a physical mass-ratio test."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/"independent_results.json"


def read(folder,name="results.json"):
    return json.loads((HERE.parent/str(folder)/name).read_bytes())


def bridge(up,down):
    return down.T@np.linalg.solve(up,down)


def first_second_probability(up,down):
    # Eigenvectors of the two Hermitian Gram matrices with the canonical kinetic metric.
    _,uv=np.linalg.eigh(up@up.conj().T)
    _,dv=np.linalg.eigh(down@down.conj().T)
    return float(abs(np.vdot(uv[:,0],dv[:,1]))**2)


def spectrum(C):
    return np.linalg.svd(C,compute_uv=False)[::-1]


def shape(s):
    return (s[2]**2-s[0]**2)/(s[1]**2-s[0]**2)


def lift(c,phase):
    return np.diag(np.exp(1j*np.array([phase,-phase,0.])))@c


def split(C,v):
    w=C@v
    return C-np.outer(w,w)/(v@w)


def compute():
    old=runpy.run_path(str(HERE/"../1702/check.py"))
    gamma=read(1709)["same_1708_point"]["two_port_Cabibbo"]["quarks"]["gamma_deg"]
    y_u=old["central"]("up")
    y_d=old["central"]("down")
    y_d[0]-=old["low"]["yukawa"]["down"][0][1]
    V=old["low_ckm"](gamma,old["angles"]("direct",-1.))
    A=np.diag(y_u)
    D=V@np.diag(y_d)
    aligned=np.diag(y_d)
    # Right 5bar-complexified change of coordinates relates aligned D to the physical one.
    S=np.linalg.solve(aligned,D)
    np.testing.assert_allclose(np.linalg.det(S),1,atol=1e-12)
    np.testing.assert_allclose(aligned@S,D,atol=1e-14)
    np.testing.assert_allclose(spectrum(aligned),spectrum(D),rtol=1e-12)
    first=first_second_probability(A,D)
    assert abs(math.sqrt(first)-.22333)<1e-11
    assert first_second_probability(A,aligned)==0
    # Test the bilinear transformation and the labelled vector Schur complement under general GL(3)^2.
    G=np.array([[1,.2j,0],[0,1,.23],[.07,0,1]],complex)
    H=np.array([[1,-.4,.12j],[0,1,.11j],[-.2j,.05,1]],complex)
    C=bridge(A,D)
    transformed=bridge(G.T@A@G,G.T@D@H)
    np.testing.assert_allclose(transformed,H.T@C@H,atol=2e-14)
    vectors=read(1688,"inputs.json")
    a=np.array(vectors["source_a"],float)
    b=np.array(vectors["source_b"],float)
    v=np.cross(a,b)
    v/=np.linalg.norm(v)
    v_h=np.linalg.solve(H,v)
    reduced=split(C,v)
    # Singlet channels N plus an extra singlet L with no L^2 term implement the same Schur complement.
    w=D@v
    heavy=np.block([[A,w[:,None]],[w[None,:],np.zeros((1,1))]])
    sources=np.vstack((D,np.zeros((1,3))))
    block_form=sources.T@np.linalg.solve(heavy,sources)
    np.testing.assert_allclose(block_form,reduced,atol=2e-15)
    np.testing.assert_allclose(split(transformed,v_h),H.T@reduced@H,atol=2e-14)
    np.testing.assert_allclose(split(C,v*cmath_phase(.831)),reduced,atol=2e-14)
    assert np.linalg.norm(reduced@v)<1e-14
    sr=spectrum(reduced)
    assert sr[0]<sr[1]*1e-10 and sr[1]>1e-9
    # The archived 1709 complex source uses rows with flavour phases. A raw null vector transforms covariantly.
    phi=math.pi*read(1685,"inputs.json")["scalar_angle_over_pi"]
    phases=np.exp(1j*phi*np.array(vectors["su5_generation_charges_5bar"]))
    dressed=np.diag(phases)@np.column_stack((b,a))
    vn=v/phases
    assert np.linalg.norm(dressed.T@vn)<1e-14
    # Determinants, singular values, and every CKM quartet are fixed under this det-one rephasing.
    inputs=read(1668,"inputs.json")
    solar=inputs["juno_NO_dm21"][0]
    atm=inputs["juno_NO_dm31"][0]
    def evaluate(t):
        d=lift(D,t)
        np.testing.assert_allclose(np.linalg.det(d),np.linalg.det(D),rtol=1e-12,atol=1e-22)
        np.testing.assert_allclose(first_second_probability(A,d),first,atol=5e-13)
        return spectrum(bridge(A,d))
    lo,hi=0.,.05
    low,high=shape(evaluate(lo)),shape(evaluate(hi))
    assert low<atm/solar<high
    for _ in range(60):
        mid=(lo+hi)/2
        if shape(evaluate(mid))<atm/solar:lo=mid
        else:hi=mid
    t=(lo+hi)/2
    s=evaluate(t)
    scale=math.sqrt(solar/(s[1]**2-s[0]**2))
    masses=s*scale
    np.testing.assert_allclose(masses[2]**2-masses[0]**2,atm,rtol=1e-12)
    # For any three row phases the largest singular value is <= K by the triangle inequality.
    # The second <= L after subtracting the largest norm individual rank-one term.
    row_weights=np.sum(abs(D)**2,axis=1)/y_u
    K=sum(row_weights)
    L=K-max(row_weights)
    det=np.prod(y_d)**2/np.prod(y_u)
    np.testing.assert_allclose(np.linalg.det(C),det,rtol=1e-12,atol=1e-25)
    floor_solar=math.sqrt(solar)*det/(K*L**2)
    # m1*sqrt(m1^2+d21)*sqrt(m1^2+d31)=scale^3*|det C| >= rhs.
    rhs=(math.sqrt(atm)/K)**3*det
    a0,b0=0.,.1
    for _ in range(75):
        x=(a0+b0)/2
        if x*math.sqrt(x*x+solar)*math.sqrt(x*x+atm)<rhs:a0=x
        else:b0=x
    assert masses[0]>max(a0,floor_solar)
    return {"quark_seed_scale":"MZ", "gamma_deg":gamma,"original_common_CKM_modulus":math.sqrt(first),
            "same_determinants_and_singular_values_control_modulus":0.,
            "SL_right_change_of_basis_det_error":float(abs(np.linalg.det(S)-1)),
            "SL_right_change_of_basis_nonunitarity_norm":float(np.linalg.norm(S.conj().T@S-np.eye(3),2)),
            "symmetric_open_tensor_det_abs":det,"source_row_weights":row_weights.tolist(),
            "extra_row_phase_neutrino_gap_ratio":{"t0":low,"t0p05":high},
            "central_two_frequency_calibration":{"target_ratio":atm/solar,"det_one_row_phase_rad":t,
                "common_normalization_eV":scale,"three_masses_eV":masses.tolist()},
            "light_mass_normal_ordering_floor_all_row_phases_eV":{"solar_only":floor_solar,
                "with_central_both_gaps":max(a0,floor_solar)},
            "extra_projective_label_reduction":{"illustrative_unphased_v":v.tolist(),
                "reduced_singular_values":sr.tolist(),"null_norm":float(np.linalg.norm(reduced@v)),
                "four_singlet_block_matching_error":float(np.linalg.norm(block_form-reduced)),
                "released_phased_source_transpose_null_residual":float(np.linalg.norm(dressed.T@vn))}}


def cmath_phase(x):
    return complex(math.cos(x),math.sin(x))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.check:
        assert result==json.loads(OUT.read_bytes())
        print("1711 independent flavour orbit, rank and mass-spectrum checks passed")
    else:
        OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

