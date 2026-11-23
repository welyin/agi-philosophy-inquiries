"""Test the joint discrete charges and a two-flavon lift through physical heavy Dirac ports."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def sign3(a):
    return (a+1)%3-1


def charge_checks():
    # Each tuple is (multiplicity, old F mod 3, chi, conventional lepton triality).
    entries=[]
    for t,f in zip([0,0,1],[2,0,1]):
        entries.extend([(6,t,-1,0),(3,t,-1,0),(1,t,-1,-1),(3,f,0,0),(2,f,0,1)])
    for f,c,l in [(0,1,-1),(2,1,-1),(0,-1,1),(6,-1,1),(0,-1,0)]:
        entries.append((1,f,c,l))
    pure=sum(n*l**3 for n,f,c,l in entries)
    mixed=[sum(n*sign3(f)**(3-k)*l**k for n,f,c,l in entries)%3 for k in [1,2]]
    assert pure==3 and mixed==[2,0]
    good_L=[j for j in itertools.combinations_with_replacement([-1,0,1],4)
            if (pure+sum(q**3 for q in j))%9==0]
    assert good_L==[(-1,-1,-1,0)]
    passes=[]
    old_axial=[]
    for f in itertools.product([-1,0,1],repeat=4):
        if sum(f)%3 or sum(q*q for q in f)%3:
            continue
        old_axial.append(f)
        j=good_L[0]
        if (2+sum(q*q*l for q,l in zip(f,j)))%3==0 and sum(q*l*l for q,l in zip(f,j))%3==0:
            passes.append(f)
    assert len(passes)==0
    # Mixed SM Lie-algebra coefficients and pure U(1) with each full family + right neutrino.
    d=np.array([6,3,3,2,1,1])
    y=np.array([1,-4,2,-3,6,0])   # 6Y in Q,u^c,d^c,L,e^c,J order
    q=np.array([1,-1,-1,-3,3,3])
    single=[sum(d*q),sum(d*q**3),sum(d*q*y**2),sum(d*q**2*y),
            2*q[0]+q[1]+q[2],3*q[0]+q[3]]
    assert not any(single)
    return {"pure_lepton_base_cubic":pure,"pure_F_F_L_and_F_L_L_residues_before_J":mixed,
            "J_lepton_charges_required":list(good_L[0]),
            "F3_assignments_passing_two_old_mixed_requirements":len(old_axial),
            "passing_all_required_F3_and_L_mixed_tests":len(passes),
            "one_family_integer_BminusL_mixed_and_pure":list(map(int,single))}


def realized(w):
    flav=read("../1685/inputs.json")
    eps=flav["scalar_magnitude_over_cutoff"]
    x=flav["real_X_over_cutoff"]
    s=eps*np.exp(1j*np.pi*flav["scalar_angle_over_pi"])
    phase=s/eps
    ell=np.array([2,0,1])
    row=np.diag(phase**ell)
    omega=phase**4
    legs=np.array([1,phase**2,w["neutral_modulus"]*phase**3*(3+2*omega)/math.sqrt(7)])
    a,b=legs*[0,1,1],legs*[1,3,1]
    major=np.outer(a,a)+w["source_r"]*np.exp(2j*np.pi/3)*np.outer(b,b)
    major*=math.sqrt(w["d21_eV2"])/np.linalg.svd(major,compute_uv=False)[1]
    target=(major.conj().T@major).conj()
    gram=row.conj().T@target@row
    ev,u=np.linalg.eigh(gram)
    nr=u[:,0]
    real=np.cross(nr.real,nr.imag)
    real/=np.sqrt(np.vdot(real, np.linalg.pinv(gram,rcond=1e-12)@real).real)
    rem=gram-np.outer(real,real)
    eig,vr=np.linalg.eigh(rem)
    z=math.sqrt(eig[-1])*vr[:,-1]
    np.testing.assert_allclose(np.outer(real,real)+np.outer(z,z.conj()),gram,atol=1e-17)
    assert max(abs(gram.imag).ravel())>1e-5
    # Preserve the two existing N heavy links, add a single gauge-allowed two-insertion H6-J3 edge.
    c=read("../1689/inputs.json")
    m,g,h=[c[i] for i in ["m_N0H0","coupling_to_s2","coupling_to_x"]]
    k=g/m
    heavy=np.array([[0,0,0,m,g*s**2],[0,0,h*s**6,0,h*x]],complex)
    wmix=s**6/x
    normalization=math.sqrt(1+abs(wmix)**2*(1+abs(k*s*s)**2))
    basis=np.eye(5,3,dtype=complex)
    basis[3,2]=k*s**2*wmix
    basis[4,2]=-wmix
    basis[:,2]/=normalization
    np.testing.assert_allclose(heavy@basis,0,atol=1e-15)
    np.testing.assert_allclose(basis.conj().T@basis,np.eye(3),atol=1e-15)
    cphase=wmix*s**2
    y=-normalization*z.imag/cphase.imag
    basic=normalization*z.real+cphase.real*y
    n0=np.array([.003,.012,.009])   # eV, independent real row couplings
    port=row@np.column_stack([real/math.sqrt(2),real/math.sqrt(2),basic,n0,s**2*(y+k*n0)])
    effective=port@basis
    np.testing.assert_allclose(effective@effective.conj().T,target,atol=1e-17)
    light=np.linalg.svd(effective,compute_uv=False)
    heavy_eV=np.linalg.svd(heavy,compute_uv=False)*c["sterile_unit_GeV"]*1e9
    rel=np.linalg.norm(port,2)**2/(heavy_eV[-1]**2-light[0]**2)
    full_dirac=np.vstack([port/(c["sterile_unit_GeV"]*1e9),heavy])
    no=np.zeros_like(full_dirac)
    symmetric=np.block([[no,full_dirac],[full_dirac.T,no]])
    adj=symmetric.conj().T
    for momentum in [1e-3,1,100]:
        lnv=np.linalg.solve(momentum**2*np.eye(len(adj))+adj@symmetric,adj)[:3,:3]
        assert np.linalg.norm(lnv)==0
    return {"R":round(w["R"],9),"target_sin2_2theta13_solar_atmos_in_same_charged_basis":w["three_angles"],
            "single_insertion_max_im_fixed_row_frame_eV2":float(f"{np.max(np.abs(gram.imag)):.10g}"),
            "full_two_insertion_oscillation_gram_error_eV2":float(f"{np.max(abs(effective@effective.conj().T-target)):.3e}"),
            "light_masses_eV":[round(mm,11) for mm in light],
            "tree_level_heavy_mixing_relative_bound":float(f"{rel:.10g}"),
            "phase_s8_over_x_degrees":round(float(np.angle(cphase,deg=True)),9),
            "light_right_null_column_norm_before_normalizing":round(normalization,11)}


def compute():
    return {"anomaly_tests":charge_checks(),
            "Dirac_two_insertion_through_H6_J_witnesses":[realized(w) for w in read("../1682/independent_results.json")["witnesses"]]}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=compute()
    if args.check:
        expected=read("independent_results.json")
        for a,b in zip(result["Dirac_two_insertion_through_H6_J_witnesses"],expected["Dirac_two_insertion_through_H6_J_witnesses"]):
            assert a["full_two_insertion_oscillation_gram_error_eV2"]<1e-16
            assert b["full_two_insertion_oscillation_gram_error_eV2"]<1e-16
            a["full_two_insertion_oscillation_gram_error_eV2"]=b["full_two_insertion_oscillation_gram_error_eV2"]
        assert result==expected
        print("PASS: pure lepton obstruction and exact full two-insertion Dirac lifts at both physical ports")
    else:
        print(json.dumps(result,indent=2,allow_nan=False))
