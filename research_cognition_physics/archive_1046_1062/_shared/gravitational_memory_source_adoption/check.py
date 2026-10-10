"""Count-zero checks of the adopted asymptotic memory/source dictionary.

No new weak-decay calculation, detector model, or gravity simulation.
Default mode recomputes and compares without writing.
"""
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
import argparse
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parents[1]

def calculate():
    old=ARCHIVE/"1056/results.json"
    prior=json.loads(old.read_text(encoding="utf-8"))
    cb=float(prior["quadrature_diagnostics_not_proof"]["72"]["dalitz"])
    pb=Fraction(prior["exact_certificate"]["record_probability"])
    # Exact energy fractions, already derived from the same three-body kernel.
    energy=2*Fraction(7,20)+Fraction(3,10)
    spin=2*Fraction(3,20)-Fraction(3,10)
    assert energy==1 and spin==0

    u,wu=np.polynomial.legendre.leggauss(32)
    phi=2*np.pi*np.arange(48)/48
    uu,pp=np.meshgrid(u,phi,indexing="ij")
    nn=np.stack([np.sqrt(1-uu**2)*np.cos(pp),
                 np.sqrt(1-uu**2)*np.sin(pp),uu],axis=-1)
    weights=np.broadcast_to(wu[:,None]*2*np.pi/48,uu.shape)
    eye=np.eye(3)
    proj=eye-nn[..., :, None]*nn[...,None,:]
    aa=np.zeros((3,3));aa[0,2]=aa[2,0]=.5
    f=nn[...,0]*nn[...,2]
    # Tangential trace-free Hessian on S^2 of f=n^T A n:
    # H[f] = 2 P A P + f P, since Tr A=0 and Delta f=-6f.
    hf=2*np.einsum("...ij,jk,...kl->...il",proj,aa,proj)+f[...,None,None]*proj
    scalar_norm=float(np.sum(weights*f*f))
    tensor_norm=float(np.sum(weights*np.einsum("...ij,...ij->...",hf,hf)))
    mean=float(np.sum(weights*f))
    transverse=float(np.max(abs(np.einsum("...ij,...j->...i",hf,nn))))
    trace=float(np.max(abs(np.trace(hf,axis1=-2,axis2=-1))))
    assert abs(mean)<1e-14
    assert abs(scalar_norm-4*np.pi/15)<1e-13
    assert abs(tensor_norm-16*np.pi/5)<1e-12
    assert transverse<1e-14 and trace<1e-14

    # Independent retarded-kernel normalization. A single null ray along z
    # is a Green-kernel diagnostic, not a complete conserved decay.
    # f_z=u^2-1/3: H_theta,theta=1-u^2, H_phi,phi=-(1-u^2).
    # Physical shear coefficient C_theta,theta=2GE(1+u), C_phi,phi=-C_theta,theta.
    # Set GE=1 only for this dimensionless normalization test.
    kernel_integral=float(2*np.pi*np.dot(wu,4*(1+u)*(1-u*u)))
    expected=16*np.pi*Fraction(2,3)
    assert abs(kernel_integral-expected)<1e-12
    # Polarization contraction and 1/2 conversion to fractional displacement
    # are kept separate; C is the physical metric shear, not canonically scaled h.
    l2_projection_factor=16*np.pi/tensor_norm
    assert abs(l2_projection_factor-5)<1e-12
    assert cb < -1/4500
    results={
        "adoption_count":0,
        "all_checks_passed":True,
        "scope":"massless-daughter leading weak-decay branch; adopted leading-G asymptotic retarded memory",
        "sphere_f_mean":mean,
        "sphere_f_squared":scalar_norm,
        "sphere_H_squared":tensor_norm,
        "sphere_H_transversality_residual":transverse,
        "sphere_H_trace_residual":trace,
        "retarded_kernel_test_integral_GE1":kernel_integral,
        "retarded_kernel_expected_GE1":float(expected),
        "l2_projection_coefficient_per_G_Q":l2_projection_factor,
        "joint_source_Cb_over_m_from_1056":cb,
        "joint_shear_functional_over_Gm":16*np.pi*cb,
        "joint_shear_l2_projection_over_Gm":5*cb,
        "conditional_record_probability":str(pb),
        "total_energy_fraction":str(energy),
        "total_spin_energy_fraction":str(spin),
        "unconditional_l2_memory_mean":0,
        "finite_radius_or_time_error_certified":False,
        "actual_gravity_instrument_certified":False,
        "full_coherent_daughter_state_certified":False,
        "nonzero_daughter_mass_error_certified":False,
        "new_empirical_data":False,
        "historical_sha256":{
            str(p.relative_to(ARCHIVE)):sha256(p.read_bytes()).hexdigest()
            for p in [ARCHIVE/"1056/proof.md",old,ARCHIVE/"research_note_1056.md"]
        },
    }
    return results

def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-12,abs_tol=2e-13),(path,a,b)
    else: assert a==b,(path,a,b)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--save",action="store_true")
    args=parser.parse_args()
    result=calculate()
    target=HERE/"results.json"
    if args.save:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        compare(json.loads(target.read_text(encoding="utf-8")),result)
    print(json.dumps(result,ensure_ascii=False))

