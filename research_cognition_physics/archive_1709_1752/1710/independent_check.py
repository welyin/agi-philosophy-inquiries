"""Resolve the anomalous adjoint phase and the split-heavy-multiplet soft response."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"independent_results.json"


def prior(round_number, filename="results.json"):
    return json.loads((ROOT.parent/str(round_number)/filename).read_bytes())


def loop(a2,b2,p2,n):
    # a is the fermion mass, b the scalar mass; real below the fermion/scalar cut.
    q,w=np.polynomial.legendre.leggauss(n)
    x=(q+1)/2
    den=x*a2+(1-x)*b2-x*(1-x)*p2
    assert min(den)>0
    return np.dot(w/2,(1-x)/den)


def chirality_slopes(mass,soft,p2,n):
    return np.array([2*m**2*loop(m**2,m**2+soft**2,p2,n) for m in mass])


def mixed_threshold_check(mass, scalar, deriv, p2, eta, n):
    # Direct unexpanded scalar mass-matrix loop; eta=F_S/S in GeV.
    mu2=np.diag(scalar**2)
    B=eta*deriv
    scalar_matrix=np.block([[mu2,B],[B.conj().T,mu2]])
    eig,V=np.linalg.eigh(scalar_matrix)
    assert eig.min()>0
    q,w=np.polynomial.legendre.leggauss(n)
    x=(q+1)/2
    ans=0j
    for i,m in enumerate(mass):
        den=x[:,None]*m**2+(1-x[:,None])*eig-x[:,None]*(1-x[:,None])*p2
        assert den.min()>0
        # The arbitrary common log constant drops out in the off-diagonal matrix element.
        kernel=(w[:,None]/2*np.log(den/m**2)).sum(axis=0)
        ans+=2*m*np.dot(V[i]*V[i+2].conj(),kernel)
    return ans


def compute():
    raw=prior(1698)["shared_witness"]
    base=prior(1685,"inputs.json")
    old=prior(1685)
    theta=prior(1389)["conditional_edm"]["theta_limit_rad"]
    weak=prior(1369,"inputs.json")
    soft=raw["soft_GeV"]
    M3=soft
    # The rematched S vev and cutoff are from 1698. A3=S^3/Lambda^2.
    vev=raw["shared_S_X_B_vev_modulus_GeV"]
    cutoff=raw["common_flavor_cutoff_GeV"]
    mix=np.array(raw["color_A3_mixing_abs_yukawas"])*vev**3/cutoff**2
    top=1000*base["cutoff_TeV"]*base["colored_M_T_over_cutoff"]
    matrix=np.array([[0,mix[0]],[-mix[1],top]])
    u,m,vt=np.linalg.svd(matrix)
    deriv=u.T@(3*np.array([[0,mix[0]],[-mix[1],0]]))@vt.T
    slopes=np.diag(deriv)/m
    np.testing.assert_allclose(m,raw["colored_fermions_GeV"],rtol=1e-9)
    np.testing.assert_allclose(slopes,old["individual_log_mass_radial_slopes"],atol=5e-8)
    mu=np.hypot(m,soft)
    np.testing.assert_allclose(mu,raw["colored_scalars_GeV"],rtol=1e-8)
    positivity=1/np.linalg.norm(deriv/(mu[:,None]*mu[None,:]),ord=2)
    batches={}
    for label,p2 in [("zero_momentum",0.),("unperturbed_gluino_shell",M3**2)]:
        h=chirality_slopes(m,soft,p2,192)
        np.testing.assert_allclose(h,chirality_slopes(m,soft,p2,96),rtol=6e-9)
        if p2==0:
            r=(mu/m)**2
            analytic=2*(r-1-np.log(r))/(r-1)**2
            np.testing.assert_allclose(h,analytic,rtol=5e-12)
        coef=float(np.dot(h,slopes))
        btest=.25j
        check=mixed_threshold_check(m,mu,deriv,p2,btest,96)/btest
        assert abs(check-coef)<2e-8
        batches[label]={"split_kernel_per_heavy_and_light":h.tolist(),
                        "individual_weighted_slopes":(h*slopes).tolist(),
                        "gluino_mass_insertion_coefficient":coef,
                        "finite_B_cross_check_at_imaginary_F_over_S_0p25_GeV":float(check.real)}
    fermion_cp=2*(4/3)*slopes*M3**2*np.array([loop(M3**2,s**2,0.,96) for s in mu])
    # The witness contains just the light colored Dirac triplet below 30 TeV. DR matching as in 1698.
    inv_alpha=4*math.pi/weak["g3"][0]**2+7*math.log(soft/weak["mu_ref_GeV"])/(2*math.pi)
    inv_alpha-=((2/3)*math.log(soft/m[-1])/(2*math.pi)+3/(12*math.pi))
    g2=4*math.pi/inv_alpha
    loopfac=g2/(16*math.pi**2)
    fbound=M3*theta/(3*loopfac*batches["zero_momentum"]["gluino_mass_insertion_coefficient"])
    # An independent flavour-neutral external SUSY-breaking spurion leaves the scalar vacuum unchanged.
    I=vev/cutoff*(vev/cutoff)**4*cmath.exp(-2j*math.pi/3)
    omega=1e-6
    coeff=math.tan(omega)/I.imag
    h=(1+coeff*(I-I.real))/math.hypot(1,coeff*I.imag)
    assert abs(h-cmath.exp(1j*omega))<1e-15
    return {"common_S_GeV":vev,"dimless_X_S4_invariant":[I.real,I.imag],
            "central_soft_and_abs_M3_GeV":M3,"g3_squared_MSplus_as_DR":g2,
            "heavy_light_fermions_GeV":m.tolist(),"soft_split_scalar_masses_GeV":mu.tolist(),
            "heavy_light_radial_slopes":slopes.tolist(),"finite_B_tachyon_free_abs_F_over_S_strict_upper_GeV":positivity,
            "adjoint_coefficient":3,"single_theta_working_limit_rad":theta,
            "soft_mass_phase_distance_to_2pi_over_3_lattice_at_tree_rad":theta/3,
            "same_modulus_gluino_tree_control":{"abs_h":abs(h),"phase_rad":cmath.phase(h),
                "real_coupling":coeff,"absolute_theta_shift_rad":3*abs(cmath.phase(h))},
            "crossed_heavy_multiplet_linear_response":batches,
            "threshold_conjugate_dirac_mass_phase_abs_coefficients_at_p0":fermion_cp.tolist(),
            "sum_dirac_phase_abs_coefficients":float(fermion_cp.sum()),
            "gluino_only_p0_one_source_abs_Im_F_over_S_GeV_budget":fbound,
            "all_orders_mark":"The F deformation is a response away from 1698's B=0 boundary, not a fit or full theta matching."}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.check:
        assert result==json.loads(OUT.read_bytes())
        print("1710 independent adjoint and colored-soft threshold checks passed")
    else:
        OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

