"""960: actual dimer spectral response matched to a declared Maxwell dipole EFT.

One leading-order stationary interface, not a real-time QED error certificate.
NumPy quadrature verifies formulas; analytic bounds are stated in note 960.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent; STAGE=HERE.parent; ROOT=HERE.parents[2]
TARGET=HERE/"material_dispersion_bridge_results.json"

def read(p): return json.loads(p.read_text("utf-8-sig"))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a): return float(np.linalg.norm(a,2))

def load956():
    path=STAGE/"956/native_material_interface.py"
    spec=importlib.util.spec_from_file_location("material956",path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def f(x): return np.exp(-x)*(1+x+x*x)
def fp(x): return np.exp(-x)*x*(1-x)

def integrate(fun,a,b,n=256):
    x,w=np.polynomial.legendre.leggauss(n)
    return float((b-a)/2*np.dot(w,fun(a+(b-a)*(x+1)/2)))

def response(r,n=256):
    """F(r), r F'(r); differentiate at fixed physical frequency."""
    if r<=1:
        def energy(t):
            y=np.tan(t);return np.cos(t)**2*f(r*y)**2
        def source(t):
            y=np.tan(t);x=r*y
            return np.cos(t)**2*2*f(x)*x*fp(x)
        return (4/math.pi*integrate(energy,0,math.pi/2,n),
                4/math.pi*integrate(source,0,math.pi/2,n))
    # x=r*y change of variable avoids a shrinking boundary layer.
    # Omitted x>40 tail is exponentially small, separately bounded below.
    return (4/(math.pi*r)*integrate(
                lambda x:f(x)**2/(1+(x/r)**2)**2,0,40,n),
            4/(math.pi*r)*integrate(
                lambda x:2*f(x)*x*fp(x)/(1+(x/r)**2)**2,0,40,n))

def run():
    old=read(STAGE/"958/capacitive_material_write_results.json")
    model=load956();U=1.;v=.1;kappa=.2
    m=model.material(U,v);H=m["H"];d=m["occupancy"];Delta=U+m["J"]
    cs=[model.annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[i for i in range(16) if i.bit_count()==2]]
    Q=sector.T@((ns[0]+ns[1]-ns[2]-ns[3])/2)@sector
    val,V=np.linalg.eigh(H);g=V[:,0];gap=val[1:]-val[0]
    weights=np.abs(V[:,1:].T@Q@g)**2
    spectral_rows=[]
    for xi in (0.,.1,Delta,10*Delta):
        spectral=float(np.sum(2*gap*weights/(gap*gap+xi*xi)))
        formula=2*Delta*d/(Delta*Delta+xi*xi)
        spectral_rows.append(dict(xi=xi,spectral=spectral,formula=formula,
                                  error=abs(spectral-formula)))
    spectral_error=max(r["error"] for r in spectral_rows)
    assert spectral_error<1e-13
    dark_error=norm(Q@m["W"]@(np.eye(4)-m["ps"]))
    assert dark_error<1e-13
    bright_error=norm((H-U*np.eye(6))@Q@g)
    assert bright_error<1e-13
    chi2=kappa*kappa*d*d/(2*Delta)
    # Independent ordinary second-order perturbation sum on actual 36 states.
    H0=np.kron(H,np.eye(6))+np.kron(np.eye(6),H)
    E,Z=np.linalg.eigh(H0);gg=Z[:,0];Vint=np.kron(Q,Q)
    coeff=float(np.sum(np.abs(Z[:,1:].T@Vint@gg)**2/(E[1:]-E[0])))
    chi_spectral=kappa*kappa*coeff
    assert abs(chi_spectral-chi2)<1e-15
    assert abs(chi2-old["exact_interface"]["second_order_chi"])<1e-15
    # Integral y=tan(t): alpha^2 dxi -> 4 d^2/Delta cos^2(t) dt.
    chi_integral=kappa*kappa/(2*math.pi)*4*d*d/Delta*integrate(
        lambda t:np.cos(t)**2,0,math.pi/2)
    assert abs(chi_integral-chi2)<1e-15
    B0=(1+3/math.e)/2; B1=6/math.e; BF=B0+B1/6
    rows=[]
    for r in (0.,.01,1.,100.,10000.):
        a,b=response(r,256);a2,b2=response(r,512)
        diff=max(abs(a-a2),abs(b-b2))
        assert diff<3e-10,(r,diff)
        assert abs(a-1)<=B0*r*r+1e-12
        assert abs(b)<=B1*r*r+1e-12
        # F_R=d chi/dR at fixed U,v,C,a,c, including kappa~R^-3.
        force_ratio=a-b/6
        assert abs(force_ratio-1)<=BF*r*r+1e-12
        if r>0:
            h=r*2e-4
            fm2=response(r-2*h,512)[0];fm1=response(r-h,512)[0]
            fp1=response(r+h,512)[0];fp2=response(r+2*h,512)[0]
            derivative=r*(fm2-8*fm1+8*fp1-fp2)/(12*h)
            deriv_error=abs(derivative-b)
            assert deriv_error<2e-9
        else:deriv_error=0.
        rows.append(dict(r=r,energy_ratio=a,logarithmic_kernel_derivative=b,
            force_ratio=force_ratio,quadrature_difference=diff,
            independent_derivative_error=deriv_error))
    # Check far-zone coefficient with exact elementary exponential moments.
    coefficients=[1,2,3,2,1] # (1+x+x^2)^2
    moment=sum(c*math.factorial(n)/2**(n+1) for n,c in enumerate(coefficients))
    assert abs(moment-13/4)<1e-15
    far_coefficient=math.pi*rows[-1]["r"]*rows[-1]["energy_ratio"]
    assert abs(far_coefficient-13)<2e-6
    # Exponential quadrature tail: int_L^infty x^n exp(-2x) dx exactly.
    def tail_moment(n,L):
        return math.exp(-2*L)*sum(math.factorial(n)/math.factorial(j)*
                                      L**j/2**(n-j+1) for j in range(n+1))
    far_quad_tail_F=4/math.pi*sum(c*tail_moment(n,40)
                                  for n,c in enumerate(coefficients))
    # For x>=40, |2 f x fp|=2 exp(-2x)(x^5-x^2)<=2 exp(-2x)x^5.
    far_quad_tail_derivative=8/math.pi*tail_moment(5,40)
    assert max(far_quad_tail_F,far_quad_tail_derivative)<1e-25
    # Domain certificate, not a fitted collection of exact-UV modes.
    rmax=.01;Y=100.
    energy_tail=12/(math.pi*math.e**2*Y**3)
    derivative_tail=24/(math.pi*math.e)*rmax*rmax/Y
    force_tail=energy_tail+derivative_tail/6
    # Geometry inherited only as a comparison, not silently as point dipoles.
    geom=old["sources"]["extra_static_geometry"]
    R=geom["R"];a=geom["a"];C=geom["coulomb_coefficient"]
    k_site=2*C*(1/R-1/math.sqrt(R*R+a*a));k_dip=C*a*a/R**3
    geometry_energy_ratio=(k_dip/k_site)**2
    assert abs(k_site-kappa)<1e-14 and geometry_energy_ratio>1.15
    # General finite-size leading-static bound via alternating Taylor remainder.
    z=(a/R)**2;site_over_dip=k_site/k_dip
    assert 1-3*z/4<=site_over_dip<=1
    chi_exact=old["exact_interface"]["conditional_energy_chi"]
    delta_chi=chi_exact-chi2
    source_files=[Path(__file__),STAGE/"956/native_material_interface.py",
        STAGE/"958/capacitive_material_write_results.json",
        STAGE/"959/gauss_material_parent_results.json",
        HERE/"drafts/material_propagation_decision.md"]
    return dict(round=960,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(U=U,hopping=v,kappa_for_958_coefficient=kappa,
            charge_weight=d,transition_gap=Delta,hbar=1.,
            r_definition="R*Delta/c",geometry="parallel dipoles perpendicular to separation"),
        actual_material_response=dict(spectral_rows=spectral_rows,
            maximum_spectral_error=spectral_error,triplet_dipole_error=dark_error,
            bright_transition_eigenstate_error=bright_error,
            chi_from_36_state_spectral_sum=chi_spectral,
            chi_from_polarizability_integral=chi_integral,
            chi_second_order_958=chi2,chi_exact_958=chi_exact,
            higher_electrostatic_order_difference=delta_chi,
            old_read_time_times_difference=old["exact_interface"]["controlled_phase_time"]*delta_chi),
        maxwell_kernel=dict(rows=rows,far_zone_moment=moment,
            far_coefficient_target=13.,far_coefficient_numeric=far_coefficient,
            x_cut_40_energy_tail_upper_for_r_at_least_1=far_quad_tail_F,
            x_cut_40_derivative_tail_upper_for_r_at_least_1=far_quad_tail_derivative),
        analytic_effective_domain=dict(r_upper=rmax,
            energy_relative_error_upper=B0*rmax*rmax,
            force_relative_error_upper=BF*rmax*rmax,
            tail_frequency_in_gap_units=Y,
            one_model_energy_tail_relative_upper=energy_tail,
            one_model_force_tail_relative_upper=force_tail,
            tail_hypothesis="|alpha_Q(i xi)|<=2*Delta*d/xi^2 above Lambda",
            two_completions_tail_difference_upper_is_sum=True,
            tails_do_not_bound_unknown_low_frequency_matching=True),
        finite_geometry_comparison=dict(R=R,a=a,coulomb_coefficient=C,
            exact_four_site_coupling=k_site,dipole_coupling=k_dip,
            dipole_to_site_leading_energy_ratio=geometry_energy_ratio,
            relative_energy_mismatch=geometry_energy_ratio-1,
            static_site_to_dipole_ratio_lower=1-3*z/4,
            static_site_to_dipole_ratio_actual=site_over_dip),
        scope=dict(same_material_leading_stationary_Maxwell_matching=True,
            record_phase_response_and_distance_force_share_coefficients=True,
            Maxwell_3plus1_and_dipole_coupling_are_physical_inputs=True,
            zero_temperature_fixed_spin_sector_stationary_domain=True,
            no_requirement_of_exact_UV_completion=True,
            physical_matching_above_material_band_not_certified=True,
            full_real_time_958_bound_transported=False,
            complete_QED_implementation_certified=False,
            full_parent_SM_Einstein_certified=False,
            six_protocols_jointly_certified=False,full_goal_completed=False,
            stop_dispersion_branch_after_this_interface=True),
        primary_sources=[
            "https://arxiv.org/abs/1111.4224",
            "https://arxiv.org/abs/1502.06129"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files})

def compare(old,new,path=""):
    if isinstance(new,dict):
        assert old.keys()==new.keys(),path
        for k in new:compare(old[k],new[k],path+"/"+k)
    elif isinstance(new,list):
        assert len(old)==len(new),path
        for i,(a,b) in enumerate(zip(old,new)):compare(a,b,path+f"/{i}")
    elif isinstance(new,float):
        assert math.isclose(old,new,rel_tol=1e-9,abs_tol=1e-12),(path,old,new)
    else:assert old==new,(path,old,new)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:
        with TARGET.open("x",encoding="utf-8") as fobj:
            json.dump(out,fobj,ensure_ascii=False,indent=2);fobj.write("\n")
    else:compare(read(TARGET),out)
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},
                     ensure_ascii=False,indent=2))
