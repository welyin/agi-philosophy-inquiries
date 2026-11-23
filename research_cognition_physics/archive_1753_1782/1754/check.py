"""Fourth-family spectator, single-Higgs charge discrimination and chiral determinant test."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())


def interval_rate(amplitude,domain):
    lo,hi=domain
    points=[lo,hi,min(hi,max(lo,-amplitude.real))]
    rates=[abs(1+x/amplitude)**2 for x in points]
    return [min(rates),max(rates)]


def forms():
    m=runpy.run_path(str(HERE/CFG["on_shell_loop_code"]))
    atlas=json.loads((HERE/CFG["loop_and_ATLAS"]).read_bytes())
    w,_,fermions=m["common_inputs"](m["read"]("inputs.json"))
    mh=atlas["higgs_GeV"]
    def F(mass): return m["loop_functions"](mh,mass)[1]
    ag=sum(F(mass) for mass,nc,q in fermions if nc==3)
    ap=m["loop_functions"](mh,w)[2]+sum(nc*q*q*F(mass) for mass,nc,q in fermions)
    return ag,ap,F,atlas


def exchange_masses():
    bg=runpy.run_path(str(HERE/CFG["physical_parent"]))["mother"](.125)
    old=json.loads((HERE/CFG["family_exchange"]).read_bytes())
    R=old["fourth_response_R"]
    results={}
    for key,mat in (("quark",bg["yW"][1,:,:3]),("lepton",bg["yW"][2,:,:3])):
        y=np.linalg.svd(mat,compute_uv=False)[::-1]
        square=y**2
        x=np.r_[(square-square[0])/(square[2]-square[0]),R]
        w=(1+x)**(-2)
        pairs=list(itertools.combinations(range(4),2))
        prob=np.array([w[i]*w[j] for i,j in pairs])
        prob/=prob.sum()
        J=np.array([sum(p*(v-x[i])*(v-x[j]) for p,(i,j) in zip(prob,pairs)) for v in x])
        np.testing.assert_array_equal(np.sign(J),[1,-1,1,1])
        frozen=old["two_sources"][key]["common_alignment_tests"]["2"]
        np.testing.assert_allclose(J,frozen["raw_selection"],rtol=3e-9)
        y4=np.sqrt(square[0]+R*(square[-1]-square[0]))
        results[key]=dict(R=R,y4_at_weak_heavy_threshold=y4,
                          assigned_pole_mass=CFG["v_GeV"]*y4/math.sqrt(2),
                          exchange_signs=np.sign(J),probability_pair_includes_fourth=sum(
                              p for p,(i,j) in zip(prob,pairs) if j==3))
    return results


def charged_data(ag,ap,atlas):
    old=runpy.run_path(str(HERE.parent/"1667"/"check.py"))
    prod=old["window"](atlas["VBF_over_gluon"],2)
    decay=old["window"](atlas["diphoton_over_ZZ"],2)
    br=CFG["bbH_SM_fraction_interval"]
    gg=[(1/r-b)/(1-b) for r,b in itertools.product(prod,br)]
    intervals=old["real_shift_intervals"](ag,1,[min(gg),max(gg)])
    bands=[[8/3+lo/3,4+hi/3] for lo,hi in intervals]
    gamma=[interval_rate(ap,band) for band in bands]
    assert max(z[1] for z in gamma)<decay[0]
    return dict(ggF_VBF_inverted_observation=prod,gamma_ZZ_observation=decay,
                extra_gluon_real_shift_intervals=intervals,
                down_only_delta_photon_intervals=bands,down_only_gamma_ZZ_intervals=gamma,
                minimum_enhancement_of_max_down_only_gamma_ZZ_to_reach_window=decay[0]/max(v[1] for v in gamma))


def generated_family(masses,ag,ap,F):
    MD=masses["quark"]["assigned_pole_mass"]
    ME=masses["lepton"]["assigned_pole_mass"]
    MU=CFG["fourth_up_pole_GeV"]
    M=CFG["vectorlike_common_mass_GeV"]
    def scalar_sum(q):
        return (F(MU)+F(MD)+q/(1+q)*F(M*(1+q))-q/(1-q)*F(M*(1-q)))
    a,b=CFG["q_interval"]
    assert scalar_sum(a)>0>scalar_sum(b)
    for _ in range(70):
        mid=(a+b)/2
        if scalar_sum(mid)>0:
            a=mid
        else:
            b=mid
    q=(a+b)/2
    assert abs(scalar_sum(q))<3e-12
    U0,D0,E0=F(MU),F(MD),F(ME)
    extra=-q/(1-q)*F(M*(1-q))+q/(1+q)*F(M*(1+q))
    dg=U0+D0+extra
    gammaU=(4*(U0+extra)+D0)/3+E0
    gammaD=(4*U0+D0+extra)/3+E0
    np.testing.assert_allclose([gammaU,gammaD],[E0-D0,E0+U0],atol=2e-12)
    # Fixed q and independent down/lepton pole matching envelopes.
    down=[F(1.2*MD),F(.8*MD)]
    lepton=[F(1.2*ME),F(.8*ME)]
    photon_shifts=[(4*(U0+extra)+d)/3+e for d,e in itertools.product(down,lepton)]
    rg=interval_rate(ag,[d-D0 for d in down])
    weighted=[1/((1-b)*r+b) for b,r in itertools.product(CFG["bbH_SM_fraction_interval"],rg)]
    rp=interval_rate(ap,[min(photon_shifts),max(photon_shifts)])
    assert min(weighted)>.70 and max(weighted)<1.88 and rp[0]>.63 and rp[1]<1.15
    # Soft derivative of log of the complete determinant, including the chiral pair.
    u=-q**2/(1-q**2)
    c1=2+2*u
    c2=c1-4-4*u*u
    assert abs(c2-(-c1*c1+5*c1-8))<1e-13
    unhidden_g=1/abs(1+(U0+D0)/ag)**2
    precision=CFG["precision_comparison"]
    def theta(a,b):
        if a==b:return np.zeros(2)
        x,y=a*a,b*b
        ell=np.log(x/y)
        return np.array([x+y-2*x*y*ell/(x-y),2*a*b*((x+y)*ell/(x-y)-2)])
    def weak(qU,qD):
        up=[M*(1+qU),M*(1-qU)]
        dn=[M*(1+qD),M*(1-qD)]
        charged=sum(sum(theta(i,j)) for i,j in itertools.product(up,dn))/2
        neutral=(sum(theta(*up))+sum(theta(*dn)))/2
        rho_vector=3*(charged-neutral)/(16*np.pi**2*CFG["v_GeV"]**2)
        rho_chiral=(3*theta(MU,MD)[0]+theta(ME,CFG["fourth_neutral_Dirac_pole_GeV"])[0])/(16*np.pi**2*CFG["v_GeV"]**2)
        return dict(vector_delta_rho=rho_vector,chiral_delta_rho=rho_chiral,
                    combined_T=precision["alpha_inverse"]*(rho_vector+rho_chiral))
    only=weak(q,0)
    np.testing.assert_allclose(weak(q,q)["vector_delta_rho"],0.,atol=1e-14)
    gtrain=precision["two_singlet_training_gamma_ratio"]
    ph=-np.sqrt(gtrain*abs(ap)**2-ap.imag**2)-ap.real
    def shift(z):return z/(1+z)*F(M*(1+z))-z/(1-z)*F(M*(1-z))
    def invert(target):
        lo,hi=0.,CFG["q_interval"][-1]
        assert shift(lo)>target>shift(hi)
        for _ in range(65):
            mid=(lo+hi)/2
            if shift(mid)>target:lo=mid
            else:hi=mid
        return (lo+hi)/2
    qu,qd=invert(ph-E0-U0),invert(-D0+E0-ph)
    both=weak(qu,qd)
    both.update(q_up=qu,q_down=qd,gluon_rate=abs(1+(U0+D0+shift(qu)+shift(qd))/ag)**2,
                photon_over_ZZ=abs(1+((4*(U0+shift(qu))+D0+shift(qd))/3+E0)/ap)**2)
    np.testing.assert_allclose([both["gluon_rate"],both["photon_over_ZZ"]],[1,gtrain])
    tc,tw=precision["Gfitter_free_STU_marginal_T"]
    assert abs(both["combined_T"]-tc)<tw and only["combined_T"]>tc+2*tw
    return dict(masses_chiral_UDEN=[MU,MD,ME,CFG["fourth_neutral_Dirac_pole_GeV"]],
                hard_vector_mass=M,q=q,direct_Yukawa_UDE=np.array([MU,MD,ME])*np.sqrt(2)/CFG["v_GeV"],
                vector_Yukawa=q*M*np.sqrt(2)/CFG["v_GeV"],
                vector_poles_including_inactive_weak_partner=M*np.array([1-q,1+q,1]),
                signed_kappas_including_chiral_UD=[1.,1.,-q/(1-q),q/(1+q),0.],
                chiral_only_VBF_over_ggF=unhidden_g,extra_gluon_amplitude=dg,
                VBF_over_grouped_ggF=1/abs(1+dg/ag)**2,
                photon_upcomp_downcomp_shifts=[gammaU,gammaD],photon_over_ZZ_upcomp_downcomp=[
                    abs(1+gammaU/ap)**2,abs(1+gammaD/ap)**2],
                finite_matching_envelope_fixed_q=dict(VBF_over_grouped=[min(weighted),max(weighted)],photon_over_ZZ=rp),
                one_up_singlet_precision_holdout=only,extra_down_singlet_precision_holdout=both,
                soft_g_coefficients_K1_K2=[c1,c2],zero_at_Higgs_background_t=1/q)


def polynomial_tests():
    # Taylor theorem in s=t^2: verify diverse chiral indices, as well as both explicit mass blocks.
    ntest=0
    for p,d in itertools.product(range(1,8),range(7)):
        a=Fraction(p,2)
        inv=[Fraction(1)]
        for i in range(1,d+2):
            inv.append(-inv[-1]*(a+i-1)/i)
        forward=[Fraction(1)]
        for i in range(1,d+2):
            forward.append(forward[-1]*(a-i+1)/i)
        terms=[sum(forward[k-i]*inv[i] for i in range(min(k,d)+1)) for k in range(d+2)]
        assert terms[0]==1 and all(q==0 for q in terms[1:d+1])
        predicted=(-1)**d*2**(d+1)*math.prod(a+j for j in range(d+1))
        assert 2**(d+1)*math.factorial(d+1)*terms[-1]==predicted
        ntest+=1
    C=np.array([[.5,-1/np.sqrt(12)],[1/np.sqrt(12),.5]])
    def matrix(t):
        return np.block([[np.eye(2),t*np.eye(2)],[t*C,np.eye(2)]])
    for t in (.1,.5,1.,1.3,2.1):
        det=t*t*np.linalg.det(matrix(t))/np.linalg.det(matrix(1.))
        np.testing.assert_allclose(det,1+(t*t-1)**3,atol=3e-13)
    m=matrix(1.)
    z=np.linalg.solve(m,matrix(1.)-matrix(0.))
    coeff=np.array([2+np.trace(z),-2-np.trace(z@z),4+2*np.trace(z@z@z)])
    np.testing.assert_allclose(coeff,[0,0,48],atol=2e-12)
    return dict(exact_rational_checks=ntest,minimal_real_linear_factor_upper_K2=-7/4,
                two_factor_real_gauge_blocks_K1_K2_K3=coeff,
                two_factor_spectrum_in_common_hard_mass_units=np.linalg.svd(m,compute_uv=False),
                two_factor_determinant_zeros_positive_real_away_from_origin=0)


def simple(v):
    if isinstance(v,dict):return {k:simple(x) for k,x in v.items()}
    if isinstance(v,(list,np.ndarray,tuple)):return [simple(x) for x in v]
    if isinstance(v,(float,np.floating)):return 0. if abs(v)<1e-12 else float(f"{v:.10g}")
    return v


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    masses=exchange_masses()
    ag,ap,F,atlas=forms()
    result=simple(dict(exchange_and_matched_masses=masses,
                       SM_one_loop_gluon_photon=[[ag.real,ag.imag],[ap.real,ap.imag]],
                       down_compensator_bound=charged_data(ag,ap,atlas),
                       family=generated_family(masses,ag,ap,F),determinant=polynomial_tests()))
    target=HERE/"results.json"
    if args.check:
        assert result==json.loads(target.read_bytes())
        print("1754 OK: four-family exchange, finite Higgs loops, charge intervals and soft determinant.")
    else:
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

