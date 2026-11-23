"""Heavy-light portals on the two lepton doublets and the coloured singlet."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
L=16*math.pi**2


def vectors():
    """Compute hard Coleman-Weinberg coefficients from separated singular values."""
    rng=np.random.default_rng(1808)
    def complex_array(shape):
        return (rng.normal(size=shape)+1j*rng.normal(size=shape))/5
    data={}
    for n in [1,2]:
        y=complex_array((3,3))
        p=complex_array((3,n))
        a=float(np.trace(p.conj().T@p).real)
        b=float(np.linalg.norm(y.conj().T@p)**2)
        def hard_trace(x):
            mat=np.block([[math.sqrt(x)*y,math.sqrt(x)*p],
                          [np.zeros((n,3)),np.eye(n)]])
            vals=np.linalg.svd(mat,compute_uv=False)[:n]**2
            return sum(vals**2*(np.log(vals)-1.5))+1.5*n
        def extract(h):
            return (hard_trace(h)+2*a*h)/h**2
        h=.003
        approx=2*extract(h/2)-extract(h)
        assert abs(approx+2*b)<1e-5
        data[str(n)]=dict(a=a,b=b,twice_negative_b=-2*b,heavy_loop_trace_coefficient=approx)
    # MSbar two-point bubble at the common natural pole; unit mass and coupling.
    nodes,weights=np.polynomial.legendre.leggauss(16)
    def kinetic(qsq):
        reduced=1-np.dot(weights/2,np.log1p(qsq*(1+nodes)/2))
        return 1+(1+qsq)*reduced
    h=.0004
    derivative=(kinetic(h)-kinetic(-h))/(2*h)
    np.testing.assert_allclose(derivative,.5,atol=1e-7)
    data["wavefunction_finite_bubble"]=float(derivative)
    return data


def comparison(p):
    masses,T=p["physical_D_T_GeV"]
    assert len(masses)==2 and masses[0]==masses[1] and len(T)==1
    M_D,M_T=masses[0],T[0]
    actual=np.array(p["squared_yukawas_at_Q"])
    S,A_ch,B_ch=p["S"],p["charged_A"],p["charged_B"]
    X=p["missing_A"]
    np.testing.assert_allclose(X,S-A_ch,rtol=1e-8)
    cd,ce=actual[3],np.mean(actual[6:8])
    other=.2+12*(ce-cd)/(5*X)
    cases={}
    for name,x in [("only_T",1.),("also_D",other)]:
        # Three additional *colour* copies of p_T but one Dirac quark.
        aT=x*X/3
        aD=(1-x)*X/2  # for one of the two different doublets
        yq=np.column_stack([np.diag(np.sqrt(actual[3:6])),[math.sqrt(aT),0,0]])
        ye=np.vstack([np.diag(np.sqrt(actual[6:9])),np.diag(np.sqrt([aD,aD,0]))[:2]])
        gq,ge=yq@yq.T,ye.T@ye
        added_b=3*(np.trace(gq@gq)-sum(actual[3:6]**2))+(np.trace(ge@ge)-sum(actual[6:]**2))
        form=X**2*(x**2/3+(1-x)**2/2)+2*X*(x*cd+(1-x)*ce)
        np.testing.assert_allclose(added_b,form,rtol=1e-12)
        lam=(B_ch+added_b)/S
        local_mass_cost=2*(3*aT*M_T**2+2*aD*M_D**2)/L
        mass_budget=1e6
        cases[name]=dict(colour_fraction=x,one_T_row_norm_squared=aT,
            each_D_row_norm_squared=aD,added_B=float(added_b),same_boundary_lambda_H=float(lam),
            conditional_unit_transport=dict(T_SM_H_partial_width_over_mass=aT/(16*math.pi),
                each_charged_D_SM_H_partial_width_over_mass=aD/(32*math.pi),
                Higgs_wavefunction_sum=X/(2*L),
                additive_natural_matching_mass_cost_GeV2=local_mass_cost,
                limit_on_common_local_over_Q_Yukawa_magnitude_if_mass_cost_below_1TeV_squared=
                    math.sqrt(mass_budget/local_mass_cost)))
    np.testing.assert_allclose(cases["only_T"]["added_B"],cases["also_D"]["added_B"],rtol=1e-12)
    first=cases["only_T"]["conditional_unit_transport"]["additive_natural_matching_mass_cost_GeV2"]
    second=cases["also_D"]["conditional_unit_transport"]["additive_natural_matching_mass_cost_GeV2"]
    return dict(Q=p["Q_GeV"],X=X,M_T=M_T,M_D=M_D,observed_passive_lambda=p["lambda_observed"],
        sum_width_relation_at_threshold="X_local=48pi Gamma_TH/M_T+32pi sum_i Gamma_DiH/M_D",
        min_largest_partial_width_over_mass_if_each_Q_to_pole_squared_gain_at_most_one=X/(112*math.pi),
        equal_moments_different_pole_signals=cases,
        ratio_of_matching_mass_costs_conditional_unit_transport=second/first)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    old=json.loads((HERE/"results.json").read_bytes())
    out={key:comparison(old[key]["two_doublets_one_triplet"]["point"])
         for key in ["common_MZ","independent_lattice"]}
    out["direct_complex_SVD_and_wavefunction"]=vectors()
    raw=json.dumps(out,indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==out
    else:
        dest.write_bytes(raw.encode())
    print("1808 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])

