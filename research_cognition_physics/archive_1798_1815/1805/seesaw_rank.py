"""Light beta-endpoint and a finite number of canonical type-I ports at the common trace scale."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VEV = 246.22


def load(name):
    return json.loads((HERE / name).read_bytes())


def bounds():
    gauge = load("../1369/results.json")
    nu = load("../1370/inputs.json")
    mass2_spread = max(nu[o]["atm_abs"][2]+nu[o]["solar"][2] for o in ("NO","IO"))
    # m_beta^2 is a convex combination of the three squared masses.
    cap_eV = math.sqrt(nu["katrin"]["m_beta_upper_eV"]**2+mass2_spread)
    time = gauge["one_loop_corner_bounds"]["t"][1]
    boundary = gauge["one_loop_corner_bounds"]["mu_GeV"][1]
    # Declared sufficient upper envelope: lambda_(Hdag H)^2 <= .25.
    # Drop negative SU(2) and charged-lepton congruence contributions to upward growth.
    integral_quark = 10+0.8*(time-10)
    growth = math.exp((6*integral_quark+2*.001*time+time)/(16*math.pi**2))
    epsilon = 2*boundary*cap_eV*1e-9*growth/VEV**2
    return dict(max_mZ_singular_eV=cap_eV,mu_max_GeV=boundary,
                max_SM_operator_growth=growth,max_second_Y_squared=epsilon,
                max_second_Y=math.sqrt(epsilon))


def critical_lambda(central, eps):
    omega, norm = central["omega"], central["S"]
    charged_a, charged_b = central["A_ch"], central["B_ch"]
    neutral_a = (norm-charged_a)/omega
    assert neutral_a>3*eps
    neutral_b = (neutral_a-2*eps)**2+2*eps**2
    return {"neutral_A":neutral_a,"minimum_neutral_B_for_at_most_three":neutral_b,
            "neutrino_B_over_A_squared_min":neutral_b/neutral_a**2,
            "lambda_lower_at_most_three_pure_I":(charged_b+omega*neutral_b)/norm,
            "lambda_at_two_equal_large_sources":(charged_b+omega*neutral_a**2/2)/norm,
            "lambda_at_three_equal_large_sources":(charged_b+omega*neutral_a**2/3)/norm}


def examples(neutral_a, mass):
    # Any unitary left orientation works; rotate into a generic complex flavour frame.
    phases=np.diag(np.exp(1j*np.array([-.31,.75,.13])))
    omega=np.exp(2j*math.pi/3)
    v=phases @ (
        omega**np.outer(np.arange(3),np.arange(3))/math.sqrt(3))
    nu=load("../1370/inputs.json")["NO"]
    light_eV=np.sqrt([0.,nu["solar"][0],nu["atm_abs"][0]])
    d=2*mass*light_eV*1e-9/VEV**2
    target=v@np.diag(d)@v.T
    # For three singlets, a complex orthogonal boost fits a large second moment.
    double_t=math.acosh((neutral_a-d[0])/(d[1]+d[2]))
    c,h=math.cosh(double_t/2),math.sinh(double_t/2)
    orth=np.eye(3,dtype=complex)
    orth[1:,1:]=[[c,1j*h],[-1j*h,c]]
    cols={3:v@np.diag(np.sqrt(d))@orth}
    # Two or three independent cancelling pairs with arbitrary light mass difference.
    for n,ids in ((4,[1,2]),(6,[0,1,2])):
        Y=[]
        a=neutral_a/len(ids)
        for j in ids:
            Y.extend([v[:,j]*math.sqrt((a+d[j])/2),1j*v[:,j]*math.sqrt((a-d[j])/2)])
        cols[n]=np.column_stack(Y)
    results={}
    for n,Y in cols.items():
        eigen=np.linalg.svd(Y,compute_uv=False)[::-1]**2
        np.testing.assert_allclose(Y@Y.T,target,atol=5e-15,rtol=1e-10)
        np.testing.assert_allclose(eigen.sum(),neutral_a,rtol=5e-12)
        if n==3:
            assert eigen[-2] <= max(d)*(1+1e-9)
        else:
            np.testing.assert_allclose(eigen@eigen,neutral_a**2/(n/2),rtol=5e-12)
        results[str(n)]={"relative_C5_error":np.linalg.norm(Y@Y.T-target)/np.linalg.norm(target),
                         "A":eigen.sum(),"B_over_A_squared":(eigen@eigen)/neutral_a**2,
                         "Y_squared_singular_values":eigen.tolist()}
    return results


def calculate():
    bound=bounds()
    main=load("results.json")
    central=main["central"]
    spectral=critical_lambda(central,bound["max_second_Y_squared"])
    outer=main["outer"]
    wlo, whi=load("inputs.json")["gauge_ratio_outer"]
    neutral_min=outer["DeltaA_lower"]/whi
    eps=bound["max_second_Y_squared"]
    assert neutral_min>3*eps
    neutral_B_min=(neutral_min-2*eps)**2+2*eps**2
    all_lam=(outer["B_ch_lower"]+wlo*neutral_B_min)/outer["S"][1]
    matches=examples(spectral["neutral_A"],load("../1369/results.json")["one_loop_central"]["mu_GeV"])
    return {"transport_and_low_energy_bounds":bound, "conditional_central_gauge_and_charged_inputs":spectral,
            "box_lambda_min_in_three_singlet_type_I":all_lam,
            "box_quartic_matching_shift_at_least":all_lam-outer["lambda_data"][1],
            "identical_mass_coefficient_and_neutral_trace_with_three_four_six_singlets":matches}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    raw=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    path=HERE/"seesaw_rank_results.json"
    if parser.parse_args().check:
        assert json.loads(raw)==json.loads(path.read_bytes())
    else:
        path.write_text(raw,encoding="utf8")
    print("1805 seesaw rank",hashlib.sha256(raw.encode()).hexdigest()[:12])
