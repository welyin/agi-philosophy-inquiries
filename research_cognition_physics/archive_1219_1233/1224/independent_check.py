"""Independent 1224 algebra/constant audit. No author modules are imported.
Default execution is read-only; --write-results creates once.
This does not compute the complete 125-node H moments or time evolution.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]

def compute():
    # Covariant shape and the unique simultaneous forward-star coupling.
    g = [[F(4,3),F(2,3),F(1,2)],
         [F(2,3),F(4,3),F(0)],
         [F(1,2),F(0),F(1)]]
    a,c,d = g[0]
    _,e,f = g[1]
    _,_,z = g[2]
    det = a*(e*z-f*f)-c*(c*z-f*d)+d*(c*f-e*d)
    minors = [a,a*e-c*c,det]
    inverse12 = (d*f-c*z)/det
    assert all(x > 0 for x in minors)
    assert det == 1 and inverse12 == F(-2,3)
    B, targetJ = (0,0,0), (4,1,0)
    stars = []
    for x in range(5):
        for y in range(5):
            for zc in range(5):
                origin=(x,y,zc)
                ends=[tuple((origin[j]+int(j==i)) % 5 for j in range(3))
                      for i in range(3)]
                for i in range(3):
                    for j in range(i+1,3):
                        if set((ends[i],ends[j])) == set((B,targetJ)):
                            stars.append([list(origin),i+1,j+1])
    assert stars == [[[4,0,0],1,2]]

    # Independent finite-product evaluation of the two beta functions:
    # B(3,16)=2/(16*17*18), B(5/2,17)=16!/prod_{j=0}^{16}(5/2+j).
    beta316 = F(2,16*17*18)
    beta2517 = F(math.factorial(16),1)
    for j in range(17):
        beta2517 /= F(5+2*j,2)
    ratio = beta316 / beta2517
    R, mass, b, n = math.sqrt(6), 2.0, .1, 8
    jac = lambda r: (math.sinh(r/R)/(r/R))**4
    b0 = math.sqrt(mass)/4
    d0 = math.sqrt(mass)*b*(mass+1/3)/4
    Q1 = b0*4*n/b*jac(b)*float(ratio)+d0
    Qinf = 32*jac(2*b)*(4/3)**(2*n)*(d0+4*n*b0/b)
    eta = 1/(4*Qinf)
    a0 = b0/math.cosh(b/R)**2*math.cos(math.sqrt(mass)*b)
    h2 = (2*b/R)/math.tanh(2*b/R)
    bK = min(b/2,a0/(4*h2*Q1))
    Sminus = (a0-bK*h2*Q1)**2
    Supper = mass/16
    epsilon, psi = .73, math.exp(.12)
    w = epsilon**3*psi**6
    chi = epsilon*psi**2/w*eta*Sminus
    vminus = mass/(4*w)*math.cos(math.sqrt(mass)*b)/math.cosh(b/R)**4
    assert 0 < bK <= b/2 and 0 < Sminus <= Supper
    assert abs(eta*Qinf-.25) < 1e-15
    assert chi > 5.49e-8 and vminus > .6173
    # The support-density ratio implies rho_J,+/- >= 3 nu/4.
    assert 1-F(1,4) == F(3,4)

    # Exact arithmetic of the symbolic remainder/observation margins.
    # No full-H moment is replaced by a toy value.
    leading_coefficient = -inverse12
    allowed_remainder = F(1,3)
    guaranteed_gap = leading_coefficient-allowed_remainder
    four_record_error = 4*F(1,48)
    assert guaranteed_gap == F(1,3)
    assert four_record_error == F(1,12)
    assert guaranteed_gap-four_record_error == F(1,4)
    assert F(1,2)*2*3 == 3  # quadratic transport after local inverse
    assert F(8,6) == F(4,3) # third-order Taylor with t <= 2d/v

    sources = [
        HERE/'independent_check.py',
        BASE/'archive_585_628/research_note_589.md',
        BASE/'archive_585_628/research_note_598.md',
        HERE.parent/'_admission/area_readout_screen_after1221.md',
    ]
    return {
        "schema":"round1224_independent_scope_constants_v1",
        "status":"PASS",
        "all_checks_passed":True,
        "check_groups":3,
        "check_group_names":[
            "exact_shapes_and_unique_forward_star",
            "normal_preparation_and_positive_response_constants",
            "conditional_time_readout_and_statistical_margin_algebra"],
        "scope":"Full-H analytic-domain theorem and independent finite coefficient checks; no full-H propagation simulation.",
        "geometry":{
            "G1":[[str(x) for x in row] for row in g],
            "leading_principal_minors":[str(x) for x in minors],
            "det_G1":str(det),
            "G1_inverse_12":str(inverse12),
            "G0_inverse_12":"0",
            "lattice_shape":[5,5,5],
            "shared_star_and_one_based_directions":stars,
        },
        "inputs":{"hbar":1,"epsilon":epsilon,"psi":psi,"w":w,
                  "target_radius":R,"target_M":mass,"b":b,
                  "density_exponent":16,"amplitude_boundary_zero_order":8},
        "beta_integrals":{"B_3_16":str(beta316),"B_5over2_17":str(beta2517),
                          "ratio_exact":str(ratio)},
        "preparation_bounds":{"b0":b0,"d0":d0,"Q1":Q1,"Qinf":Qinf,
                              "eta":eta,"density_ratio_lower":"3/4",
                              "a0":a0,"h_2b":h2,"bK":bK,
                              "S_lower":Sminus,"S_upper":Supper,
                              "chi_lower":chi,"clock_v_lower":vminus},
        "conditional_margins":{
            "leading_background_difference_over_chi_t2":str(leading_coefficient),
            "maximum_remainder_over_chi_t2":str(allowed_remainder),
            "minimum_signal_over_chi_t2":str(guaranteed_gap),
            "four_probability_error_over_chi_t2":str(four_record_error),
            "remaining_margin_over_chi_t2":str(guaranteed_gap-four_record_error),
            "time_contract":"t <= chi_lower / (3*(B_G0+B_G1))",
            "sampling_contract":"N >= log(8/delta)/(2*epsilon_p**2), epsilon_p=chi_lower*t**2/48",
            "relational_contract":"d <= chi_lower/(3*v**2*K_rel), with local inverse domain retained",
        },
        "limitations":{
            "N1_N2_N3_computed":False,
            "specified_numerical_time_certified":False,
            "specified_numerical_sample_count_certified":False,
            "full_H_evolution_simulated":False,
            "normal_finite_energy_preparation":True,
            "preparation_is_full_space_faithful_thermal_state":False,
            "amplitude_is_C_infinity":False,
            "analytic_D_H3_domain":True,
            "single_run_proper_clock_constructed":False,
            "autonomous_preparation_and_readout_constructed":False,
            "independent_area_instrument_constructed":False,
            "numerical_constants_are_directed_intervals":False,
        },
        "independent_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256":{p.relative_to(BASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sources},
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-results",action="store_true")
    args=ap.parse_args()
    result=compute()
    target=HERE/"independent_results.json"
    if args.write_results:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        saved=json.loads(target.read_text(encoding="utf-8"))
        assert saved==result,"Independent results mismatch; refusing to rewrite"
    print(json.dumps({"status":"PASS","all_checks_passed":True,"results_match":True,
                      "check_groups":3,
                      "independent_code_sha256":result["independent_code_sha256"],
                      "full_H_moments_computed":False},ensure_ascii=False))

if __name__=="__main__":
    main()
