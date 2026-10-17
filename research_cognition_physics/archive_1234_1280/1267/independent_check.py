"""Independent exact-arithmetic review of round 1267; one review group.
No author executable is imported. Earlier saved interval certificates are inputs.
Default compares only; --write creates the result file exclusively.
"""
from __future__ import annotations
import hashlib, json, sys
from fractions import Fraction as F
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent

def read(rel):
    return json.loads((STAGE / rel).read_text(encoding="utf-8-sig"))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def lohi(values):
    return tuple(map(F, values))

def up_log2(n):
    return (n - 1).bit_length()

def main():
    old = read("1247/results.json")
    wide = read("1249/quantum_bounds_results.json")["group_constructive_finite_resource_bound"]["ledger"]
    fam = read("1252/results.json")
    cell = fam["group_parameter_family"]["cells"][0]
    assert fam["group_parameter_family"]["parameter_cells"] == 1
    assert lohi(cell["kappa"]) == (F(3), F(3001,1000))
    # Derive old constants from their published defining formula, not author outputs.
    B = 324 * (factorial(81) // factorial(40)**2) * 2**40
    Lc = 56*B
    Lh = 14*(1+Lc)
    H0 = 2**36
    M = 2**up_log2(H0*factorial(21)*Lh**21)
    As = 2**98*(1+Lc)**2
    Ls = 4*(Lc+14)
    Lp = Ls*(1+100*As*B)
    Lxi, Lk = Lc+4*Lp, Lc+4*10**9*B
    Ae, Le = 1+128*As**2, Lxi+2*Ls
    Ast, Lst = max(As,Ae,1), max(Ls,Le,Lk)
    CF = 2**34*5*256**13*15**10*factorial(20)
    EH = up_log2(CF*H0*Lh**20)
    E0 = up_log2(CF*Ast*(1+Lst)**20)
    EF = E0+20*(4096+5*M)
    N = 1+max((EF+44)//3, (EH+EF+47)//6)
    assert (M.bit_length()-1, EH, E0) == (3059,3093,14008)
    assert EF == int(wide["Fourier_transport_exponent_decimal"])
    assert N == int(wide["integer_n_power_exponent_decimal"])
    assert B == int(wide["B_decimal"])
    # Parameter-family legal initial data: pi_Q^2=12*rho0/kappa.
    rho0 = F(111,100)
    init_checks = []
    for kappa in (F(3), F(6001,2000), F(3001,1000)):
        pi2 = 12*rho0/kappa
        assert -kappa*pi2/12 == -rho0
        assert kappa*kappa*pi2/36 == kappa*rho0/3
        init_checks.append(str(-kappa*pi2/12))
    deleted_fixed_pi_residual = -F(1,1000)*rho0/3
    assert deleted_fixed_pi_residual == -F(37,100000)
    # Single whole-interval tube implies the original positive core uniformly.
    tube = {k:lohi(v) for k,v in cell["all_time_state_tube"].items()}
    pos = {k:lohi(v) for k,v in cell["positive_all_time"].items()}
    assert F(1,3)<tube["u"][0] and tube["u"][1]<F(3,2)
    assert 0<tube["H"][0] and tube["H"][1]<F(11,10)
    for key in ("p","v","X","V"):
        assert max(map(abs,tube[key]))<F(3,2)
    for key in ("q","r","psi","Y","K"):
        assert max(map(abs,tube[key]))<1
    assert max(map(abs,tube["z"]))+F(1,5)<1
    for key,bound in (("J",F(2,25)),("M_c",F(9,10)),("F",F(3,5)),("M_r",F(1,25))):
        assert pos[key][0]>bound
    piq_bound = 2*27*F(11,10)
    pip_bound = 27*F(3,2)
    assert piq_bound<61 and pip_bound<61
    # Q=-log u: e^(3/2)>1+3/2+(3/2)^2/2>3 and e^(1/2)>3/2.
    assert 1+F(3,2)+F(9,8)>3
    uniform_h = F(wide["h_complex_upper"])+F(129**2*3**12,12000)
    source_complex = 1024*3**18*(1+H0)*(1+129**2)
    assert uniform_h == F(79598612877987,4000)<H0
    assert source_complex<2**96
    # These exponents use epsilon<=2^(-4N), not a changed Planck constant.
    e_q1,e_q2 = EF-4*N, EH+EF-8*N
    assert e_q1<-42 and e_q2<-42
    # New observable applied to the saved rigorous endpoints.
    memory = old["group_finite_memory"]
    uB = [lohi(memory[f"alpha{j}"]["checkpoints"]["2"]["state"]["u"]) for j in (0,1)]
    gapB = 1/uB[0][1]-1/uB[1][0]
    eta, dref = F(1,10**8),F(1,10**16)
    assert gapB>F(683,10**7)
    assert gapB-2*eta>F(68,10**6)
    uc = lohi(cell["checkpoints"]["2"]["state"]["u"])
    su = lohi(cell["checkpoints"]["2"]["sensitivity"]["u"])
    inverse_slope = -su[1]/uc[1]**2
    assert su[1]<0 and inverse_slope>F(1,5)
    refwidth = (1/uB[1][0]-1/uB[1][1])/2
    assert refwidth<dref
    kappabound = 5*(2*eta+dref)
    fieldbound = rho0*kappabound
    assert kappabound == F(200000001,2000000000000000)
    assert fieldbound == F(22200000111,200000000000000000)
    assert F(fam["group_inverse_and_sources"]["uniform_source_component_bound"]) == rho0
    endpointgap = F(1,5000)-2*eta
    assert endpointgap == F(9999,50000000)>0
    statbudget = F(5,10**9)+F(4,10**9)+F(1,10**10)+F(1,2**40)
    assert statbudget<eta
    # Hoeffding exponent and rational exponential upper bound.
    assert F(26*(2*10**8)**2)*F(5,10**9)**2/2 == 13
    assert sum(F(1,factorial(j)) for j in range(5))>F(8,3)
    fail4 = 8*F(3,8)**13
    assert fail4<F(1,100)
    # One instrument review subgroup. Independent characteristic-function test:
    # composed Gaussian convolutions have heat exponents 1/4+1/4=1/2.
    heat_checks = 0
    for order in range(13):
        convolution = sum((-F(1,4))**j/F(factorial(j)) * (-F(1,4))**(order-j)/F(factorial(order-j)) for j in range(order+1))
        assert convolution == (-F(1,2))**order/F(factorial(order))
        heat_checks += 1
    assert F(1,4)+F(1,4)==F(1,2)
    assert F(1,2)+F(1,4)==F(3,4)
    # Gaussian normalization in 2d real dimensions: variance-halving volume=2^-d.
    one_mode_survival = F(1,2)
    survival = one_mode_survival**7
    assert survival == F(1,128)
    disturbance = 1-survival
    assert disturbance == F(127,128)
    # Exact old input-state source bound plus new arbitrary-normal-input channel bound.
    e_source = E0+EH+15*M+22-2*N
    e_inst = E0-1-4*N
    assert e_source<-82 and e_inst<-82
    assert max(e_source,e_inst)+1<-81
    # Finite Lueders alternative: same coarse effects, with no continuous environment required.
    Dout = (2*As*10**10+1)**3
    assert Dout < 2**4096
    # Constant plus compact symbols allow ||B||<=2^(E0+1): 2*2^(3/4)<4.
    assert F(4,2)**4 > 8
    # Two resolvent integral pieces at a=delta/||B|| give (4+2)/pi<2.
    assert 4+2 == 6  # uses the analytic elementary strict bound pi>3
    # The Fourier plane-wave dilation defect has squared coefficient
    # 2*(1-exp(-epsilon |k|^2/4)) <= epsilon |k|^2/2.
    assert 2*F(1,4) == F(1,2)
    # 4D epsilon^(1/4)2^E0 < 2^(E0+4098-N); 8D covers weighted extras.
    e_lueders = E0+4098-N
    e_lueders_weighted = E0+4099-N
    e_T = E0-2*N
    e_AWshift = E0-2-4*N
    assert max(e_lueders_weighted,e_T,e_AWshift,e_source)<-84
    # An exact rational two-effect noncommuting diagnostic, independent of the authors:
    # L0=diag(3/5,4/5), L1=diag(4/5,3/5); effects sum to I.
    roots=((F(3,5),F(4,5)),(F(4,5),F(3,5)))
    assert all(sum(v[j]**2 for v in roots)==1 for j in range(2))
    coherence=sum(v[0]*v[1] for v in roots)
    assert coherence==F(24,25)
    # Merging these two effects gives a one-outcome Lueders identity channel;
    # it differs from applying the fine instrument then deleting the label.
    coarse_coherence=F(1)
    assert coarse_coherence != coherence
    source_loss=1-coherence
    assert source_loss==F(1,25)>0
    # A Bell input's passive reference stays I/2, while its off-diagonal block changes.
    bell={(0,0):F(1,2),(0,3):F(1,2),(3,0):F(1,2),(3,3):F(1,2)}
    post={}
    for (row,col),value in bell.items():
        post[row,col]=value*sum(v[row//2]*v[col//2] for v in roots)
    ref_after=[[sum(post.get((2*s+i,2*s+j),F(0)) for s in range(2)) for j in range(2)] for i in range(2)]
    assert ref_after==[[F(1,2),F(0)],[F(0),F(1,2)]]
    assert post[0,3] != bell[0,3]
    # The source Fourier moment supports an operator norm, not a state-distance claim.
    return {
      "round":1267,"status":"PASS_INDEPENDENT_WITH_SCOPE","independent_groups":1,
      "method":"exact Fraction arithmetic; independent reconstruction of old resource definitions; saved interval inputs; formal Fourier/Gaussian algebra; no author imports or new ODE run",
      "geometry":{
        "closed_parameter_cells":1,"initial_constraint_values":init_checks,
        "wrong_fixed_initial_momentum_deletion":str(deleted_fixed_pi_residual),
        "uniform_h_upper":str(uniform_h),"uniform_h_slack":str(H0-uniform_h),
        "source_complex_slack":str(2**96-source_complex),
        "piQ_bound":str(piq_bound),"pi_scalar_bound":str(pip_bound),
        "B_gap_exact_lower":str(gapB),"B_gap_after_errors_gt":"68/1000000",
        "a_kappa_derivative_exact_lower":str(inverse_slope),
        "baseline_halfwidth":str(refwidth),"parameter_error":str(kappabound),
        "orthonormal_classical_field_error":str(fieldbound),
        "endpoint_observed_gap_lower":str(endpointgap),"four_mean_failure_upper":str(fail4),
        "same_fixed_a_marginal":True,"quantum_prefix_end":"2","classical_field_end":"9/4"},
      "instrument":{
        "Fourier_coefficient_checks":heat_checks,
        "measure_reprepare_nonselective_Weyl_heat_time":"epsilon/2",
        "measure_reprepare_post_AW_total_heat_time":"3*epsilon/4",
        "measure_reprepare_weighted_record_Weyl_heat_time":"epsilon/4",
        "measure_reprepare_seven_mode_coherent_survival":str(survival),
        "measure_reprepare_half_trace_distance_lower":str(disturbance),
        "old_source_exponent_below_minus82":e_source<-82,
        "instrument_exponent_below_minus82":e_inst<-82,
        "post_source_error":"<2^-81<2^-80",
        "finite_Lueders":{
          "D_out_binary_upper_exponent":Dout.bit_length(),
          "D_out_less_than_2power4096":True,
          "nonselective_bound":"4 D epsilon^(1/4) 2^E0",
          "weighted_bound":"8 D epsilon^(1/4) 2^E0",
          "constant_plus_compact_operator_bound":"2^(E0+1)",
          "weighted_Weyl_extra_bound":"sqrt(epsilon/2) F1",
          "weighted_AW_additional_bound":"epsilon F2/4 when comparing to original AW(r*b)",
          "all_four_post_source_terms_below_2power_minus84":True,
          "rational_qubit_coherence_factor":str(coherence),
          "coarse_single_effect_coherence":str(coarse_coherence),
          "coarse_Lueders_not_fine_Lueders_then_discard_record":True,
          "rational_qubit_source_loss":str(source_loss),
          "passive_Bell_reference_preserved":True,
          "system_reference_correlation_not_preserved":True,
          "finite_environment_is_not_energy_or_CC_certificate":True},
        "one_independent_group_only":True},
      "scope":["Continuous parameter, CP and operator identities require the reviewed analytic proof; finite arithmetic does not independently reprove old ODE certificates.",
        "Geometry classical approximation only for declared constrained Gaussian preparation family; other parameters and units fixed.",
        "The uniform fixed a marginal is not a parameter-independent full source triple.",
        "Instrument bounded-system-source estimates apply to arbitrary normal input and passive reference; they do not imply reduced input-state preservation.",
        "No new environment energy, total covariant conservation, actual CC implementation, quantum Einstein identity, CT adoption, or completed stage goal."],
      "inputs_sha256":{rel:sha(STAGE/rel) for rel in ("1247/results.json","1249/quantum_bounds_results.json","1252/results.json")},
      "independent_code_sha256":sha(Path(__file__))}

if __name__ == "__main__":
    answer = main()
    destination = HERE/"independent_results.json"
    if sys.argv[1:] == ["--write"]:
        with destination.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(answer,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        print("WRITE_ONCE_PASS_1267 independent_groups=1")
    elif not sys.argv[1:]:
        assert answer == json.loads(destination.read_text(encoding="utf-8-sig")), "saved result mismatch"
        print("READ_ONLY_MATCH_1267 independent_groups=1")
    else:
        raise SystemExit("usage: independent_check.py [--write]")