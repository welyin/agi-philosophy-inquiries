"""1210: role/permission zeroth law; exact finite arithmetic audit.
Default is read-only, compares the complete JSON with saved results.
--write creates or refreshes this round's results only. No old code imported.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import itertools
import json
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
SOURCES = ['1201/results.json', '1202/results.json', '1204/results.json',
           '1208/proof.md', '1209/proof.md']


def load(name):
    return json.loads((BASE / name).read_text(encoding='utf-8'))


def qstr(x):
    return str(x.numerator) + '/' + str(x.denominator)


def decimal_enclosure(lo, hi, digits=12):
    scale = 10 ** digits
    l = (lo.numerator * scale) // lo.denominator
    u = -((-hi.numerator * scale) // hi.denominator)
    return [qstr(F(l, scale)), qstr(F(u, scale))]


def log_bounds(x, terms=64):
    assert x > 0
    if x == 1:
        return F(0), F(0)
    if x < 1:
        lo, hi = log_bounds(1/x, terms)
        return -hi, -lo
    z = (x-1)/(x+1)
    power = z
    total = F(0)
    for k in range(terms):
        total += 2*power/(2*k+1)
        power *= z*z
    tail = 2*power/((2*terms+1)*(1-z*z))
    return total, total+tail


def exp_positive_bounds(x, degree=60):
    assert 0 <= x < degree+2
    term = F(1)
    total = term
    for k in range(1, degree+1):
        term *= x/k
        total += term
    next_term = term*x/(degree+1)
    tail = next_term/(1-x/(degree+2))
    return total, total+tail


def f(x):
    return 32000*x*x*x-31200*x*x-2095*x-17


def matvec(M, p):
    return [sum((x*y for x,y in zip(row,p)), F(0)) for row in M]


def td(p, q):
    return sum((abs(x-y) for x,y in zip(p,q)), F(0))/2


def permute_tuple(s, i, j):
    t=list(s)
    t[i],t[j]=t[j],t[i]
    return tuple(t)


def probability(s, states):
    x=F(1)
    for i,j in enumerate(s):
        x *= states[i][j]
    return x


def swap_audit(states, i, j):
    d=len(states[0])
    basis=list(itertools.product(range(d), repeat=len(states)))
    diffs=[]
    for s in basis:
        t=permute_tuple(s,i,j)
        # Identical-material swaps conserve the complete bare energy for
        # every possible physical spectrum, by preserving each level count.
        assert sorted(s)==sorted(t)
        assert permute_tuple(t,i,j)==s
        diffs.append(abs(probability(s,states)-probability(t,states)))
    return sum(diffs,F(0))/2, len(basis)


def run():
    j1,j2,j4=load('1201/results.json'),load('1202/results.json'),load('1204/results.json')
    M=[[F(x) for x in row] for row in j1['contacts']['600']['K']]
    inv=[[F(x) for x in row] for row in j1['inverse_rational']]
    for i in range(4):
        for j in range(4):
            assert sum(inv[i][k]*M[k][j] for k in range(4)) == int(i==j)
    linf=max(sum(abs(x) for x in row) for row in inv)
    em=F(j1['contacts']['600']['entry_error_upper'])
    eta=F(j1['finite_error_contract']['population_error_upper'])
    assert linf*(F(1,200)+em)<eta
    assert linf/(1-linf*4*em)<F(314,100)

    # Scalar-energy mediator and a non-scalar, partner-dependent two-port one.
    warm=[F(2,3),F(1,3)]      # beta = log 2, H = diag(0,1)
    cold=[F(4,5),F(1,5)]      # beta = log 4, same H and units
    scalar=[F(1,2),F(1,2)]    # H = 0: compatible with every beta
    assert warm[0]/warm[1]==2 and cold[0]/cold[1]==4
    # A-C energy-one degeneracy has unequal populations, excluding any Gibbs.
    cross_gap_difference=abs(warm[0]*cold[1]-warm[1]*cold[0])
    assert cross_gap_difference==F(2,15)
    d_ab, count=swap_audit([warm,warm,cold,cold],0,1)
    d_bc,_=swap_audit([warm,warm,cold,cold],2,3)
    d_ac,_=swap_audit([warm,warm,cold,cold],0,3)
    assert d_ab==d_bc==0 and d_ac==F(2,15)
    assert all(scalar[j]==scalar[0] for j in range(2))
    # Negative beta remains mathematically allowed; no positivity is inferred.
    inverted=list(reversed(warm))
    assert inverted[0]/inverted[1]==F(1,2)

    # Fixed-tolerance 'indistinguishability' is not an equivalence relation.
    approx=[[F(1,2),F(1,2)],[F(3,5),F(2,5)],[F(7,10),F(3,10)]]
    assert td(approx[0],approx[1])==td(approx[1],approx[2])==F(1,10)
    assert td(approx[0],approx[2])==F(1,5)

    # Native independently prepared pool marginals: not exact Gibbs states.
    rC=list(map(F,j2['one_marginal_rational']))
    rH=list(map(F,j4['new_reference']['one_marginal_rational']))
    assert sum(rC)==sum(rH)==1 and min(rC+rH)>0
    states=[rC,rC,rH,rH]
    n_ab,native_count=swap_audit(states,0,1)
    n_bc,_=swap_audit(states,2,3)
    n_ac,_=swap_audit(states,0,3)
    assert n_ab==n_bc==0 and n_ac>F(1916,10000)
    local=td(rC,rH)
    assert local>F(1898,10000)
    yC,yH=matvec(M,rC),matvec(M,rH)
    gaps=[abs(a-b) for a,b in zip(yC,yH)]
    event=max(range(4),key=lambda i:gaps[i])
    assert event==3
    model_lower=gaps[event]-2*em
    assert model_lower>F(11613,100000)
    # For TWO old thermometry groups, each total probability error <= .005.
    empirical_lower=model_lower-F(1,100)
    reported_CI_lower=model_lower-F(2,100)
    assert empirical_lower>F(10613,100000)
    assert reported_CI_lower>F(9613,100000)
    # These are old-readout predictions, not a new autonomous swap experiment.

    # Re-certify root brackets of h with rational signs and disjoint intervals.
    roots=[F(x) for x in j1['eigenvalues']]
    er=F(1,10**60)
    brackets=[]
    for i,x in enumerate(roots):
        if i==2:
            brackets.append((F(1),F(1)))
        else:
            lo,hi=x-er,x+er
            assert f(lo)*f(hi)<0
            brackets.append((lo,hi))
    assert all(brackets[i][1]<brackets[i+1][0] for i in range(3))
    wlo=brackets[3][0]-brackets[0][1]
    whi=brackets[3][1]-brackets[0][0]
    assert wlo>1

    # Gibbs controls are only calibrated controls from 1204, NOT pool states.
    calibration={}
    declared={'hot':(F(100,1000),F(578,1000)),
              'cold':(F(891,1000),F(1585,1000))}
    for name in ['hot','cold']:
        rec=j4['temperature_comparisons'][name]
        tau=list(map(F,rec['gibbs_population_rational']))
        et=F(rec['tau_entry_error_upper'])
        assert abs(sum(tau)-1)<=4*et
        assert min(tau)>2*eta+et
        # A real data estimate can already move by eta. The interval then
        # adds another eta; this 2 eta is essential for an ALL-valid-data claim.
        ratio_lo=(tau[0]-et-2*eta)/(tau[3]+et+2*eta)
        ratio_hi=(tau[0]+et+2*eta)/(tau[3]-et-2*eta)
        assert 1<ratio_lo<ratio_hi
        loglo,_=log_bounds(ratio_lo)
        _,loghi=log_bounds(ratio_hi)
        lower,upper=loglo/whi,loghi/wlo
        dl,du=declared[name]
        assert dl<lower<upper<du
        calibration[name]={
            'beta_control': rec['beta'],
            'gibbs_control_entry_error':qstr(et),
            'all_valid_data_beta_interval_outer': [qstr(dl),qstr(du)],
            'tighter_certified_outer':decimal_enclosure(lower,upper),
            'gibbs_control_population':rec['gibbs_population_rational'],
            'no_claim_pool_equals_control':True}
    sep=declared['cold'][0]-declared['hot'][1]
    assert sep==F(313,1000)

    # Non-scalar alone does not give uniform finite readout sensitivity.
    gap=F(1,1000)
    elo,ehi=exp_positive_bounds(gap)
    tiny_lo=1/(1+ehi)
    tiny_hi=1/(1+elo)
    tiny_dlo=F(1,2)-tiny_hi
    tiny_dhi=F(1,2)-tiny_lo
    assert 0<tiny_dlo<=tiny_dhi<F(1,4000)
    # A completely blind allowed readout cannot calibrate even a large gap.
    blind=[[F(1,2),F(1,2)],[F(1,2),F(1,2)]]
    assert matvec(blind,warm)==matvec(blind,cold)

    e8lo,_=exp_positive_bounds(F(8))
    failure=F(16)/e8lo
    assert failure<F(5368,10**6)
    return {
      'schema':'zeroth_law_role_permission_1210_v1',
      'status':'PASS_CONDITIONAL_ZEROTH_LAW_FINITE_CERTIFICATE',
      'source_sha256':{n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in SOURCES},
      'abstract_scope':{
        'joint_role':'1209 full finite energy-conserving internal three-budget comparison',
        'non_scalar_qualified_objects':'unique real beta; exact joint compatibility is an equivalence relation',
        'scalar_maximally_mixed_beta_set':'all real beta; cannot be a temperature-separating mediator',
        'no_positive_beta_selection':True,
        'common_physical_energy_units_required':True,
        'all_n_classification_is_analytic_not_finite_test':True},
      'finite_exact_permission_counterexample':{
        'qubit_warm_population':list(map(qstr,warm)),
        'qubit_cold_population':list(map(qstr,cold)),
        'middle_energy_levels':[0,1,1,2],
        'middle_state':'warm B1 tensor cold B2; distinct active port for each partner',
        'basis_states_checked':count,
        'AB1_swap_trace_distance':qstr(d_ab),
        'B2C_swap_trace_distance':qstr(d_bc),
        'AC_swap_trace_distance':qstr(d_ac),
        'energy_one_population_mismatch':qstr(cross_gap_difference),
        'endpoint_energy_preservation':'each basis level count preserved',
        'only_material_marginals_fixed_not_unknown_purifications':True},
      'native_permission_readout_bridge':{
        'cold_population_source':'1202/results.json:one_marginal_rational',
        'hot_population_source':'1204/results.json:new_reference.one_marginal_rational',
        'four_material_basis_states_checked':native_count,
        'AB1_swap_trace_distance':qstr(n_ab),
        'B2C_swap_trace_distance':qstr(n_bc),
        'AC_swap_trace_distance':qstr(n_ac),
        'A_population_trace_distance':qstr(local),
        'readout_event_zero_based':event,
        'cold_pointer_centers':list(map(qstr,yC)),
        'hot_pointer_centers':list(map(qstr,yH)),
        'pointer_center_sum_residuals':[qstr(sum(yC)-1),qstr(sum(yH)-1)],
        'single_event_true_gap_lower':'11613/100000',
        'all_budget_valid_empirical_gap_lower':'10613/100000',
        'all_budget_valid_reported_CI_gap_lower':'9613/100000',
        'extra_swap_implementation_error_not_silently_included':True,
        'not_claimed_exact_gibbs_or_full_role':True},
      'actual_readout_temperature_interface':{
        't':600,'g':'1/10',
        'entry_error':qstr(em),'inverse_center_inf_norm':qstr(linf),
        'actual_inverse_inf_upper':'157/50','population_error':qstr(eta),
        'confidence_interval_formula':'log((rhat_0-eta)/(rhat_3+eta))/w <= beta <= log((rhat_0+eta)/(rhat_3-eta))/w',
        'positive_endpoint_populations_required':True,
        'all_data_certificate_uses_double_eta':True,
        'controls':calibration,
        'separation_between_declared_outer_intervals':qstr(sep),
        'root_bracket_radius':qstr(er),
        'spectral_width_outer':decimal_enclosure(wlo,whi),
        'unknown_reference_thermal_role_not_certified_by_finite_fit':True},
      'finite_error_deletions':{
        'fixed_population_tolerance':'1/10',
        'AB_distance':'1/10','BC_distance':'1/10','AC_distance':'1/5',
        'only_additive_path_error_bound_not_fixed_tolerance_transitivity':True,
        'small_gap':'1/1000','beta_difference':'1',
        'small_gap_binary_probability_difference_upper':'1/4000',
        'small_gap_difference_outer':decimal_enclosure(tiny_dlo,tiny_dhi),
        'blind_readout_examples_equal':True},
      'statistics_per_two_reference_readout_comparison':{
        'fresh_iid_trials_each':250000,'groups':2,
        'statistical_error_per_probability':'1/250',
        'systematic_error_per_probability':'1/1000',
        'union_failure_formula':'16 exp(-8)',
        'union_failure_upper':'671/125000',
        'not_an_experiment':True,
        'four_groups_if_both_distinct_demonstrations_run_jointly_failure_upper':'671/62500'},
      'not_claimed':[
        'new adopted cognitive axiom or unconditional zeroth law',
        'original QQ controls every 1209 comparison',
        'actual finite pool marginals are exact Gibbs references',
        'all finite data certify all-n stability',
        'new autonomous swap controller or success-flag-free factory sampling',
        'unknown material-reference correlations return unchanged',
        'third law, geometry or all four laws completed']}


if __name__=='__main__':
    result=run()
    path=HERE/'results.json'
    text=json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+'\n'
    if '--write' in sys.argv:
        path.write_text(text,encoding='utf-8',newline='\n')
    else:
        assert path.exists(), 'run --write to save this new round first'
        assert json.loads(path.read_text(encoding='utf-8'))==result, 'saved result mismatch'
    print(text,end='')
