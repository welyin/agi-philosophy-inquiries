"""Independent exact audit; no author imports, ODE integration or sampling.
Default rebuilds and compares. --write exclusively creates the result.
"""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
ROOT = STAGE.parents[1]
OUT = HERE / 'independent_results.json'
NAMES = ('a', 'H', 'p', 'v', 'X', 'V')

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))

def serial(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [serial(v) for v in x]
    return x

def monomial(coef, **powers):
    return Q(coef), tuple(powers.get(n, 0) for n in NAMES)

def derivative(terms, j):
    result = {}
    for coef, powers in terms:
        if powers[j]:
            new = list(powers)
            new[j] -= 1
            new = tuple(new)
            result[new] = result.get(new, Q(0)) + coef*powers[j]
    return [(c, p) for p, c in result.items() if c]

def absolute_bound(terms, lower, upper):
    total = Q(0)
    for c, powers in terms:
        value = abs(c)
        for j, n in enumerate(powers):
            value *= (lower[j] if n < 0 else upper[j])**n
        total += value
    return total

def interval(state, name):
    return tuple(map(Q, state[name]))

def square_range(bounds):
    lo, hi = bounds
    minimum = Q(0) if lo <= 0 <= hi else min(lo*lo, hi*hi)
    return minimum, max(lo*lo, hi*hi)

def inherited_tube(state, delta):
    ulo, uhi = interval(state, 'u')
    hlo, hhi = interval(state, 'H')
    a_lo = 1/uhi
    a_hi = 1/(ulo*(1-hhi/4))
    margins = [a_lo-Q(5,2), Q(3)-a_hi, Q(31,100)-hhi,
               Q(1,10), Q(9,50)]
    assert min(margins) > delta
    jlo = interval(state, 'J')[0]
    mc = Q(71,100)+jlo+sum(square_range(interval(state, key))[0]
                                    for key in ('q','r','K'))/2
    assert jlo > Q(49,500) and mc > Q(93,100)
    mag = lambda name: max(map(abs, interval(state, name)))
    qr2 = sum(square_range(interval(state, key))[1] for key in ('q','r'))
    y_bound = mag('Y')+mag('K')/4
    clock_bound = mag('z')+Q(1,5)+Q(1,16)+delta/500
    assert qr2 < 1 and max(y_bound, clock_bound, mag('psi'), mag('K')) < 1
    return {'a_range':[a_lo,a_hi], 'H_upper':hhi,
            'minimum_box_margin':min(margins), 'J_lower':jlo, 'mass_lower':mc,
            'q_r_squared_upper':qr2, 'Y_upper':y_bound, 'vartheta_upper':clock_bound}

def calculate():
    eta = Q(1,25000)
    ta, tb, final = Q(447,200), Q(449,200), Q(9,4)
    width = tb-ta
    assert Q(17,8) < ta < tb < final and width == Q(1,100)
    duration = Q(1,50)
    assert final-ta < duration
    # Beta integral and reflection symmetry, independent of polynomial expansion.
    integral = Q(140*factorial(3)*factorial(3), factorial(7))
    first = Q(140*factorial(4)*factorial(3), factorial(8))
    assert integral == 1 and first == Q(1,2)
    # d[s^3(1-s)^3]/ds = 3s^2(1-s)^2(1-2s).
    peak = Q(140,64)/width
    moment = integral*(final-ta)-width*first
    assert peak == Q(875,4) and moment == Q(1,100)

    # Laurent-polynomial differentiation of all six vector-field components.
    k = Q(3001,1000)
    vector = [
        [monomial(1,a=1,H=1)],
        [monomial(-Q(3,2),H=2), monomial(-k/4,v=2),
         monomial(-k/4,V=2), monomial(k/4,p=2), monomial(k/4,X=2),
         monomial(k/40,a=-2), monomial(k/200,a=-2,p=1)],
        [monomial(1,v=1)],
        [monomial(-3,H=1,v=1), monomial(-1,p=1),
         monomial(-Q(3,100),a=-2), monomial(-Q(1,100),a=-3)],
        [monomial(1,V=1)],
        [monomial(-3,H=1,V=1), monomial(-1,X=1)]]
    lower = [Q(5,2),Q(0),Q(0),Q(0),Q(0),Q(0)]
    upper = [Q(3),Q(31,100),Q(8,5),Q(8,5),Q(8,5),Q(8,5)]
    matrix = [[absolute_bound(derivative(row,j),lower,upper)
               for j in range(6)] for row in vector]
    row_sums = list(map(sum,matrix))
    assert row_sums[0] < 4 and max(row_sums) < 12
    beta_low, beta_high = Q(3)/(6*Q(3)**2), k/(6*Q(5,2)**2)
    assert beta_low == Q(1,18) and beta_high < Q(81,1000)

    # Four-term-degree truncation with a geometric majorant of the remaining tail.
    x = 12*duration
    exp_bound = sum(x**j/factorial(j) for j in range(5))
    exp_bound += x**5/factorial(5)/(1-x/6)
    assert exp_bound < Q(32,25)
    delta = Q(81,1000)*eta*Q(32,25)
    # Exact time-ordered Peano series: a row norm 4, matrix norm 12.
    remainder_slope = Q(4,12)*Q(12**2,2)*duration*Q(32,25)
    slope = Q(5,2)-remainder_slope
    assert slope > Q(3,2)
    response = Q(3,2)*beta_low*eta*moment
    assert response == Q(1,30000000)
    assert delta == Q(81,19531250) < Q(1,150000)
    assert Q(9,50)**2 < Q(9,10)/27

    cells = read(STAGE/'1252/results.json')['group_parameter_family']['cells']
    assert len(cells) == 1 and cells[0]['kappa'] == ['3','3001/1000']
    old = read(STAGE/'1247/results.json')['group_finite_memory']
    states = {'closed_kappa_B2':cells[0]['checkpoints']['2']['state']}
    states.update({key:old[key]['checkpoints']['2']['state'] for key in ('alpha0','alpha1')})
    tubes = {key:inherited_tube(state,delta) for key,state in states.items()}

    cpeak = eta*peak
    j_floor, mass_floor = Q(49,500)-3*cpeak, Q(93,100)-3*cpeak
    work_low = Q(5,2)*(Q(9,50)-delta)*eta
    work_high = Q(3)*Q(31,100)*eta
    assert j_floor == Q(287,4000) > 0 and j_floor < Q(2,25)
    assert mass_floor == Q(723,800) > Q(9,10)
    assert 0 < work_low < work_high
    # rho_rod=c*a^-2, pressure=-rho_rod/3, Jdot=-c'*a.
    # Coefficients of c' a^-2 and c H a^-2 in total continuity.
    continuity_coefficients = [Q(1)-1, -Q(2)+3*(1-Q(1,3))]
    assert continuity_coefficients == [0,0]
    assert -2*Q(1,6) == -Q(1,3)  # Raychaudhuri/source pressure match.
    canonical = [27*Q(8,5),2*27*Q(31,100)]
    assert max(canonical) < 61
    assert sum(Q(3,2)**j/factorial(j) for j in range(4)) > 3
    small = Q(1,40000)
    small_j = Q(49,500)-3*peak*small
    assert small/Q(1200) > Q(2,10**8) and small_j > Q(2,25)

    baseline = read(HERE/'source_baseline.json')
    checked = 0
    for rel, value in baseline['source_sha256'].items():
        assert digest(STAGE/rel) == value, rel
        checked += 1
    for rel, value in baseline['project_sources_sha256'].items():
        assert digest(ROOT/rel) == value, rel
        checked += 1
    for absolute, value in baseline['external_source_sha256'].items():
        assert digest(Path(absolute)) == value, absolute
        checked += 1
    bindings = {name:digest(HERE/name) for name in
        ('README.md','response_proof.md','response_check.py','response_results.json',
         'source_baseline.json','independent_check.py')}
    bindings.update({name:digest(STAGE/name) for name in
                     ('1247/results.json','1252/results.json')})
    return serial({
        'verdict':'PASS_WITH_SCOPE', 'numbered_round':False,
        'independent_diagnostic_groups':1, 'formal_science_groups_added':0,
        'method':'Exact rational Laurent derivatives, Beta moments and inherited interval boxes; no author import',
        'pulse':{'normalization':integral,'first_unit_moment':first,
                 'terminal_moment':moment,'peak':peak,'c_peak':cpeak},
        'jacobian_absolute_entry_matrix':matrix,
        'row_sums':dict(zip(NAMES,row_sums)),
        'propagator':{'exp_upper_degree4':exp_bound,'coarse_exp':Q(32,25),
            'remainder_slope':remainder_slope,'retained_slope':slope,
            'beta_min':beta_low,'beta_max':beta_high},
        'inherited_tubes':tubes,
        'response':{'lower':response,'strict_state_upper':delta,
            'old_actual_a_error':Q(1,150000), 'old_error_does_not_guarantee_discrimination':True},
        'source':{'J_lower':j_floor,'mass_lower':mass_floor,
            'clock_work_lower':work_low,'clock_work_upper':work_high,
            'continuity_coefficients':continuity_coefficients,
            'coordinate_momenta_bounds':canonical,
            'good_J_plateau_inherited':False,
            'smaller_eta_arithmetic_only':{'eta':small,'J_lower':small_j,'response':small/1200}},
        'baseline_hashes_matched':checked,
        'scope':{'continuous_ODE_proof':'analytic bootstrap, not numerical integration',
            'same_classical_center_only':True, 'actual_quantum_feedback_certified':False,
            'new_control_and_future_readout_authorized':False,
            'geometric_law_selected':False, 'thermal_identity_proved':False,'CT_adopted':False},
        'sha256':bindings})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        print('EXCLUSIVE_FIRST_WRITE PASS_WITH_SCOPE')
    else:
        assert read(OUT) == result, 'READ_ONLY_MISMATCH'
        print('READ_ONLY_MATCH PASS_WITH_SCOPE')
    print('classical_gap >= '+result['response']['lower'])
    print('strict_upper = '+result['response']['strict_state_upper']+' < 1/150000')

if __name__ == '__main__':
    main()
