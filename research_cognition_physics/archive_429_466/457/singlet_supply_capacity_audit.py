"""Round 457: invariant supply capacity and a static exchange concentration witness.

Only Python/NumPy.  --dry-run recomputes without writing; default uses mode x.
"""
from __future__ import annotations
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ROUND = 457
BASELINE = 455
TARGET = Path(__file__).with_name('singlet_supply_capacity_audit_results.json')


def spin_one_multiplicities(n):
    rows = [{0: 1}]
    for _ in range(n):
        row = {}
        for j, count in rows[-1].items():
            for l in ([1] if j == 0 else (j-1, j, j+1)):
                row[l] = row.get(l, 0) + count
        rows.append(row)
    return rows


def halfspin_mult(n, j):
    # Number of copies of spin j among 2n spin-1/2 factors.
    if j < 0 or j > n:
        return 0
    return math.comb(2*n, n-j) - (math.comb(2*n, n-j-1) if n-j-1 >= 0 else 0)


def capacity(n, p):
    p = F(p)
    q = (1-p)/3
    rows = spin_one_multiplicities(n)
    answer = F(0)
    coefficients = [0]*(n+1)
    by_sector = []
    for j in range(n+1):
        remaining = halfspin_mult(n-1, j)
        cells = [(p**(n-k)*q**k, k, math.comb(n, k)*rows[k].get(j, 0))
                 for k in range(n+1)]
        used = []
        for value, k, count in sorted(cells, reverse=True):
            take = min(count, remaining)
            answer += (2*j+1)*take*value
            coefficients[k] += (2*j+1)*take
            remaining -= take
            if take:
                used.append([k, take])
        assert remaining == 0
        by_sector.append({'j': j, 'capacity': halfspin_mult(n-1, j), 'selected': used})
    return answer, coefficients, by_sector


def unrestricted(n, p):
    p = F(p); q = (1-p)/3; remaining = 4**(n-1); answer = F(0)
    for value, count in sorted([(p**(n-k)*q**k, math.comb(n,k)*3**k)
                                for k in range(n+1)], reverse=True):
        take = min(count, remaining); answer += take*value; remaining -= take
    return answer


def swap(n, a, b):
    index = np.arange(2**n)
    flip = ((index >> (n-1-a)) ^ (index >> (n-1-b))) & 1
    target = index ^ (flip << (n-1-a)) ^ (flip << (n-1-b))
    return np.eye(2**n, dtype=np.int64)[:, target]


def tensor(*items):
    out = np.array([[1]])
    for item in items:
        out = np.kron(out, item)
    return out


SINGLET2 = np.array([[0,0,0,0],[0,1,-1,0],[0,-1,1,0],[0,0,0,0]], dtype=np.int64)


def rho_pairs(n, p=.7):
    tau = (1-p)/3*np.eye(4) + (p-(1-p)/3)*SINGLET2/2
    return tensor(*([tau]*n))


def collective(n):
    d = 2**n
    plus = np.zeros((d,d))
    z = np.zeros(d)
    for k in range(d):
        z[k] = n/2-k.bit_count()
        for a in range(n):
            bit = 1 << a
            if k & bit:
                plus[k^bit, k] = 1
    return plus, np.diag(z)


def invariant_optimizer(n, p):
    """Construct a common-SU(2) unitary; NOT a static 2-body synthesis claim."""
    raw = 2*n; d = 2**raw
    jp, jz = collective(raw)
    rho = rho_pairs(n, p)
    effect = tensor(SINGLET2/2, np.eye(4**(n-1)))
    unitary = np.zeros((d,d), complex)
    for j in range(n+1):
        indices = np.flatnonzero(np.diag(jz) == j)
        _, singular, vh = np.linalg.svd(jp[:,indices], full_matrices=True)
        rank = int(np.sum(singular > 1e-9))
        highest = np.zeros((d, len(indices)-rank))
        highest[indices,:] = vh[rank:,:].T
        rv, rb = np.linalg.eigh(highest.T@rho@highest)
        ev, eb = np.linalg.eigh(highest.T@effect@highest)
        left = highest@eb[:,np.argsort(ev)[::-1]]
        right = highest@rb[:,np.argsort(rv)[::-1]]
        for m in range(j, -j-1, -1):
            unitary += left@right.conj().T
            if m > -j:
                scale = math.sqrt((j+m)*(j-m+1))
                left = jp.T@left/scale
                right = jp.T@right/scale
    return unitary, rho, effect, jp, jz


def rounded_divide(array, denominator):
    # Nearest integer, exact ties upward. Error of each real entry <= 1/2.
    return (array + denominator//2)//denominator


def dyadic_unitary_certificate(h5, time=93, squarings=12, degree=14, bits=100):
    """Integer fixed-point scaling/squaring with a rational norm error bound.

    H=h5/5. Horner evaluates exp(-itH/2^s), and every matrix entry is rounded
    with explicit <=1/(2B) real/imaginary error. No floating-point in certificate.
    """
    n = h5.shape[0]; scale = 1 << bits
    integer_h = h5.astype(object)
    real = np.eye(n, dtype=object)*scale
    imag = np.zeros((n,n), dtype=object)
    x = F(time*int(np.max(np.sum(np.abs(h5), axis=1))), 5*(1 << squarings))
    assert x < F(1,2)
    round_bound = F(n, scale)  # Frobenius bounds simultaneous complex rounding.
    error = F(0)
    for k in range(degree,0,-1):
        denominator = 5*(1 << squarings)*k
        next_real = rounded_divide(time*(integer_h@imag), denominator)
        next_imag = rounded_divide(-time*(integer_h@real), denominator)
        next_real += np.eye(n,dtype=object)*scale
        real, imag = next_real, next_imag
        error = x*error/k + round_bound
    tail = 2*x**(degree+1)/math.factorial(degree+1)
    error += tail
    for _ in range(squarings):
        nr = rounded_divide(real@real-imag@imag, scale)
        ni = rounded_divide(real@imag+imag@real, scale)
        real, imag = nr, ni
        # ||A^2-U^2|| <= (2+e)e when U is unitary.
        # Round UP to a fixed rational grid to prevent denominator explosion.
        bound = (2+error)*error+round_bound
        egrid = 10**35
        error = F((bound.numerator*egrid+bound.denominator-1)//bound.denominator, egrid)
    return real, imag, scale, error, x


def certificate_probability(real, imag, scale):
    # p=7/10: tau=(I+3*SINGLET2)/10, all integer numerator.
    r1 = np.eye(4,dtype=object)+3*SINGLET2.astype(object)
    rho_n = tensor(r1,r1,r1)
    # Tr[(P_s tensor I) U rho U*], P_s=SINGLET2/2.
    out_real = real@rho_n@real.T+imag@rho_n@imag.T
    effect2 = tensor(SINGLET2.astype(object),np.eye(16,dtype=object))
    numerator = sum((effect2*out_real.T).flat)
    return F(int(numerator), 2000*scale*scale)


def entropy(matrix):
    values = np.linalg.eigvalsh(matrix)
    values = values[values > 1e-13]
    return float(-np.sum(values*np.log2(values)))


def h4(p):
    q=(1-p)/3
    return sum(-x*math.log2(x) for x in [p,q,q,q] if x > 0)


def run_checks():
    checks = []
    # 1. Exact representation capacities, independent dimension identities.
    for n in range(1,9):
        rows=spin_one_multiplicities(n)
        for k,row in enumerate(rows):
            assert sum((2*j+1)*m for j,m in row.items()) == 3**k
        for j in range(n+1):
            assert sum(math.comb(n,k)*rows[k].get(j,0) for k in range(n+1)) == halfspin_mult(n,j)
        assert sum((2*j+1)*halfspin_mult(n-1,j) for j in range(n+1)) == 4**(n-1)
        for p in (F(0),F(1,10),F(1,4),F(7,10),F(1)):
            value,_,_=capacity(n,p)
            assert value <= unrestricted(n,p)
    checks.append({'name':'exact_sector_multiplicity_and_target_capacity','passed':True,'pairs_checked':8})

    # 2. First nontrivial non-Abelian charge penalty.
    p=F(7,10);q=(1-p)/3
    for n, expected in [(1,[1,0]),(2,[1,3,0]),(3,[1,9,6,0]),(4,[1,12,44,7,0])]:
        value,coeff,_=capacity(n,p);assert coeff == expected
        if n < 4:assert value == unrestricted(n,p)
    f4,_,sectors=capacity(4,p);loose=unrestricted(4,p)
    assert loose-f4 == 7*p*q*q*(p-q)>0
    zero_resource_bound,_,_=capacity(5,F(0))
    assert zero_resource_bound == F(184,243) and unrestricted(5,F(0)) == 1
    f3,_,_=capacity(3,p)
    assert f3-p == p*(1-p)*(4*p-1)/3
    checks.append({'name':'first_charge_penalty_and_three_pair_improvement','passed':True,
                   'p':str(p),'three_pair_optimum':str(f3),'four_pair_optimum':str(f4),
                   'four_pair_unrestricted':str(loose),'charge_penalty':str(loose-f4),
                   'five_triplet_only_pairs_invariant_optimum':str(zero_resource_bound),
                   'five_triplet_only_pairs_unrestricted_optimum':'1',
                   'four_pair_sectors':sectors})

    # 3. Full 256-dimensional invariant optimizer, not a 2-body compilation.
    u,rho,e,jp,jz=invariant_optimizer(4,float(p))
    errors={'unitarity':float(np.linalg.norm(u.conj().T@u-np.eye(256))),
            'Jplus_commutator':float(np.linalg.norm(u@jp-jp@u)),
            'Jz_commutator':float(np.linalg.norm(u@jz-jz@u))}
    attained=float(np.trace(e@u@rho@u.conj().T).real)
    assert max(errors.values()) < 2e-10 and abs(attained-float(f4))<1e-11
    checks.append({'name':'explicit_invariant_optimizer_attains_sector_bound','passed':True,
                   'attained_probability':attained,'matrix_errors':errors})

    # 4. A single fixed, positive, two-body exchange supply: rigorous integer certificate.
    h5=10*swap(6,0,1)+5*swap(6,2,3)+5*swap(6,4,5)
    cross=[(0,2),(0,3),(0,5),(1,2),(1,4),(1,5)]
    for a,b in cross:h5+=swap(6,a,b)
    r,i,scale,err,x=dyadic_unitary_certificate(h5)
    probability=certificate_probability(r,i,scale)
    prob_error=2*err+err*err
    lower=probability-prob_error
    assert err < F(1,10**19) and lower > F(809,1000)
    # Only 6 cross terms fail to commute with P01; ||[S,P]||<=1.
    # On |t-93|<=1/200 the probability drops by at most 6/(5*200)=.006.
    window_lower=lower-F(3,500)
    assert window_lower > F(4,5)
    hv,hb=np.linalg.eigh(h5/5)
    us=(hb*np.exp(-93j*hv))@hb.T
    rawrho=rho_pairs(3);output=us@rawrho@us.conj().T
    actual=float(np.trace(tensor(SINGLET2/2,np.eye(16))@output).real)
    assert abs(actual-float(probability))<1e-11
    checks.append({'name':'static_two_body_supply_integer_certificate_and_window','passed':True,
                   'edges_strong':[[0,1,2],[2,3,1],[4,5,1]],
                   'edges_one_fifth':[list(z) for z in cross], 'time':93,
                   'probability':actual,'certified_probability_greater_than':'809/1000',
                   'dyadic_bits':100,'taylor_degree':14,'squarings':12,'small_step_norm_bound':str(x),
                   'unitary_error_bound':str(err),'probability_error_bound':float(prob_error),
                   'window':['18599/200','18601/200'],'window_probability_greater_than':'4/5'})

    # 5. Resource ledger: total entropy kept, remainder necessarily changes.
    reshaped=output.reshape(4,16,4,16)
    reduced=np.trace(reshaped,axis1=1,axis2=3)
    rem=np.trace(reshaped,axis1=0,axis2=2)
    assert np.linalg.norm(reduced-rho_pairs(1,actual))<1e-11
    hs=entropy(reduced); hb_entropy=entropy(rem); total=entropy(output)
    mutual=hs+hb_entropy-total
    r1=np.eye(4,dtype=object)+3*SINGLET2.astype(object)
    rho_n=tensor(r1,r1,r1)
    energy=F(int(np.trace(h5.astype(object)@rho_n)),5000)
    energy2=F(int(np.trace(h5.astype(object)@h5.astype(object)@rho_n)),25000)
    assert abs(float(np.trace((h5/5)@output).real)-float(energy))<1e-10
    assert abs(total-3*h4(.7))<1e-10
    assert hb_entropy-2*h4(.7) >= h4(.7)-h4(actual)-1e-10 > 0
    assert mutual>0
    for n in range(1,9):
        f,_,_=capacity(n,p)
        assert f <= 1-(2*n+1)*q**n and f<1
    checks.append({'name':'closed_entropy_ledger_and_finite_charge_obstruction','passed':True,
                   'input_pair_entropy_bits':h4(.7),'target_entropy_bits':hs,
                   'remainder_entropy_bits':hb_entropy,'remainder_entropy_increase_bits':hb_entropy-2*h4(.7),
                   'minimum_required_increase_bits':h4(.7)-h4(actual),
                   'target_remainder_mutual_information_bits':mutual,
                   'constant_total_energy':str(energy),'constant_energy_second_moment':str(energy2),
                   'constant_energy_variance':str(energy2-energy*energy),
                   'total_entropy_error':abs(total-3*h4(.7))})

    # 6. Delivery contract: target mixed pair vs pure singlet, including internal remainder.
    s=np.array([0,1,-1,0])/np.sqrt(2)
    conditional=np.einsum('a,abcd,c->bd',s,reshaped,s)/actual
    ideal=tensor(SINGLET2/2,conditional)
    half_trace=.5*np.sum(np.abs(np.linalg.eigvalsh(output-ideal)))
    assert 1-actual-1e-10 <= half_trace <= math.sqrt(1-actual)+1e-10
    reduced_error=.5*np.sum(np.abs(np.linalg.eigvalsh(reduced-SINGLET2/2)))
    assert abs(reduced_error-(1-actual))<1e-11
    checks.append({'name':'unknown_subject_independent_delivery_and_retained_remainder','passed':True,
                   'reduced_supply_half_diamond_error':1-actual,
                   'closed_supply_half_trace_error':float(half_trace),
                   'closed_supply_upper_bound':math.sqrt(1-actual),
                   'old_unknown_subject_action':'identity during supply; whole GLR arbitrary',
                   'automatic_supply_to_conversion_handoff_implemented':False})
    return {'round':ROUND,'scientific_baseline_round':BASELINE,'test_count':len(checks),
            'checks':checks,'all_passed':all(c['passed'] for c in checks),
            'scope':{'exact_any_n_SU2_invariant_unitary_capacity':True,
                     'strict_nonabelian_charge_cost_at_four_pairs':True,
                     'one_static_positive_pair_exchange_supply_witness':True,
                     'unknown_subject_and_internal_remainder_preserved':True,
                     'static_exchange_reaches_global_optimum_claimed':False,
                     'perfect_singlet_from_finite_full_rank_supply':False,
                     'all_resources_regenerated':False,
                     'joint_autonomous_handoff_completed':False,
                     'new_cognitive_axiom_derived':False,
                     'full_GR_goal_completed':False,'phase_closure_triggered':False}}


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args();report=run_checks()
    if not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as file:
            json.dump(report,file,ensure_ascii=False,indent=2);file.write('\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))
