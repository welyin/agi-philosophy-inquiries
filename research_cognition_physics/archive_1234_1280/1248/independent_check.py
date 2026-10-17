"""Independent 1248 certificate: exact hat weights and full-source bounds.

No author checker or interval classes are imported and no old evolution is
rerun.  The 1247 endpoint enclosures are read as inherited certificates.
--write creates a result exclusively; default execution only compares it.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / '1247'
RESULT = HERE / 'independent_results.json'
N, L = 128, Q(1, 4)
DS, DT = Q(1, N), L / N
ETA, ETA_K = Q(1, 10**8), Q(1, 10**10)


def packed(x):
    if isinstance(x, Q):
        return {'exact': str(x), 'decimal': float(x)}
    if isinstance(x, dict):
        return {k: packed(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [packed(v) for v in x]
    return x


def integral(poly, lo, hi):
    return sum((c * (hi**(i+1) - lo**(i+1)) / (i+1)
                for i, c in enumerate(poly)), Q(0))


def times_linear(poly, a, b):
    out = [Q(0)] * (len(poly) + 1)
    for i, c in enumerate(poly):
        out[i] += a*c
        out[i+1] += b*c
    return out


def w(s):
    return 30*s*s*(1-s)**2


def cdf(s):
    return 10*s**3 - 15*s**4 + 6*s**5


def first_moment(s):
    return Q(15, 2)*s**4 - 12*s**5 + 5*s**6


def weights():
    """Two independently expressed exact integration formulae agree."""
    wa, wb = [], []
    poly = [Q(0), Q(0), Q(30), Q(-60), Q(30)]
    p2 = [Q(60), Q(-360), Q(360)]
    for j in range(N+1):
        s = Q(j, N)
        # Integrate the derivative weight by parts against the linear hat.
        # Boundary terms use w(0)=w(1)=w'(0)=w'(1)=0.
        if j == 0:
            a = w(DS)/DS/L**2
        elif j == N:
            a = w(1-DS)/DS/L**2
        else:
            a = (w(s-DS) - 2*w(s) + w(s+DS))/DS/L**2
        b = Q(0)
        aa, bb = Q(0), Q(0)
        if j:
            lo = s-DS
            b += ((first_moment(s)-first_moment(lo))
                  - lo*(cdf(s)-cdf(lo)))/DS
            aa += integral(times_linear(p2, -lo/DS, 1/DS), lo, s)/L**2
            bb += integral(times_linear(poly, -lo/DS, 1/DS), lo, s)
        if j < N:
            hi = s+DS
            b += (hi*(cdf(hi)-cdf(s))
                  - (first_moment(hi)-first_moment(s)))/DS
            aa += integral(times_linear(p2, hi/DS, -1/DS), s, hi)/L**2
            bb += integral(times_linear(poly, hi/DS, -1/DS), s, hi)
        assert a == aa and b == bb and b > 0
        wa.append(a)
        wb.append(b)
    assert sum(wa) == 0 and sum(wb) == 1
    assert sum(wa[j]*Q(j, N) for j in range(N+1)) == 0
    assert sum(wb[j]*Q(j, N) for j in range(N+1)) == Q(1, 2)
    assert sum(abs(a) for a in wa) < 384
    # ||phi''||_1 = (40 sqrt(3)/3)/L^2 < 384, since sqrt(3)<9/5.
    assert Q(3) < Q(9, 5)**2
    assert Q(40, 3)/L**2 * Q(9, 5) == 384
    # A quadratic realizes the standard interval interpolation envelope h^2/8
    # when its second derivative is 1 (sign independent).
    assert DT**2/4 - DT**2/8 == DT**2/8
    eta_b = (3*4*ETA + Q(3, 10)*ETA + 4*ETA**2)/2
    e_full = DT**2*(384+20)/8 + 384*ETA + eta_b
    e_minus = DT**2*(384+21)/8 + 384*ETA + eta_b + Q(1, 10**10)
    assert e_full > 384*ETA and e_minus > e_full
    return {
        'status': 'PASS',
        'nodes': N+1, 'mesh': DT,
        'normalization': sum(wb),
        'derivative_hat_l1': sum(abs(a) for a in wa),
        'continuous_derivative_l1_exact': '640*sqrt(3)/3',
        'continuous_derivative_l1_strict_upper': Q(384),
        'b_noise': eta_b,
        'E_full': e_full, 'E_omitted': e_minus,
        'fixed_noise_floor': 384*ETA + eta_b,
        'weights_sha256': hashlib.sha256(json.dumps(
            [[str(x) for x in wa], [str(x) for x in wb]],
            separators=(',', ':')).encode()).hexdigest(),
        'scope': 'Exact finite-element integration weights and conditional error bounds; not sampled observations.'
    }


def domain_and_source():
    old = json.loads((OLD/'results.json').read_text(encoding='utf-8'))
    checkpoints = []
    for key in ('alpha0', 'alpha1'):
        cp = old['group_finite_memory'][key]['checkpoints']['2']
        s, e = cp['state'], cp['energy']
        lo = lambda a: Q(a[0])
        hi = lambda a: Q(a[1])
        assert lo(s['u']) > Q(187, 500)
        assert hi(s['u']) < Q(3, 8)
        assert 0 < lo(s['H']) <= hi(s['H']) < Q(27, 100)
        assert 0 < lo(e['rho']) <= hi(e['rho']) < Q(3, 40)
        assert 0 < lo(e['M_c']) <= hi(e['M_c']) < Q(117, 125)
        assert max(abs(lo(s['K'])), abs(hi(s['K']))) < Q(3, 1000)
        checkpoints.append({'alpha': key, 'rho_upper': hi(e['rho']),
                            'u_lower': lo(s['u']), 'K_absolute_upper': max(abs(lo(s['K'])), abs(hi(s['K'])))})

    # The inherited positivity/bootstrap gives positive field, rod, dust
    # contributions and |P|<=rho.  Monotone H and rho retain the endpoint bounds.
    rho, H, u, z = Q(3,40), Q(27,100), Q(3,8), Q(2,5)
    assert z*z > 2*rho
    u_lower = Q(187,500)*(1-H*L)
    assert u_lower > Q(1,3)
    a, ap = Q(3), Q(81,100)
    Hp = 3*rho
    up = H*u
    upp = (H*H+Hp)*u
    assert Hp == Q(9,40) and up < Q(11,100) and upp < Q(3,25)
    # Use the same public rational relaxations, but rederive each by the
    # complete independent chain rule (rods have F=1+p/5).
    up, upp = Q(11,100), Q(3,25)
    vp = 3*H*z+z+Q(3,100)*u*u+Q(1,100)*u**3
    Vp = 3*H*z+z
    assert vp < Q(3,4) and Vp < Q(3,4)
    vpp = 3*Hp*z+3*H*Q(3,4)+z+Q(3,50)*u*up+Q(3,100)*u*u*up
    Vpp = 3*Hp*z+3*H*Q(3,4)+z
    assert vpp < Q(13,10) and Vpp < Q(13,10)
    F = 1+z/5
    rodp = Q(3,100)*z*u*u+Q(3,10)*F*u*up
    rodpp = Q(3,100)*Q(3,4)*u*u+Q(3,25)*z*u*up+Q(3,10)*F*(up*up+u*upp)
    assert rodp < Q(2,125) and rodpp < Q(3,100)
    Pp = 2*z*(Q(3,4)+z)+Q(2,125)/3
    Ppp = 2*(Q(3,4)**2+z*Q(13,10)+z*z+z*Q(3,4))+Q(3,100)/3
    assert Pp < 1 and Ppp < Q(16,5)
    rhop = 3*H*2*rho
    assert rhop < Q(1,8)
    rhopp = 3*Hp*2*rho+3*H*(Q(1,8)+1)
    assert rhopp < Q(11,10)
    R, Rp, Rpp = 4*rho, Q(1,8)+3, Q(11,10)+3*Q(16,5)
    assert R == Q(3,10) and Rp == Q(25,8) and Rpp < 11
    app = a*(Hp+H*H)
    assert app < 1
    bpp = (R+2*ap*Rp+3*11)/2
    assert bpp < 20

    # Retained K is queried only at the declared terminal time.  All 129
    # samples and their source/clock/measurement deviations belong to the
    # *new*, unimplemented total menu contract, not to the old K permission.
    mk = Q(11,1000)
    lost_bpp = mk*(up*up+u*upp)
    assert lost_bpp < Q(1,1000)
    assert bpp+lost_bpp < 21
    k_abs = Q(3,1000)
    delta_m = k_abs*ETA_K + ETA_K**2/2
    a_observed_min = Q(8,3)-ETA
    subtraction_noise = delta_m/(2*a_observed_min**2) + mk*ETA/a_observed_min**3
    assert a_observed_min > 0 and subtraction_noise < Q(1,10**10)
    eta_b = (3*4*ETA+Q(3,10)*ETA+4*ETA**2)/2
    e = DT**2*(384+20)/8+384*ETA+eta_b
    em = DT**2*(384+21)/8+384*ETA+eta_b+Q(1,10**10)
    gap = Q(1,100)/2*Q(1,3)**2
    assert gap == Q(1,1800)
    assert gap > e+em
    # Deleting positive pressureless memory from the account changes only
    # rho, so b_minus-b=-m/(2*a^2), without changing the common geometry.
    aa, rr, pp, mm = Q(14,5), Q(7,100), Q(-1,200), Q(10003,10**6)
    assert aa*((rr-mm/aa**3)+3*pp)/2-aa*(rr+3*pp)/2 == -mm/(2*aa**2)
    return {
        'status': 'PASS', 'inherited_endpoint_checks': checkpoints,
        'hold_u_lower': u_lower,
        'chain_rule_upper_bounds_before_relaxation': {
            'v_prime': vp, 'V_prime': Vp, 'v_second': vpp, 'V_second': Vpp,
            'rod_prime': rodp, 'rod_second': rodpp,
            'pressure_prime': Pp, 'pressure_second': Ppp,
            'rho_prime': rhop, 'rho_second': rhopp,
            'a_second': app, 'b_second': bpp,
            'omitted_term_second': lost_bpp,
        },
        'subtraction_readout_error_upper': subtraction_noise,
        'subtraction_error_allowance': Q(1,10**10),
        'source_gap_strict_lower': gap,
        'E_full_plus_E_omitted': e+em,
        'separation_margin_strict_lower': gap-e-em,
        'scope': 'One weak acceleration test; omission changes the bookkeeping, not the solution. Initial constraints and full equations are not reconstructed from this test.'
    }


def bindings(require_all=True):
    paths = {
        'admission.md': HERE/'admission.md', 'proof.md': HERE/'proof.md',
        'check.py': HERE/'check.py', 'results.json': HERE/'results.json',
        'NEXT.md': HERE/'NEXT.md', 'report': HERE.parent/'research_note_1248.md',
        '1247/proof.md': OLD/'proof.md', '1247/results.json': OLD/'results.json',
        'independent_check.py': Path(__file__),
    }
    if require_all:
        assert all(p.is_file() for p in paths.values())
    return {k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items() if p.is_file()}


def run(draft=False):
    first, second = weights(), domain_and_source()
    agreement = 'not yet available in draft'
    if (HERE/'results.json').is_file():
        peer = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        group = peer['group_weights_and_error']
        peer_weights = hashlib.sha256(json.dumps(
            [group['a_weights'], group['b_weights']],
            separators=(',', ':')).encode()).hexdigest()
        assert peer_weights == first['weights_sha256']
        errors = peer['group_common_domain_and_separation']['errors']
        assert Q(errors['total_full']['rational']) == first['E_full']
        assert Q(errors['total_omitted']['rational']) == first['E_omitted']
        assert Q(errors['separation_margin']['rational']) == second['separation_margin_strict_lower']
        agreement = 'All 258 rational weights, both budgets and separation margin agree after independent derivation.'
    return packed({
        'round':1248, 'status':'PASS', 'scientific_groups':2,
        'group_1_exact_interpolation_and_noise':first,
        'group_2_inherited_domain_full_source_and_deletion':second,
        'author_result_crosscheck':agreement,
        'source_sha256':bindings(not draft),
        'CT_adopted':False, 'new_menu_implemented':False,
        'whole_GR_completed':False,
        'method':'Independent Fraction polynomial integration, chain-rule bounds and inherited endpoint checks; no author code import and no repeated 0--2 integration.'
    })


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--draft',action='store_true')
    args=parser.parse_args()
    assert not (args.write and args.draft)
    data=run(args.draft)
    if args.draft:
        print(json.dumps(data,ensure_ascii=False,indent=2))
    elif args.write:
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(data,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('INITIAL_SAVE PASS: 2 independent scientific groups')
    else:
        assert data==json.loads(RESULT.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH PASS: 2 independent scientific groups')
