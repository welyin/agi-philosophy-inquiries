"""Fixed-summary CKM outer-set check. Default is read-only; --write is exclusive.

This is neither a joint likelihood nor a new statistical theorem. Decimal
outward padding controls arithmetic; publication rounding is a separate input
envelope and does not bound missing systematics or lattice-scale rematching.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
from math import erfc, sqrt
import argparse
import json

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / 'public_inputs.json').read_text(encoding='utf-8'))
PAD = D('1e-38')


def show(value):
    return format(value, '.24g')


def outward(value, upper):
    return str(value.quantize(D('1e-24'), rounding=ROUND_CEILING if upper else ROUND_FLOOR))


def row_region(box):
    """Positive CKM magnitudes; ratio strip plus three-component row norm."""
    xl, xu = box['x']; yl, yu = box['y']; rl, ru = box['r']
    zl, zu = box['z']
    assert min(xl, yl, rl, zl) >= 0 and ru > 0
    a = max(xl, yl / ru)
    b = min(xu, yu / rl) if rl else xu
    if a > b:
        return {'feasible': False, 'reason': 'ratio_empty',
                'x_projection': [show(a), show(b)]}
    low = (a, max(yl, rl*a))
    high = (b, min(yu, ru*b))
    smin = sum(v*v for v in low)
    smax = sum(v*v for v in high)
    nmin, nmax = smin+zl*zl, smax+zu*zu
    result = {'feasible': nmin <= 1 <= nmax,
              'x_projection': [show(a), show(b)],
              'Smin': outward(smin, False), 'Smax': outward(smax, True),
              'norm_min': outward(nmin, False), 'norm_max': outward(nmax, True),
              'deficit_at_max': outward(1-nmax, False)}
    if nmin > 1:
        result['reason'] = 'norm_above_one'
    elif nmax < 1:
        result['reason'] = 'norm_below_one'
    else:
        # The monotone segment in the convex ratio strip covers all S values.
        target = max(smin, 1-zu*zu)
        target = min(max(target, 1-D(DATA['z']['mean'])**2),
                     min(smax, 1-zl*zl))
        dx, dy = high[0]-low[0], high[1]-low[1]
        aa = dx*dx+dy*dy
        bb = 2*(low[0]*dx+low[1]*dy)
        delta = target-smin
        if not delta:
            t = D(0)
        elif aa:
            t = 2*delta/(bb+(bb*bb+4*aa*delta).sqrt())
        else:
            t = delta/bb
        x, y = low[0]+t*dx, low[1]+t*dy
        assert x > 0
        witness = {'x': x, 'y': y, 'r': y/x,
                   'z': max(D(0), 1-target).sqrt()}
        assert all(l-PAD <= witness[k] <= u+PAD for k,(l,u) in box.items())
        assert abs(witness['x']**2+witness['y']**2+witness['z']**2-1) < PAD
        assert abs(witness['y']-witness['x']*witness['r']) < PAD
        result.update(reason='outer_set_nonempty',
                      witness={k: show(v) for k,v in witness.items()})
    return result


def make_box(contract, rounded, zmode):
    box = {}
    for key in ('x', 'y', 'r', 'z'):
        item = DATA[key]
        mean = D(item['mean'])
        half = D(item['rounding_half_unit']) if rounded else D(0)
        if 'sigma_components' in item:
            components = [D(s)+half for s in item['sigma_components']]
            sigma = sum(s*s for s in components).sqrt()
            if key == 'z' and zmode == 'linear_widening':
                sigma = sum(components)
        else:
            sigma = D(item['sigma'])+half
        if contract == 'supplement_three_two_sided' and key == 'r':
            # No r measurement is used in this subtest. Finite bound 1 is
            # nonbinding throughout these data boxes (verified below).
            box[key] = (D(0), D(1))
            continue
        info = DATA['contracts'][contract]
        q = D(info['q']) if 'q' in info else D(
            info['ratio_two_sided_q'] if key == 'r' else info['upper_q'])
        upper = mean+half+q*sigma+PAD
        lower = max(D(0), mean-half-q*sigma-PAD)
        if contract == 'supplement_directional' and key != 'r':
            lower = D(0)
        box[key] = (lower, min(D(1), upper))
    if contract == 'supplement_three_two_sided':
        assert box['y'][1] < box['x'][0]
    return box


def solver_checks():
    cases = [
        ({'x':('0.8','0.9'),'y':('0.1','0.2'),'r':('0.5','0.6'),'z':('0','0.1')}, 'ratio_empty'),
        ({'x':('0.9','1'),'y':('0.6','0.7'),'r':('0.6','0.9'),'z':('0','0.1')}, 'norm_above_one'),
        ({'x':('0.3','0.4'),'y':('0.1','0.2'),'r':('0.2','0.8'),'z':('0','0.1')}, 'norm_below_one'),
        ({'x':('0.6','0.6'),'y':('0.8','0.8'),'r':('1','1.5'),'z':('0','0')}, 'outer_set_nonempty'),
    ]
    for b, reason in cases:
        answer = row_region({k: tuple(map(D,v)) for k,v in b.items()})
        assert answer['reason'] == reason
    return len(cases)


def calculate(precision=75):
    with localcontext() as ctx:
        ctx.prec = precision
        qtail = lambda q: erfc(q/sqrt(2))/2
        coverage = {
            'main_four_two_sided': 1-8*qtail(2.5),
            'supplement_three_two_sided': 1-6*qtail(2.4),
            'supplement_directional': 1-3*qtail(2.25)-2*qtail(2.5),
        }
        assert all(v > .95 for v in coverage.values())
        results = {}
        for mode in DATA['z_sigma_modes']:
            for rounded in (False, True):
                tag = f'zsigma_{mode}_publication_rounding_{str(rounded).lower()}'
                results[tag] = {}
                for contract in DATA['contracts']:
                    box = make_box(contract, rounded, mode)
                    result = row_region(box)
                    result['box'] = {k: [outward(l,False),outward(u,True)] for k,(l,u) in box.items()}
                    results[tag][contract] = result
                assert results[tag]['main_four_two_sided']['feasible']
                assert not results[tag]['supplement_directional']['feasible']
                assert results[tag]['supplement_three_two_sided']['feasible'] == rounded
        return {
            'round': 1061, 'all_checks_passed': True,
            'scope': 'retrospective conditional summary outer regions, not a common nuisance fit',
            'coverage_lower_bounds_under_adopted_marginals':
                {k: format(v,'.16g') for k,v in coverage.items()},
            'arithmetic_outward_padding': str(PAD),
            'boundary_solver_cases': solver_checks(),
            'results': results,
            'main_four_summary_region_rejected': False,
            'full_shared_nuisance_or_scale_matching_passed': False,
            'direction_or_subset_selection_adjusted': False,
            'new_physics_discovery_claimed': False,
            'whole_roadmap_complete': False,
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    assert result == calculate(100)
    path = HERE / 'results.json'
    if args.write:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        assert result == json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps({'round':1061,'all_checks_passed':True,
        'contracts':12,'main_outer_set_nonempty':True,
        'publication_rounding_changes_three_summary_subtest':True}))
