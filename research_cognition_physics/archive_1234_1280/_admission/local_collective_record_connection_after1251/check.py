"""Exact algebra/support certificate; no ODE, implementation or new research round."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def run():
    project = next(p for p in Path(__file__).resolve().parents
                   if (p / 'research_cognition_physics').is_dir())
    stage = project / 'research_cognition_physics/archive_1234_'
    sources = {rel: hashlib.sha256((stage / rel).read_bytes()).hexdigest()
               for rel in ('1247/proof.md', '1249/canonical.md', '1251/proof.md')}
    # Q=(Y1+Y2,Y1-Y2), P=((K1+K2)/2,(K1-K2)/2).
    A = [[F(1), F(1)], [F(1), F(-1)]]
    B = [[F(1,2), F(1,2)], [F(1,2), F(-1,2)]]
    canonical = [[sum(A[i][k]*B[j][k] for k in range(2))
                  for j in range(2)] for i in range(2)]
    assert canonical == [[1,0],[0,1]]
    # (K1^2+K2^2)/2 = Kbar^2 + R^2, including the mixed coefficient.
    energy_coefficients = [sum(v*v for v in column)/2
                           for column in ([F(1),F(1)], [F(1),F(-1)])]
    mixed_coefficient = F(1)*F(1)+F(1)*F(-1)
    assert energy_coefficients == [1,1] and mixed_coefficient == 0

    d, s = F(1), F(1,16)
    mixtures = {'A': [(F(1,2),d), (F(1,2),-d)],
                'B': [(F(3,4),F(0)), (F(1,8),2*d), (F(1,8),-2*d)]}
    moments, supports, probabilities = {}, {}, {}
    for name, mixture in mixtures.items():
        assert sum(w for w, _ in mixture) == 1
        assert all(w > 0 for w, _ in mixture)
        moments[name] = {degree: sum(w*c**degree for w,c in mixture)
                         for degree in (1,2,4)}
        supports[name] = []
        probability = F(0)
        for weight, centre in mixture:
            lo, hi = centre-2*s, centre+2*s
            inside = -d/2 < lo and hi < d/2
            outside = hi < -d/2 or lo > d/2
            assert inside or outside  # No integration of unspecified bump is needed.
            probability += weight if inside else 0
            supports[name].append({'weight': str(weight), 'lo': str(lo),
                                   'hi': str(hi), 'inside_local_event': inside})
        probabilities[name] = probability
    assert s < d/8
    assert moments['A'][1] == moments['B'][1] == 0
    assert moments['A'][2] == moments['B'][2] == d*d
    assert moments['B'][4]-moments['A'][4] == 3*d**4
    assert probabilities == {'A': F(0), 'B': F(3,4)}
    # For a real even common bump, mean0,var=v_s,fourth=m4_s:
    # translated R^4 has coefficients (1,6*Ecentre^2,Ecentre^4).
    fourth_coeff = {name: [F(1), 6*moments[name][2], moments[name][4]]
                   for name in mixtures}
    assert fourth_coeff['A'][:2] == fourth_coeff['B'][:2]

    mu, coupling, fmax = F(1), F(1,2), F(1)
    remainder = 1/mu-coupling**2*fmax**2
    assert remainder == F(3,4) and remainder > 0
    # Exact completion of square in K,k_d, checked coefficient by coefficient.
    assert [F(1,2), coupling*fmax, coupling**2*fmax**2/2+remainder/2] == [F(1,2), F(1,2), F(1,2)]
    return {
        'status': 'PASS_EXACT_INTERFACE_SCREEN_ONLY',
        'sources_sha256': sources,
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'canonical_pair_matrix': [[str(x) for x in row] for row in canonical],
        'energy_coefficients_Kbar2_R2': [str(x) for x in energy_coefficients],
        'd': str(d), 'bump_half_width_s': str(s),
        'centre_moments': {name: {str(k): str(v) for k,v in moment.items()}
                           for name,moment in moments.items()},
        'fourth_moment_coefficients_m4_v_1': {name: [str(x) for x in values]
                                             for name,values in fourth_coeff.items()},
        'local_supports': supports,
        'exact_local_probabilities': {name: str(p) for name,p in probabilities.items()},
        'qnd_quadratic_remainder': str(remainder),
        'scope': 'Normal-state proof uses C_c infinity bump in Markdown. No complete parent embedding, source identity or local instrument is certified.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true', help='Exclusive first save only.')
    args = parser.parse_args()
    result = run()
    destination = Path(__file__).with_name('results.json')
    if args.write:
        with destination.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('SAVED_EXACT_INTERFACE_SCREEN_ONLY')
    else:
        assert result == json.loads(destination.read_text(encoding='utf-8')), 'Saved result mismatch'
        print('READ_ONLY_MATCH_EXACT_INTERFACE_SCREEN_ONLY')
