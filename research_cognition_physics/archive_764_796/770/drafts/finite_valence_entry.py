"""770 working calibration: one-loop finite valence, not loop dynamics."""
import argparse
import itertools
import json
from pathlib import Path

TARGET = Path(__file__).with_name('finite_valence_entry_results.json')


def run():
    rows = []
    for external in range(1, 5):
        # At one loop sum(n-2)V_n=E: no valence exceeds E+2.
        valences = list(range(3, external+3))
        configurations = []
        for counts in itertools.product(range(external+1), repeat=len(valences)):
            if sum((n-2)*v for n, v in zip(valences, counts)) != external:
                continue
            vertices = sum(counts)
            twice_internal = sum(n*v for n, v in zip(valences, counts))-external
            assert twice_internal % 2 == 0
            internal = twice_internal//2
            assert internal-vertices+1 == 1
            configurations.append({str(n): v for n, v in zip(valences, counts) if v})
        rows.append(dict(external_legs=external, configurations=configurations))
    assert rows[0]['configurations'] == [{'3': 1}]
    assert rows[1]['configurations'] == [{'4': 1}, {'3': 2}]
    assert [len(r['configurations']) for r in rows] == [1, 2, 3, 5]
    return dict(working_round=770, completed_round=769,
                new_numbered_scientific_tests=0, rows=rows,
                scope='Necessary connected graph valence identities; no vertex tensor, allowed species contraction, loop integral or anomaly computed.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--write-results', action='store_true')
    args = p.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('working_round', 'completed_round', 'new_numbered_scientific_tests')}))
