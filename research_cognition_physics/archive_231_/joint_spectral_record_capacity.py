"""596: spectral hardware boundary and energy-constrained record capacity.

Actual geometry modes are reproduced from the frozen entry. The n^2
spectrum below is a diagnostic example, not the original H0 spectrum.
General spectral and capacity statements are proved in the note.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_spectral_record_capacity_results.json'
spec = importlib.util.spec_from_file_location(
    'spectral_entry596', HERE / 'round596_drafts/spectral_hardware_entry.py')
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)


def geometry_check():
    got = entry.run()
    assert got == json.loads(entry.TARGET.read_text('utf8'))
    return dict(frozen_entry_reproduced=True, geometry_modes=got['geometry_mode_rows'],
                inherited_clock_bounds=got['finite_window_clock_bounds'])


def distinguishability_check():
    # Equal-prior labels j: a shared ground state and a rare distinct mode.
    # These states all have mean relative energy 1, but are not orthogonal.
    rows = []
    for m in (8, 20, 40):
        energies = np.arange(m + 1, dtype=float)**2
        states = np.zeros((m, m + 1))
        for j in range(1, m + 1):
            p = 1 / energies[j]
            states[j-1, 0] = 1-p
            states[j-1, j] = p
        assert np.max(abs(states.sum(axis=1)-1)) < 1e-14
        assert np.max(abs(states @ energies - 1)) < 1e-14
        # All states commute. Basis measurement followed by MAP is optimal.
        assignment = np.argmax(states, axis=0)
        effects = np.zeros_like(states)
        effects[assignment, np.arange(m+1)] = 1
        success = float(np.sum(effects * states) / m)
        expected = (1 - 1/m**2 + sum(1/j**2 for j in range(1, m+1))) / m
        assert abs(success-expected) < 1e-14
        cuts = []
        for R in (8., 25., 64.):
            mask = energies <= R
            rank = int(mask.sum())
            tails = states[:, ~mask].sum(axis=1)
            low_success = float(np.sum(effects[:, mask]*states[:, mask])/m)
            assert low_success <= rank/m + 1e-14
            assert np.max(tails) <= 1/R + 1e-14
            assert success <= rank/m + 2/np.sqrt(R) + 1e-14
            denominator = 1 - .1 - 2/np.sqrt(R)
            cuts.append(dict(cut=R,low_rank=rank,success_upper_bound=rank/m+2/np.sqrt(R),
                             capacity_bound_at_error_point_one=rank/denominator
                             if denominator > 0 else None))
        rows.append(dict(labels=m,mean_relative_energy=1.,optimal_success=success,cuts=cuts))
    return dict(diagnostic_spectrum='E_n=n^2, n=0,...,m; not original H0',rows=rows,
                infinitely_many_distinct_states_not_reliable_infinite_capacity=True)


def run():
    geometry = geometry_check()
    codes = distinguishability_check()
    deps = ('research_note_245.md','research_note_591.md','research_note_595.md',
            'round596_drafts/spectral_hardware_entry.py',
            'round596_drafts/spectral_hardware_entry_results.json')
    return dict(round=596,tests_run=2,failures=0,errors=0,geometry=geometry,codes=codes,
                scope=dict(fixed_original_finite_graph=True,
                           original_eigenvalues_and_optimal_capacity_not_computed=True,
                           no_exact_all_time_embedding_is_analytic=True,
                           finite_precision_capacity_is_analytic=True,
                           finite_window_memory_and_large_system_limits_not_excluded=True),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps})


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(dict(round=596,tests=2,all_passed=True)))
