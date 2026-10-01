"""Non-numbered method bridge: doubled record kernels of the frozen round 586.

This checks a representation of that instrument, not a new matter model,
full path-integral evolution, continuum limit, or a gravity derivation.
Run without options to audit; --write saves an exclusively new result file.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_source_compression as prior

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'closed_time_path_bridge_entry_results.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kernels(s, separation):
    left_f, _, left_c, _ = prior.instruments(s + separation / 2)
    right_f, _, right_c, _ = prior.instruments(s - separation / 2)
    ff, cc = left_f * right_f, left_c * right_c
    # Keep the actual fine CP map after reporting only parity.
    grouped = np.stack((ff[..., 0] + ff[..., 3], ff[..., 1] + ff[..., 2]), axis=-1)
    return grouped, cc


def run():
    s = np.linspace(-1.2, 1.2, 97)
    f, c = kernels(s, 0.)
    diagonal_error = float(np.max(abs(f-c)))
    assert diagonal_error < 5e-16
    curvature = []
    for step in (.1, .05, .025, .0125):
        f, c = kernels(s, step)
        f, c = f.sum(axis=-1), c.sum(axis=-1)
        a = -2*(f-c)/step**2
        log_a = -2*(np.log(f)-np.log(c))/step**2
        error = float(np.max(abs(a-prior.delta_a(s))))
        log_error = float(np.max(abs(log_a-prior.delta_a(s))))
        curvature.append(dict(separation=step, kernel_curvature_error=error,
                              log_kernel_curvature_error=log_error))
    assert curvature[-1]['kernel_curvature_error'] < 3e-6
    assert curvature[-1]['log_kernel_curvature_error'] < 3e-6
    assert curvature[-1]['kernel_curvature_error'] < curvature[0]['kernel_curvature_error']/40

    # Integrate the curvature of the new representation against the old packet;
    # compare to saved 586 energy, without rerunning old propagation tests.
    pp = prior.prior_packet
    q, qw = np.polynomial.legendre.leggauss(80)
    x, y = np.meshgrid(q, q, indexing='ij')
    h = pp.HCENTER + pp.HRADIUS*x
    s = pp.SCENTER + pp.SRADIUS*y
    chi = np.exp(-1/(1-x*x)-1/(1-y*y))
    prob = qw[:, None]*qw[None, :]*pp.HRADIUS*pp.SRADIUS*h**3*chi**2
    prob /= prob.sum()
    kss = (prior.M-(h*h+s*s)/6)*(1-s*s/(6*prior.M))
    saved = json.loads((HERE/'joint_record_source_compression_results.json').read_text('utf8'))
    expected = saved['evidence']['original_compact_Gauss_packet']['rows'][-1]['energy_gap']
    energy_rows = []
    for step in (.05, .025, .0125, .00625):
        f, c = kernels(s, step)
        a = -2*(f.sum(axis=-1)-c.sum(axis=-1))/step**2
        energy = prior.HBAR**2/(2*prior.WEIGHT)*float(np.sum(prob*kss*a))
        energy_rows.append(dict(separation=step, energy_from_kernel_curvature=energy,
                               error_against_frozen_round586=abs(energy-expected)))
    assert energy_rows[-1]['error_against_frozen_round586'] < 3e-7
    assert energy_rows[-1]['error_against_frozen_round586'] < energy_rows[0]['error_against_frozen_round586']/40
    protected = ['research_note_586.md', 'joint_record_source_compression.py',
                 'joint_record_source_compression_results.json', 'research_round_586_checks.json',
                 'unified_physics_condition_ledger_586.md', 'joint_curved_quantum_source.py',
                 'joint_quantum_measure_records.py', 'joint_local_lapse_source.py']
    protected += [p.relative_to(HERE).as_posix() for p in sorted((HERE/'round587_drafts').iterdir()) if p.is_file()]
    return dict(kind='non-numbered method bridge', date='2026-10-01', numbered_round=586,
                cumulative_numbered_checks_unchanged=2989, checks_passed=True,
                diagonal_probability_error=diagonal_error, curvature=curvature,
                frozen_round586_energy_gap=expected, energy_from_curvature=energy_rows,
                full_graph_path_integral_not_simulated=True, dynamic_geometry_not_derived=True,
                source_hashes={name: digest(HERE/name) for name in protected})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    elif TARGET.exists():
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('checks_passed', 'diagonal_probability_error',
                     'curvature', 'energy_from_curvature')}, indent=2))
