"""596 entry: actual compact geometry modes and finite-window resources.

The original full spectrum is not computed. Geometry Fourier modes are
normal Gauss states of the original joint model, but are not H0 eigenstates.
The clock tail estimate is reused from round 401, not a new experiment.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_geometry_work_noise as original
import persistent_prefix_history_audit as clock
TARGET=HERE/'spectral_hardware_entry_results.json'


def run():
    hbar=original.HBAR
    inertia=1.3  # Declared admissible geometric inertia, not inferred physics.
    theta=2*np.pi*np.arange(128)/128
    rows=[]
    for K in (1,2,4,8):
        n=np.arange(-K,K+1)
        modes=np.exp(1j*theta[:,None]*n[None,:])
        gram=modes.conj().T@modes/len(theta)
        derivative=1j*n[None,:]*modes
        kinetic=hbar*hbar/(2*inertia)*(derivative.conj().T@derivative)/len(theta)
        expected=hbar*hbar/(2*inertia)*np.diag(n*n)
        error=float(np.linalg.norm(gram-np.eye(len(n))))
        kinetic_error=float(np.linalg.norm(kinetic-expected))
        assert error<2e-14 and kinetic_error<2e-13
        # H_m(theta) and V_g expectations cancel against the same base
        # state because each Fourier multiplier has unit modulus.
        mean=float(np.trace(kinetic).real/len(n))
        analytic=hbar*hbar*K*(K+1)/(6*inertia)
        assert abs(mean-analytic)<2e-14
        rows.append(dict(max_mode=K,orthogonal_modes=len(n),orthogonality_error=error,
                         kinetic_matrix_error=kinetic_error,mean_energy_above_base_packet=mean,
                         analytic_mean_increment=analytic,
                         largest_energy_increment=hbar*hbar*K*K/(2*inertia)))
    tail=[]
    for scaled_end,length in ((1.5,18),(3.,30),(5.,45)):
        bound=clock.truncation_bound(length,scaled_end)
        assert 0<bound<1e-5
        tail.append(dict(scaled_time_window=scaled_end,clock_length=length,
                         existing_round401_vector_error_bound=bound))
    deps=('research_note_591.md','research_note_595.md','research_note_401.md',
          'research_note_578.md','research_note_579.md','joint_geometry_work_noise.py')
    return dict(status='596 spectral compatibility entry; not complete numbered round',checks_passed=True,
                geometry_mode_rows=rows,finite_window_clock_bounds=tail,
                declared_inertia=inertia,hbar=hbar,full_joint_base_energy_not_computed=True,
                no_original_eigenvalue_or_optimal_capacity_claim=True,
                spectral_no_embedding_and_recurrence_are_analytic=True,
                added_hardware_not_embedded_into_original_matter=True,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in deps})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
