"""744: original joint reflection/fermion phase selects an actual readout.

Full-time POVM conclusions are analytic consequences of the exact symmetry.
Numerics check the original Hamiltonian coefficients, curved target and the
actual encoded force. No bosonic time propagation is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_native_quantum_record as prior
old=prior.old;matter=old.matter;geom=matter.original
TARGET=HERE/'native_instrument_selection_entry_results.json'

def run():
    rng=np.random.default_rng(744);reflection=np.diag([1.,1.,1.,1.,-1.])
    mass_error=[];geometry_error=[]
    for _ in range(12):
        phi=rng.normal(size=5)*.35;chi=rng.normal(size=5)*.35
        h,d=matter.mass_matrices(phi);hh,dd=matter.mass_matrices(reflection@phi)
        mass_error.extend([np.max(abs(hh-h)),np.max(abs(dd+d))])
        geometry_error.extend([abs(geom.node_potential(phi)-geom.node_potential(reflection@phi)),
            abs(geom.distance_squared(phi,chi)-geom.distance_squared(reflection@phi,reflection@chi)),
            np.max(abs(geom.metric(reflection@phi)-reflection@geom.metric(phi)@reflection))])
    assert max(mass_error+geometry_error)<3e-13
    # Native complete32-mode coefficient, projected onto the actual even code.
    h,d=old.mass_x(np.eye(5)[4]);phase=np.exp(1j*np.angle(matter.Y['s']))
    basis=[{0:1.},{(1<<30)|(1<<31):phase}]
    Q=np.array([[old.dot(a,old.car.quadratic(b,h,d)) for b in basis] for a in basis])
    X=np.array([[0,1],[1,0]],complex);Z=np.diag([1.,-1.])
    assert np.max(abs(Q-abs(matter.Y['s'])*X))<2e-15
    assert np.max(abs(Z@Q@Z+Q))<2e-15
    packet=prior.packet(80);w=1.
    effect_t2=-packet['inherited_sin_readout_force_weight']*Q/(8*w)
    assert abs(effect_t2[0,1])>1e-3
    assert np.max(abs(np.diag(effect_t2)))<2e-15
    deps=('research_note_577.md','research_note_598.md','research_note_717.md','research_note_743.md',
          'joint_fermion_gauss_completion.py','joint_curved_quantum_source.py','joint_native_quantum_record.py')
    return dict(entry_round=744,latest_completed_round=743,formal_test_count_unchanged=3424,
        original_complete_species_mass_symmetry_error=float(max(mass_error)),
        original_potential_edge_metric_symmetry_error=float(max(geometry_error)),
        encoded_Majorana_matrix_real=Q.real.tolist(),encoded_Majorana_matrix_imag=Q.imag.tolist(),
        symmetric_ready_packet=packet,
        actual_plus_effect_t2_real=effect_t2.real.tolist(),actual_plus_effect_t2_imag=effect_t2.imag.tolist(),
        exact_effect_relation='Z M_plus(t) Z = I - M_plus(t) for symmetric ready state and original global singlet reflection with exp(i pi N/2).',
        exact_equal_probabilities_for_encoded_empty_and_pair='1/2 for every t; conditional on symmetry and fixed original sin-s effect.',
        no_claim_that_all_observables_or_all_preparations_are_population_blind=True,
        no_claim_of_exact_occupation_instrument_or_autonomous_final_readout=True,
        numerical_time_evolution_performed=False,
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
