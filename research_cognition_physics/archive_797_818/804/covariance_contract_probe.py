"""804 working: map the old occupied covariance to the new ordered convention.

The exact CAR dictionary is general. A constant original mass subblock provides
an independent finite Fock check; it is not the nonstationary continuum state.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'803'))
import quasifree_pairing_probe as old
vertex=old.vertex
TARGET=HERE/'covariance_contract_probe_results.json'

def run():
    ann=old.car(4);fields=ann+[a.conj().T for a in ann]
    def quad(k):return sum((fields[i].conj().T@fields[j]*k[i,j]/2
                           for i in range(8) for j in range(8)),np.zeros((16,16),complex))
    selected=[24,25,30,31];indices=selected+[i+32 for i in selected]
    phi=np.array([0.,.55,0.,0.,.4])
    mass=vertex.mass(phi)[np.ix_(indices,indices)]
    energy,vec=np.linalg.eigh(quad(mass));vacuum=vec[:,0]
    expect=lambda a:vacuum.conj()@a@vacuum
    greater=np.array([[expect(a@b.conj().T) for b in fields] for a in fields])
    occupied=np.array([[expect(b.conj().T@a) for b in fields] for a in fields])
    ev,unitary=np.linalg.eigh(mass)
    negative=unitary[:,ev<0]@unitary[:,ev<0].conj().T
    errors=dict(canonical_anticonmutator=float(np.max(abs(greater+occupied-np.eye(8)))),
        occupied_negative_spectrum=float(np.max(abs(occupied-negative))),
        greater_positive_spectrum=float(np.max(abs(greater-(np.eye(8)-negative)))),
        ground_energy=float(abs(expect(quad(mass))-.5*np.trace(mass@occupied))))
    rows=[]
    f=np.eye(8)[:,2];cf=np.eye(8)[:,6]
    d=-np.outer(f,f)+np.outer(cf,cf)
    record=np.eye(16)-ann[2].conj().T@ann[2]
    for a in range(5):
        k=-vertex.dmass(phi,np.eye(5)[a])[np.ix_(indices,indices)]
        direct=expect(record@quad(k))-expect(record)*expect(quad(k))
        correct=np.trace(occupied@d@greater@k)/2
        wrong=np.trace(greater@d@occupied@k)/2
        assert abs(direct-correct)<1e-12
        rows.append(dict(original_scalar_direction=a,
            correct_ordered_covariance=[float(correct.real),float(correct.imag)],
            occupied_as_greater_defect=float(abs(direct-wrong))))
    assert max(errors.values())<2e-12
    wrong_energy=float((.5*np.trace(mass@greater)).real)
    assert abs(wrong_energy-energy[0])>.1
    return dict(working_round=804,all_checks_passed=True,
        original_mass_subblock_indices=selected,
        mapping='P_803_greater = I - N_730_occupied, in the same [c,c_dagger] basis',
        maximum_residuals=errors,actual_ground_energy=float(energy[0]),
        wrong_energy_if_occupied_and_greater_confused=wrong_energy,
        five_original_mass_derivatives=rows,
        same_state_preserved=True,new_reference_selected=False,
        continuous_original_PDE_response_evaluated=False,
        finite_diagnostic_is_original_W=False,formal_round_completed=False,new_numbered_test_groups=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
