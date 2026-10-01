"""575 entry probe, not a completed research round.
Random phase components of a Gaussian channel are not actual pointer-Q outcomes.
Einstein-constraint completion below is conditional classical initial data.
"""
from pathlib import Path
import json
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_gauss_einstein_initial_data as old
import joint_curved_quantum_source as quantum

TARGET=HERE/'momentum_instrument_entry_probe_results.json'


def run():
    q=old.make_source(16);psi,_,_=old.solve_hamiltonian(q);dx=q['dx']
    x,y=np.moveaxis(q['grid'],-1,0)[:2]
    rows=[]
    for label,weight in (('collective',np.ones_like(x)),('fixed_profile',1+.4*np.cos(x+y))):
        direction=np.zeros_like(q['p']);direction[...,:4]=weight[...,None]*q['phi'][...,:4]
        # No singlet or electric kick. All Higgs kicks are radial.
        DX=q['Dphi'];deltaM=np.einsum('...a,...ia->...i',direction,DX)
        total=dx**3*np.sum(deltaM,axis=(0,1,2))
        hstar=np.sqrt(old.old.PAR['h2'])
        expected=np.zeros(3) if label=='collective' else old.VOL*.01*hstar*hstar*np.array([1.,1.,0.])
        assert np.max(abs(total-expected))<1e-12
        P_complex=direction[...,:2]+1j*direction[...,2:4]
        qw,q0=old.old.charges(q['f']['X'],P_complex)
        assert max(np.max(abs(qw)),np.max(abs(q0)))<1e-14
        Kinv=quantum.inverse(q['phi'])
        G=np.einsum('...i,...ij,...j->...',direction,Kinv,direction)
        quadratic=float(dx**3*np.sum(psi**-6*G)/2)
        linear=float(dx**3*np.sum(psi**-6*np.einsum('...i,...ij,...j->...',q['p'],Kinv,direction)))
        row=dict(profile=label,total_momentum_per_phase=total.tolist(),analytic_total=expected.tolist(),
                 energy_linear_in_phase=linear,energy_quadratic_in_phase=quadratic,
                 nonzero_phase_flat_periodic_CMC_obstruction=label=='fixed_profile')
        if label=='collective':
            trials=[]
            for phase in (-.05,0.,.05):
                new=dict(q);new['p']=q['p']+phase*direction;new['mom']=q['mom']+phase*deltaM
                new['pKp']=np.einsum('...i,...ij,...j->...',new['p'],Kinv,new['p'])
                fitted,_,info=old.solve_hamiltonian(new)
                assert info['original_Hamiltonian_residual']<3e-8
                trials.append(dict(phase=phase,max_psi_change=float(np.max(abs(fitted-psi))),**info))
            row['separate_classical_completions']=trials
        rows.append(row)
    return dict(status='entry probe only; no new completed round',source='same 572 source',
        rows=rows,phase_labels_are_not_pointer_position_outcomes=True,
        individual_Gauss_legal_does_not_imply_same_geometry_or_CMC_compatibility=True,
        apparatus_recoil_and_stress_not_included=True,
        no_rule_assigning_geometry_to_quantum_unravelings=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    a=p.parse_args();result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
