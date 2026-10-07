"""810: action/Hamiltonian source sign and source-splice balance calibration.

Uses original mass matrices and the 809 CAR reflection. The scalar temporal
balance check is an operator-identity diagnostic, not the full gravity PDE.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_source_contract_results.json'
import readout_source_probe as early
old=early.old
def run():
    evidence=early.run()
    N=old.occupied(np.zeros(3));P=old.I-N
    f=np.zeros(64,complex);f[30]=1/np.sqrt(2);f[31]=1j/np.sqrt(2)
    cf=old.CHARGE@f.conj();F=np.outer(f,f.conj())+np.outer(cf,cf.conj())
    R=old.I-2*F;deltaP=(R@P@R-P)/2;deltaN=-deltaP
    expanded=-F@P-P@F+2*F@P@F
    finite_rank=int(np.linalg.matrix_rank(deltaP,tol=1e-12));assert finite_rank<=4
    assert np.max(abs(deltaP-expanded))<1e-13
    # At this auxiliary zero-momentum block S_F density has -Psi* M Psi/2.
    # The old working label "mass_source_changes" contains +dH derivatives.
    def relative_energy(phi):return float((np.trace(old.vertex.mass(phi)@deltaN)/2).real)
    rows=[];step=2e-5
    for a,d in enumerate(np.eye(5)):
        numerical=(relative_energy(old.PHI+step*d)-relative_energy(old.PHI-step*d))/(2*step)
        analytic=evidence['original_five_mass_source_changes'][a]
        assert abs(numerical-analytic)<2e-10
        rows.append(dict(direction=a,Hamiltonian_derivative=analytic,
            action_mass_density_derivative=-analytic,finite_difference_error=abs(numerical-analytic)))
    # J_pre and J_post separately obey d_t J=0. Their smooth splice does not.
    jump=evidence['diagnostic_auxiliary_free_energy_change']
    widths=[.5,.25,.125];balance=[]
    for width in widths:
        t=np.linspace(-width,width,1001);x=t/width
        # cubic smooth transition on this interval, constant outside.
        chi=.5+.75*x-.25*x**3;dchi=(.75-.75*x*x)/width
        J=chi*jump;defect=dchi*jump
        integrated=float(np.trapezoid(defect,t))
        exact=jump
        assert abs(integrated-exact)<3e-7
        balance.append(dict(width=width,integrated_balance_defect=integrated,
            endpoint_jump=float(J[-1]-J[0]),max_defect=float(max(abs(defect)))))
    return dict(round=810,all_checks_passed=True,
        original_Nambu_reflection_covariance_rank=finite_rank,
        covariance_formula_residual=float(np.max(abs(deltaP-expanded))),
        source_sign_checks=rows,temporal_splice_calibration=balance,
        old_working_source_label_interpretation='The five stored values are Hamiltonian mass derivatives; action-density derivatives have the opposite sign.',
        full_original_spacetime_source_numerically_computed=False,
        source_Ward_and_relative_Cauchy_connection_proof_in_note=True,
        instrument_controller_or_closed_resources_implemented=False,
        new_numbered_test_groups=1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
