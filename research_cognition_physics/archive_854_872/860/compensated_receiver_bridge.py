"""860: the compensated shared-reference branch, constrained symbol and code.
Continuum existence/stability is proved in research_note_860.md. Finite causal
kernels below are explicit diagnostics, not numerical gravity solutions.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
import compensated_probe_kernel_probe as probe
sys.path.insert(0,str(HERE.parent/'853'))
import native_receiver_content as receiver
TARGET=HERE/'compensated_receiver_bridge_results.json'

def extended_complex():
    with ResearchRuntime(Layout()).installed():
        import joint_covariant_gauge_complex as old
        phi=np.array([.13,.64,-.11,.08,.37]);weights=np.repeat([.73,1.17,.41],[8,3,1])
        rows=[]
        for covector in ((2.,.1,.3,-.2),(1.,1.,0.,0.),(.2,.3,-.7,1.)):
            H,B,K,Ka,P,k2=old.principal(covector,phi,weights,kappa=.81)
            H6=np.zeros((64,64));H6[:63,:63]=H;H6[63,63]=1.
            K6=np.vstack((K,np.zeros((1,16))))
            Ka6=-np.linalg.solve(B,K6.T@H6)
            P6=np.zeros((64,64));P6[:63,:63]=P;P6[63,63]=k2
            D6=P6-K6@Ka6;D0=-Ka6@K6
            errors=dict(full_wave=float(np.max(abs(D6-k2*np.eye(64)))),gauge_wave=float(np.max(abs(D0-k2*np.eye(16)))),noether=float(np.max(abs(P6@K6))),adjoint=float(np.max(abs(H6@P6-P6.T@H6))))
            assert max(errors.values())<8e-15
            rows.append(dict(covector=list(covector),lorentz_square=k2,errors=errors))
        eig=np.linalg.eigvalsh(H6);signature=[int(sum(eig>1e-10)),int(sum(eig< -1e-10))]
        assert signature==[48,16]
        # The new scalar adds a zero-order diffeomorphism row, not a principal
        # gauge derivative. Its adjoint is a genuine source contribution.
        q,psi,_,_,_=old.bg.completed(12,1.)
        point=(0,3,1);epsilon=.02;grad=np.array([0.,epsilon/psi[point]**2,0.,0.])
        gauge0=np.zeros((64,16));gauge0[63,:4]=grad
        adjoint0=np.linalg.solve(B,gauge0.T@H6)
        expected=2*.81*old.ETA@grad
        error=float(np.max(abs(adjoint0[:4,63]-expected)));assert error<1e-15
        assert np.max(abs(adjoint0[:4,63]))>0
    return dict(field_components=64,gauge_components=16,field_pairing_signature=signature,physical_scalar_components=6,rows=rows,probe_gauge_adjoint_error=error,probe_gauge_adjoint_diffeomorphism_column=adjoint0[:4,63].tolist(),scope='Extended original principal symbol and a zero-order probe adjoint calibration; not the full lower-order PDE or the new 859 background integration.')

def content_bridge():
    previous=receiver.run();assert previous==json.loads(receiver.TARGET.read_text('utf-8'))
    # Original CAR probabilities, not independently chosen classical input bits.
    p0=list(map(F,previous['original_code_joint_probabilities']['0']))
    p1=list(map(F,previous['original_code_joint_probabilities']['3/5']))
    signs=list(receiver.itertools.product((-1,1),repeat=3))
    triple=sum(F(int(np.prod(a)))*(y-x) for a,x,y in zip(signs,p0,p1));assert triple==F(27,50)
    leading0=F(previous['exact_hbar_cubed_leading_coefficient']);assert leading0==F(-9,1960)
    base=probe.kernel(F(0))[0];rows=[]
    for epsilon in (F(1,50),F(1,25),F(1,10)):
        full,omitted,_=probe.kernel(epsilon);ratio=full/base
        coefficient=leading0*ratio**3;wrong=leading0*(omitted/base)**3
        assert coefficient!=0 and wrong!=coefficient
        assert probe.kernel(-epsilon)[0]/base==ratio
        # In the native receiver's exact finite formula, Gaussian damping is
        # 1+O(hbar), while each sin starts at hbar: covariance enters after hbar^3.
        amplitudes=np.array([float(F(x)) for x in previous['source_amplitudes']])*float(ratio)
        lambdas=np.array([float(F(x)) for x in previous['lambdas']]);etas=np.array([float(F(x)) for x in previous['etas']])
        noise_rows=[]
        for variances in (np.array([.2,.3,.4]),np.array([4.,5.,6.])):
            h=.00001
            value=-float(triple)*np.prod(np.exp(-etas**2*h*variances/2)*np.sin(etas*lambdas*amplitudes*h))
            normalized=value/h**3
            assert abs(normalized-float(coefficient))<3e-7
            noise_rows.append(dict(variances=variances.tolist(),finite_hbar=h,coefficient_after_hbar_cubed=float(normalized)))
        rows.append(dict(epsilon=str(epsilon),response_ratio=str(ratio),exact_leading_coefficient=str(coefficient),leading_float=float(coefficient),omitted_reference_source_coefficient=str(wrong),noise_calibrations=noise_rows,even_in_background_preparation=True))
    return dict(original_source_content_difference=str(triple),original_853_coefficient=str(leading0),rows=rows,actual_curved_receiver_coefficient_numerically_evaluated=False,finite_diagnostic_not_full_continuum=True)

def run():
    working=probe.run();assert working==json.loads(probe.TARGET.read_text('utf-8'))
    return dict(round=860,all_checks_passed=True,fresh_test_groups=1,
        earlier_working_results_reproduced=True,extended_gauge_complex=extended_complex(),shared_kernel_and_original_record=content_bridge(),
        analytic_scope='For sufficiently small nonzero prepared probe backgrounds, a declared old-material calibrated compensation keeps local new reference rank and common linearized cones, while the original six-coupling leading record coefficient remains nonzero. Positive local formal branches use the full mixed free process and all sources; no fixed full-state equivalence or finite-coupling convergence is asserted.',
        newly_declared_compensation_profile_and_action=True,existing_probe_and_receivers_used=True,
        parameter_continuity_and_constraint_surjectivity_are_analytic_not_numerical=True,
        bosonic_Hadamard_family_not_needed_for_leading_record_coefficient=True,
        exact_cones_of_full_nonlinear_or_quantum_process_proved=False,
        original_graph_to_continuum_proved=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in r.items() if k!='shared_kernel_and_original_record'},ensure_ascii=False,indent=2))
