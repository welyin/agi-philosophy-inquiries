"""580: same H5 matter, scalar one-loop local terms, and metric variation.

Continuum background-field quantization is an input, not an established limit
of the graph. Euclidean action derivatives are not Lorentz stress tensors.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_scalar_effective_geometry_results.json'
spec=importlib.util.spec_from_file_location('scalar_effective_entry',HERE/'round580_drafts/scalar_effective_action_entry_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
ALPHA=1/27
BETA=1/54


def periodic_nonclosure_check():
    # Unit transverse volume, periodic x in [0,2pi]. q is target geodesic length.
    x=np.linspace(0,2*np.pi,4096,endpoint=False)
    amplitude=.24
    rows=[]
    for n in (1,2,3,5):
        q=amplitude*np.sin(n*x)
        dq=amplitude*n*np.cos(n*x)
        qprime4=float(2*np.pi*np.mean(dq**4)/18)
        exact=np.pi*amplitude**4*n**4/24
        # Nonconstant representative U and K restricted to the same geodesic.
        potential=float(2*np.pi*np.mean(np.exp(q)+q**6))
        kinetic=float(2*np.pi*np.mean((1+q*q)*dq*dq))
        assert abs(qprime4-exact)<2e-15
        rows.append(dict(n=n,four_gradient_integral=qprime4,analytic_integral=float(exact),
                         example_potential_integral=potential,example_kinetic_integral=kinetic))
    assert max(abs(r['example_potential_integral']-rows[0]['example_potential_integral']) for r in rows)<1e-14
    assert max(abs(r['example_kinetic_integral']/r['n']**2-rows[0]['example_kinetic_integral']) for r in rows)<1e-14
    # Three distinct frequencies already preclude equality with c0+c2*n^2.
    first_two=np.array([[1.,r['n']**2] for r in rows[:2]])
    fitted=np.linalg.solve(first_two,np.array([r['four_gradient_integral'] for r in rows[:2]]))
    residual=rows[2]['four_gradient_integral']-float(np.array([1.,9.])@fitted)
    assert residual>1e-3
    return dict(amplitude=amplitude,rows=rows,third_frequency_residual_after_two_parameter_fit=residual,
                analytic_no_absorption_holds_for_arbitrary_fixed_U_and_K=True,
                off_shell_test_family_not_high_frequency_physical_prediction=True)


def four_gradient_density(inverse_metric,B):
    S=float(np.trace(inverse_metric@B))
    Q=float(np.trace(inverse_metric@B@inverse_metric@B))
    return float((ALPHA*S*S+BETA*Q)/np.sqrt(np.linalg.det(inverse_metric)))


def four_gradient_variation_check():
    rng=np.random.default_rng(580014)
    rows=[]
    for index in range(3):
        Z=rng.normal(size=(4,4))*.15
        ginv=np.eye(4)+Z@Z.T
        metric=np.linalg.inv(ginv)
        X=rng.normal(size=(4,5))*.21
        B=X@X.T
        T=rng.normal(size=(4,4));T=(T+T.T)/2
        S=float(np.trace(ginv@B));Q=float(np.trace(ginv@B@ginv@B))
        L4=ALPHA*S*S+BETA*Q
        V=-metric*L4+4*ALPHA*S*B+4*BETA*B@ginv@B
        volume=1/np.sqrt(np.linalg.det(ginv))
        analytic=float(volume*np.sum(V*T)/2)
        step=2e-5
        finite=(four_gradient_density(ginv+step*T,B)-four_gradient_density(ginv-step*T,B))/(2*step)
        assert abs(finite-analytic)<2e-9
        trace=float(np.trace(ginv@V))
        assert abs(trace)<1e-15
        rows.append(dict(case=index,analytic=analytic,finite_difference=float(finite),
                         error=abs(float(finite)-analytic),four_dimensional_trace=trace))
    return dict(rows=rows,definition='V_mn=2/sqrt(g) delta I4/delta g^mn; Euclidean')


def curvature_gradient_variation_check():
    # Exact 4D conformal geometry, with all profiles depending only on x.
    x=np.linspace(0,2*np.pi,8192,endpoint=False)
    amplitude=.3;n=2
    B=amplitude**2*n**2*np.cos(n*x)**2
    Bp=-amplitude**2*n**3*np.sin(2*n*x)
    Bpp=-2*amplitude**2*n**4*np.cos(2*n*x)
    sigma=.12*np.cos(x)+.04*np.sin(2*x)
    sp=-.12*np.sin(x)+.08*np.cos(2*x)
    spp=-.12*np.cos(x)-.16*np.sin(2*x)
    f=.37*np.cos(2*x)+.16*np.sin(x)+.21*np.cos(4*x)
    fp=-.74*np.sin(2*x)+.16*np.cos(x)-.84*np.sin(4*x)
    fpp=-1.48*np.cos(2*x)-.16*np.sin(x)-3.36*np.cos(4*x)
    integral=lambda values:float(2*np.pi*np.mean(values))
    # sqrt(g) R S = -6 B (sigma''+sigma'^2); no discretized curvature.
    def action(t):return integral(-6*B*(spp+t*fpp+(sp+t*fp)**2))
    direct=integral(-6*B*(fpp+2*sp*fp))
    finite=(action(1e-4)-action(-1e-4))/(2e-4)
    # W trace=6 Box S only when the metric variation of S is retained.
    weighted_box=Bpp-2*spp*B-2*sp*Bp
    via_variation=integral(-6*f*weighted_box)
    scalar_R=-6*np.exp(-2*sigma)*(spp+sp*sp)
    S=np.exp(-2*sigma)*B
    without_RB=via_variation+integral(2*np.exp(4*sigma)*f*scalar_R*S)
    assert abs(finite-direct)<1e-10 and abs(via_variation-direct)<1e-12
    assert abs(without_RB-direct)>1e-3
    # R=0 does not imply vanishing first variation of int R S.
    flat_direction=np.cos(2*n*x)
    flat_variation=integral(-6*flat_direction*Bpp)
    exact_flat=12*np.pi*amplitude**2*n**4
    assert abs(flat_variation-exact_flat)<1e-12
    return dict(conformal_metric_variation=dict(direct=direct,finite_difference=finite,
        from_tensor_trace=via_variation,error=abs(finite-direct),
        incorrect_value_if_metric_dependence_of_S_is_omitted=without_RB),
        flat_background=dict(action_value=0.,nonzero_variation=flat_variation,
            exact_variation=float(exact_flat),heat_kernel_mixed_coefficient=-1/9),
        not_a_Lorentz_pressure_or_solved_Einstein_system=True)


def run():
    candidate=entry.run()
    assert candidate==json.loads((HERE/'round580_drafts/scalar_effective_action_entry_results.json').read_text('utf8'))
    assert candidate['checks_passed']==4
    evidence={k:candidate[k] for k in ('inherited_potential','covariant_trace','derivative_structures','frozen_frame_check')}
    evidence['periodic_nonclosure']=periodic_nonclosure_check()
    evidence['four_gradient_variation']=four_gradient_variation_check()
    evidence['curvature_gradient_variation']=curvature_gradient_variation_check()
    deps=('joint_curved_quantum_source.py','research_note_553.md','research_note_574.md',
          'research_note_579.md','research_round_579_checks.json',
          'round580_drafts/scalar_effective_action_entry_probe.py',
          'round580_drafts/scalar_effective_action_entry_results.json')
    return dict(round=580,tests_run=len(evidence),failures=0,errors=0,checks=list(evidence),
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='fixed four-dimensional Euclidean Einstein metric, original H5 target and U, scalar determinant only, zero gauge curvature, local one-loop terms and their metric variation; off-shell closure obstruction in fixed field variables, not full-model net divergence, independent physical operator count, Lorentz observables, or graph-continuum renormalization')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
