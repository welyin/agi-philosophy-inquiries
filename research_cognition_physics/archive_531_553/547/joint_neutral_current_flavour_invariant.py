"""547: flavour-independent neutral-current conditions on the same candidate.

Tree-level closed heavy channel, exact protected mixing; no G_mu substitution,
precision electroweak fit, or identification of open-heavy events as invisible.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_singlet_common_mass_rg as rg
import joint_singlet_spectrum_and_weights as spec

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_neutral_current_flavour_invariant_results.json'


def current_blocks(c,z,rho):
    _,_,M,_=spec.neutral_spectrum(c,z,rho)
    _,U=np.linalg.eigh(M.conj().T@M)
    C=U[:3,:].conj().T@U[:3,:]
    return C


def effective_count(alpha):
    return 2+(1-alpha)**2


def charged_observables(alpha,w):
    """w=sin² theta_W; charged masses neglected, fixed tree couplings."""
    gl=w-.5;gr=w;den=gl*gl+gr*gr
    asym=(gl*gl-gr*gr)/den
    ratio=effective_count(alpha)*.25/den
    return ratio,asym


def run():
    rng=np.random.default_rng(547);checks=[];current_error=0.;trace_rows=[]
    for S,z,rho in ((.3,.4,.7),(.1,.2,1.5),(.6,.1,.2)):
        alpha=S/(S+2*z*rho)
        for _ in range(8):
            c=rng.normal(size=3)+1j*rng.normal(size=3);c*=math.sqrt(S)/np.linalg.norm(c)
            C=current_blocks(c,z,rho)
            LL=float(np.linalg.norm(C[:3,:3],'fro')**2)
            LH=2*float(np.linalg.norm(C[:3,3:],'fro')**2)
            HH=float(np.linalg.norm(C[3:,3:],'fro')**2)
            expected=[effective_count(alpha),2*alpha*(1-alpha),alpha*alpha]
            current_error=max(current_error,float(max(abs(np.array([LL,LH,HH])-expected))))
            assert np.allclose(C@C,C,atol=2e-14)
            assert np.allclose([LL,LH,HH],expected,atol=2e-14)
            assert abs(LL+LH+HH-3)<2e-14
        trace_rows.append(dict(alpha=alpha,light_current=LL,mixed_current=LH,heavy_current=HH,
                               massless_all_channels_control=LL+LH+HH))
    assert current_error<1e-13
    checks.append('complex_flavour_independent_full_current_partition_and_massless_open_channel_control')

    # Exact charge algebra followed by independent asymmetry reconstruction.
    reconstruction_error=0.
    for w in (F(1,10),F(23,100),F(1,4),F(2,5)):
        gl=w-F(1,2);gr=w;x=1-4*w;D=1-4*w+8*w*w
        assert 4*(gl*gl+gr*gr)==D==(1+x*x)/2
        assert (gl*gl-gr*gr)/(gl*gl+gr*gr)==2*x/(1+x*x)
        for alpha in (0.,.3,.8):
            ratio,asym=charged_observables(alpha,float(w))
            inferred=ratio/(1+math.sqrt(max(0.,1-asym*asym)))
            reconstruction_error=max(reconstruction_error,abs(inferred-effective_count(alpha)))
    assert reconstruction_error<2e-14
    # Outside |x|<=1 the opposite algebraic branch is needed.
    wrong_w=.8;ratio,asym=charged_observables(.3,wrong_w)
    wrong=ratio/(1+math.sqrt(1-asym*asym))
    assert abs(wrong-effective_count(.3))>.1
    checks.append('charged_width_and_asymmetry_eliminate_common_normalization_with_branch_counterexample')

    for alpha in (F(0),F(1,100),F(1,3),F(4,5),F(1)):
        N=2+(1-alpha)**2;delta=3-N
        assert delta==2*alpha-alpha*alpha and F(0)<=delta<=1
        assert math.isclose(1-math.sqrt(float(1-delta)),float(alpha),abs_tol=1e-15)
        if alpha not in (0,1):
            r=F(7,3);beta=3*(1-alpha)/(3+r*alpha);ZH=r*beta/(3+r*beta)
            assert N==2+((3+r)*ZH/r)**2
    checks.append('exact_deficit_inverse_and_common_bare_Higgs_current_relation')

    # A diagonal-only width is basis dependent for degenerate massless light modes.
    alpha=.6;C=np.diag([1.,1.,1-alpha]).astype(complex)
    U=np.fft.fft(np.eye(3))/math.sqrt(3);rotated=U.conj().T@C@U
    true=float(np.trace(C@C).real);diag_only=float(np.sum(np.abs(np.diag(rotated))**2))
    assert abs(np.trace(rotated@rotated).real-true)<1e-14
    assert true-diag_only>.1
    checks.append('offdiagonal_light_channels_are_required_for_basis_independent_width')

    saved=json.loads((HERE/'joint_top_boundary_identifiability_results.json').read_text('utf8'))
    gy,gw,_=rg.gauge_squared(-saved['endpoint_log_scale'],rg.XSTAR)
    weak=gy/(gy+gw);assert 0<weak<.5
    rows=[]
    for old in saved['inverse_examples']:
        if not old['vacuum']['strict_double_vev']:continue
        s=old['joint_spectrum'];alpha=s['active_weight_in_heavy_pair']
        rZ=4*s['neutral_heavy_squared_over_h_squared']/(gy+gw)
        assert rZ>1
        N=effective_count(alpha);ratio,asym=charged_observables(alpha,weak)
        extracted=ratio/(1+math.sqrt(1-asym*asym))
        assert abs(N-extracted)<2e-14
        rows.append(dict(synthetic_top_target=old['synthetic_target_q'],q0=old['q0'],
            active_weight=alpha,heavy_over_Z_mass_squared=rZ,heavy_Z_channels_closed=True,
            light_neutral_effective_count=N,deficit_from_three=3-N,
            charged_lepton_width_ratio=ratio,charged_lepton_asymmetry=asym,
            reconstructed_count=extracted))
    assert len(rows)==3 and min(r['deficit_from_three'] for r in rows)>.3
    checks.append('same_unique_top_solutions_have_closed_Z_channels_and_flavour_independent_deficits')

    deps=('research_note_545.md','research_note_546.md','joint_singlet_spectrum_and_weights.py',
          'joint_top_boundary_identifiability_results.json','research_round_546_checks.json')
    return dict(round=547,tests_run=len(checks),failures=0,errors=0,checks=checks,
        current_partition_max_error=current_error,asymmetry_reconstruction_max_error=reconstruction_error,
        current_partition_examples=trace_rows,
        diagonal_only_counterexample=dict(correct_inclusive_count=true,incorrect_diagonal_only_count=diag_only),
        frozen_endpoint_weak_sine_squared=float(weak),unique_boundary_examples=rows,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(tree_level_rank_one_protected_model=True,heavy_channels_must_be_closed_for_invisible_formula=True,
           no_other_invisible_decay_channels_assumed=True,charged_lepton_masses_neglected_at_Z=True,
           overall_Gmu_normalization_not_assumed=True,asymmetry_weak_angle_branch_required=True,
           open_heavy_event_selection_included=False,current_experimental_fit=False,
           all_top_targets_or_the_full_unification_excluded=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','current_partition_max_error')},ensure_ascii=False))
