"""Reproduce a bounded public-data compatibility test, not a new SM prediction.

Default mode is read-only. --write exclusively creates the result file.
"""
from pathlib import Path
from statistics import NormalDist
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / 'weak_two_scale_results.json'
HISTORY = (
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_531_553/research_note_547.md',
    'archive_1009_1043/research_note_1014.md',
    'archive_1009_1043/research_note_1016.md',
    'archive_1009_1043/research_note_1043.md',
)


def compute():
    data = json.loads((HERE / 'weak_public_inputs.json').read_text(encoding='utf8'))
    b, a = data['belle_ii'], data['atlas']
    x = b['published_effective_amplitude_ratio']
    sx = b['published_effective_amplitude_ratio_sigma']
    y = a['published_R_W']
    sy = math.sqrt(sum(a[k]**2 for k in
        ('published_R_W_stat', 'published_R_W_syst', 'published_R_W_external')))
    # This reproduces a published conversion from rounded inputs, not extra data.
    wz, z, sz = a['R_WZ'], a['external_R_Z'], a['external_R_Z_sigma']
    reconstructed_w = wz * math.sqrt(z)
    ext_error = wz * sz / (2 * math.sqrt(z))
    assert abs(reconstructed_w-y) < 1e-4
    assert abs(ext_error-a['published_R_W_external']) < 1e-5
    assert abs(sy-a['published_R_W_sigma_rounded']) < 1e-4
    reconstructed_tau = math.sqrt(b['R_tau']/b['published_SM_R_tau_rounded'])
    assert abs(reconstructed_tau-x) < 6e-5

    # Correlation-agnostic joint coverage, assuming each Gaussian marginal.
    q = NormalDist().inv_cdf(1-0.05/(2*2))
    xi, yi = [x-q*sx,x+q*sx], [y-q*sy,y+q*sy]
    assert xi[0] > 0 and yi[0] > 0
    eta = x/math.sqrt(y)-1
    eta_interval = [xi[0]/math.sqrt(yi[1])-1,
                    xi[1]/math.sqrt(yi[0])-1]
    common_interval = [max(xi[0],math.sqrt(yi[0])),
                       min(xi[1],math.sqrt(yi[1]))]
    fixed_in_rectangle = xi[0] <= 1 <= xi[1] and yi[0] <= 1 <= yi[1]
    assert fixed_in_rectangle and common_interval[0] < common_interval[1]
    assert eta_interval[0] < 0 < eta_interval[1]

    def chi(c):
        return ((c-x)/sx)**2+((c*c-y)/sy)**2
    # All stationary points plus boundary; quartic grows at infinity.
    roots = np.roots([2*sx*sx,0,sy*sy-2*sx*sx*y,-sy*sy*x])
    candidates = [0.0]+[float(r.real) for r in roots
                          if abs(r.imag)<1e-10 and r.real>0]
    best = min(candidates,key=chi)
    chi0, chifit = chi(1.0),chi(best)
    eta_sigma = math.sqrt((sx/math.sqrt(y))**2+(x*sy/(2*y**1.5))**2)
    derivative = 2*(best-x)/sx**2+4*best*(best**2-y)/sy**2
    assert abs(derivative) < 1e-6 and chifit <= chi0
    # Saturated two-parameter reduced amplitude model: no remaining GOF.
    sat_a, sat_eta = math.sqrt(y), eta
    assert abs(sat_a*(1+sat_eta)-x) < 1e-14
    assert abs(sat_a**2-y) < 1e-14
    # Without external Z normalization, any positive W value has a matching Z.
    z_for_w = {str(w): (w/wz)**2 for w in (.98,1.0,1.02)}
    assert all(abs(wz*math.sqrt(z_for_w[str(w)])-w)<1e-14
               for w in (.98,1.0,1.02))
    return {
        'round':1052, 'date':'2026-10-08', 'all_checks_passed':True,
        'new_project_empirical_groups':1, 'new_theorems_claimed':0,
        'new_cognitive_axioms':0, 'cognitive_ontology_empirically_selected':False,
        'roadmap_complete':False,
        'historical_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HISTORY},
        'source_conversions':{
            'R_W_from_rounded_WZ_and_Z':reconstructed_w,
            'external_Z_error_propagated':ext_error,
            'W_component_quadrature_sigma':sy,
            'tau_amplitude_from_rounded_SM_normalization':reconstructed_tau,
            'conversions_not_additional_independent_measurements':True},
        'fixed_adopted_SM_branch':{
            'tau_marginal_standardized_residual':(x-1)/sx,
            'W_marginal_standardized_residual':(y-1)/sy,
            'inside_primary_joint_rectangle':fixed_in_rectangle},
        'primary_marginal_gaussian_bonferroni':{
            'nominal_joint_coverage_lower_bound':0.95,
            'requires_independence':False, 'z_per_two_sided_interval':q,
            'tau_amplitude_interval':xi, 'W_ratio_interval':yi,
            'common_amplitude_allowed_interval':common_interval,
            'relative_transport_residual_point':eta,
            'relative_transport_residual_interval':eta_interval,
            'not_distribution_free':True},
        'secondary_independent_gaussian_diagnostics':{
            'fixed_SM_chi_square_two_observables':chi0,
            'fixed_SM_nominal_p_two_dof':math.exp(-chi0/2),
            'common_amplitude_best_fit':best,
            'common_amplitude_min_chi_square':chifit,
            'common_amplitude_nominal_p_one_dof':math.erfc(math.sqrt(chifit/2)),
            'transport_delta_method_sigma':eta_sigma,
            'transport_delta_method_standardized_residual':eta/eta_sigma,
            'not_collaboration_joint_likelihood':True},
        'freedom_diagnostics':{
            'saturated_amplitude':sat_a,'saturated_transport_residual':sat_eta,
            'saturated_model_remaining_goodness_of_fit_dof':0,
            'external_Z_values_reproducing_same_WZ':z_for_w,
            'unique_SMEFT_coefficient_identified':False}
    }


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(u,v) in enumerate(zip(a,b)):compare(u,v,path+f'[{i}]')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=compute()
    if args.write:
        with RESULT.open('x',encoding='utf8') as f:json.dump(out,f,indent=2,ensure_ascii=False);f.write('\n')
    else:compare(out,json.loads(RESULT.read_text(encoding='utf8')))
    print(json.dumps(out,ensure_ascii=False))


if __name__=='__main__':main()
