"""Independent-order finite record derivative diagnostics for922; no error bound."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'mixed_derivative_validation_results.json'
def run():
    q=HERE/'background_readout_variation_results.json';d=json.loads(q.read_text('utf-8'))
    direct=np.array(d['centered_background_derivatives'][-1]['discrete_record_derivative'])
    reverse=[np.array(a['derivative']) for a in d['independent_reverse_order_derivatives']]
    rich=(4*reverse[1]-reverse[0])/3
    gap=rich-direct
    scale=np.max(abs(direct))
    assert np.max(abs(gap))<1e-4*scale
    weights=[np.array(a['weighted_source']['real']) for a in d['centered_background_derivatives']]
    relative=np.max(abs(weights[1]-weights[0]))/np.max(abs(weights[1]))
    assert relative<1e-6
    rows=[d['baseline']]+d['deformed_readouts']
    ward=max(a['Ward_identity_error'] for a in rows)
    assert ward<1e-16
    assert d['baseline915_reproduction_error']==0.
    assert len(rows)==5 and all(a['weight']['minimum_sample_calibration_density']>0 for a in rows)
    return dict(round=922,date='2026-10-06',all_diagnostic_checks_passed=True,
        reverse_order_second_order_extrapolation=rich.tolist(),
        extrapolated_reverse_minus_direct=gap.tolist(),extrapolated_order_difference_max=float(np.max(abs(gap))),
        extrapolated_order_relative_difference=float(np.max(abs(gap))/scale),
        raw_reverse_order_step_difference=(reverse[1]-reverse[0]).tolist(),
        weighted_half_step_relative_difference=float(relative),maximum_Ward_array_identity_error=ward,
        mixed_derivative_comparison_uses_independent_differentiation_order=True,
        Richardson_extrapolation_is_numerical_diagnostic_not_error_certificate=True,
        all_material_points_source_kernels_and_original_reference_recomputed=True,
        all_sampled_original_receiver_denominators_positive=True,
        continuum_or_full_domain_bound_certified=False,full_physical_error_certified=False,
        source_hashes={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),q)})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))
