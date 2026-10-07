"""865 working: translated-loop anchor smearing versus a self-contact.
This is the universal flat equal-time short-distance geometric kernel diagnostic.
It is NOT a full mixed, field-dependent relational quantum calculation and does
not prove that a complete covariant renormalized loop cannot exist.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'loop_self_contact_probe_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def geometric_integral(R,a,N):
    delta=(np.arange(N)+.5)*2*np.pi/N
    # Stable chord length; avoid cancellation in 1-cos(delta).
    return float(2*np.pi*(2*np.pi/N)*np.sum(R*R*np.cos(delta)/(4*R*R*np.sin(delta/2)**2+a*a)))

def run():
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        casimir=sum(T@T for T in old.T)
    error=float(np.max(abs(casimir-(4/3)*np.eye(3))))
    assert error<1e-14
    rows=[]
    for R in (.3,.7):
        for eta in (.2,.1,.05,.025,.0125):
            a=R*eta
            exact=2*np.pi**2*((a*a+2*R*R)/(a*np.sqrt(a*a+4*R*R))-1)
            q1=geometric_integral(R,a,4096);q2=geometric_integral(R,a,8192)
            rel=abs(q2/exact-1);difference=abs(q2-q1)
            assert rel<1e-11 and difference<2e-9
            leading=2*np.pi**2*R/a
            remainder=q2-leading
            rows.append(dict(radius=R,cutoff=a,relative_cutoff=eta,geometric_integral=q2,
                exact_integral=float(exact),relative_quadrature_error=rel,
                quadrature_refinement_difference=difference,
                leading_perimeter_term=float(leading),subtracted_remainder=float(remainder),
                remainder_error_to_flat_limit=float(abs(remainder+2*np.pi**2))))
    for start in (0,5):
        errors=[r['remainder_error_to_flat_limit'] for r in rows[start:start+5]]
        assert all(x>y for x,y in zip(errors,errors[1:]))
        assert rows[start+4]['geometric_integral']>10*rows[start]['geometric_integral']
    return dict(kind='round_865_working_loop_self_contact',formal_reports=864,new_numbered_scientific_groups=0,
        all_checks_passed=True,rows=rows,original_SU3_fundamental_Casimir_error=error,
        fixed_path_same_anchor_integral='translation invariant kernel: normalized common-anchor integral is exactly the same self-contact',
        flat_cutoff_asymptotic='pi*perimeter/a - 2*pi^2 + O(a/R)',
        raw_internal_contraction_is_not_Wick_normal_ordered_loop=True,
        full_original_relational_self_contraction_computed=False,
        cancellation_with_reference_and_gravity_sources_excluded=False,
        quantum_Ward_preserving_path_counterterms_constructed=False,
        full_no_go_claim=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
