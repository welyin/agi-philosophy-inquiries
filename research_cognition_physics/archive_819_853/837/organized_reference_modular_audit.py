"""837 working: the same low-degree sources need not give the same entropy slope.

The reference is faithful on the full ten-mode CAR factor. No finite-factor
entropy is identified with spatial-region or gravitational area entropy.
"""
from pathlib import Path
import argparse,itertools,json,math,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'organized_reference_modular_audit_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code


def entropy(epsilon,r):
    values=[(1-epsilon)*(1+r)/2,(1-epsilon)*(1-r)/2]
    rest=epsilon/1022
    return -sum(p*math.log(p) for p in values)-1022*rest*math.log(rest)


def qubit_relative_to_trace(r):
    return .5*(1+r)*math.log1p(r)+.5*(1-r)*math.log1p(-r)


def run():
    gamma,_,_,_,comp,_,_=code.code_data()
    ids=(0,1,4,5,12,14)
    logical=code.mul((0,0,1),code.product(gamma[i] for i in ids))
    assert code.adj(logical)==logical and code.mul(logical,logical)==code.I
    assert comp(logical)[0]=='logical'
    counts=[]
    for degree in range(5):
        count=0
        for subset in itertools.combinations(range(20),degree):
            word=code.product(gamma[i] for i in subset)
            # Tr(P_c L word)=0 if its compressed word is zero or traceless logical.
            kind,_=comp(code.mul(logical,word));assert kind!='scalar';count+=1
        counts.append(dict(degree=degree,exact_zero_directional_responses=count))
    assert comp(code.mul(logical,logical))==('scalar',0)
    rows=[]
    for epsilon in (.1,.4,511/512):
        for r in (.0,.3,.6):
            h=1e-5
            derivative=(entropy(epsilon,r+h)-entropy(epsilon,r-h))/(2*h)
            expected=-(1-epsilon)*math.atanh(r)
            assert abs(derivative-expected)<1e-9
            finite_change=entropy(epsilon,r)-entropy(epsilon,0)
            relative=(1-epsilon)*qubit_relative_to_trace(r)
            assert abs(finite_change+relative)<3e-15
            quadratic_second=-(1-epsilon)/(1-r*r)
            rows.append(dict(outside_code_weight=epsilon,logical_bias=r,
                             entropy_directional_derivative=expected,
                             numerical_derivative_residual=abs(derivative-expected),
                             entropy_second_derivative=quadratic_second,
                             entropy_change_from_logically_unbiased_reference=finite_change,
                             relative_entropy_to_same_unbiased_organized_reference=relative,
                             all_quadratic_and_quartic_source_derivatives=0,
                             six_Majorana_readout_derivative=1-epsilon,
                             minimum_full_factor_eigenvalue=min(epsilon/1022,(1-epsilon)*(1-r)/2)))
    assert all(z['minimum_full_factor_eigenvalue']>0 for z in rows)
    return dict(round=837,status='working_not_formal',formal_test_groups_added=0,
                all_working_checks_passed=True,logical_witness_Majorana_indices=list(ids),
                witness_Hermitian_phase='i times the ordered six-Majorana product',
                exact_low_degree_responses=counts,faithful_full_factor_reference_rows=rows,
                first_law_low_degree_source_representability_on_all_logical_directions_iff_unbiased=True,
                known_round830_Gaussian_entropy_deficit_not_counted_as_new=True,
                native_instantaneous_quadratic_source_mapping_analytic=True,
                full_spatial_region_modular_mapping_proven=False,
                gravitational_area_or_Einstein_response_computed=False,
                scope='Around a faithful organized reference with biased logical contents, every degree-at-most-four source direction is blind while the entropy first derivative is not. A six-Majorana observable sees the missing logical direction. This tests a candidate reference/source identification, not quantum first-law validity or the Einstein equation.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
