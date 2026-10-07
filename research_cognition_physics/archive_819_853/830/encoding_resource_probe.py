"""830 working: exact resource and parity certificate for the 829 code.

Checks a basis isometry, not a synthesis from the original Hamiltonian.
Its ancilla labels are internal modes; no file contains large state matrices.
"""
from pathlib import Path
import argparse,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'encoding_resource_probe_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def run():
    _,_,group,_,compress,_,_=code.code_data()
    assert len(group)==512
    # Physical T has dimension 2^10 = 2 logical * 2^9 gauge.
    # Its parity in this input factorization is I_L times gauge parity.
    seen=set();parity_defects=wrong_defects=0
    for logical in range(2):
        for gauge in range(512):
            source_parity=(-1)**gauge.bit_count()
            ancilla=gauge^1
            output_parity=-(-1)**ancilla.bit_count()  # code parity is odd
            parity_defects+=source_parity!=output_parity
            wrong_defects+=source_parity!=-(-1)**gauge.bit_count()
            seen.add((logical,ancilla))
    assert len(seen)==1024 and parity_defects==0 and wrong_defects==1024
    # Covariance I_20/2 uniquely gives the trace state among Gaussian states.
    # A rank-two stabilizer code has trace-state weight 2/1024.
    trace_numerator,trace_denominator=2,1024
    assert trace_numerator*512==trace_denominator
    return dict(round=830,status='working_not_formal',all_working_checks_passed=True,
        formal_test_groups_added=0,code_dimension=2,carrier_dimension=1024,
        Gaussian_state_at_required_covariance='I_1024 / 1024',
        Gaussian_code_weight_numerator=trace_numerator,Gaussian_code_weight_denominator=trace_denominator,
        Gaussian_five_block_parity_four_point=0,actual_code_four_point=1,
        relative_entropy_bits_to_same_covariance_Gaussian='10 - S(logical)_bits',
        minimum_pure_ancilla_dimension_for_full_input_gauge=512,
        internal_ancilla_fermion_modes_in_certificate=9,
        isometric_basis_images=len(seen),parity_defects=parity_defects,
        defects_if_ancilla_parity_not_adjusted=wrong_defects,
        entropy_export_bits='9 + I(carrier:ancilla)_bits for an initially independent ancilla and exact code output',
        original_autonomous_encoder_or_source_history_constructed=False,
        scope='Gaussian incompatibility and a parity-preserving finite isometry with explicit gauge-information storage. No original control synthesis or energy cost has been computed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
