"""702 entry: nondimensionalize the existing699 certificate, no new integral."""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_flow_resolution_limit as previous
TARGET=HERE/'source_scale_entry_results.json'


def run():
    assert previous.run()==json.loads(previous.TARGET.read_text('utf8'))
    receipt=json.loads((ARCHIVE/'research_round_699_checks.json').read_text('utf8'))
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest()==digest
    old=json.loads((ARCHIVE/'joint_offdiagonal_haar_certificate_results.json').read_text('utf8'))['exact_integral']
    b=F(old['inherited_B00']);u=F(old['inherited_B11_upper']);l=F(old['exact_B10_bounds'][0])
    delta=F(67,10000)
    surplus=l*l-(1+delta)**2*b*u
    assert b>0 and u>0 and l>0 and surplus>0
    deps=('research_note_675.md','research_note_693.md','research_note_699.md','research_note_701.md',
          'joint_offdiagonal_haar_certificate_results.json','research_round_699_checks.json',
          'research_round_701_checks.json')
    return dict(date='2026-10-02',entry_round=702,latest_formal_round=701,new_formal_round=False,
        exact_certificate=dict(B00=str(b),B11_upper=str(u),B10_lower=str(l),
            delta=str(delta),squared_ratio_surplus=str(surplus),
            normalized_ratio_lower_decimal=math.sqrt(float(l*l/(b*u))),
            normalized_negative_upper=str(-delta),
            minimum_uniform_entry_error_for_PSD=str(delta/2)),
        scope=dict(no_new_Haar_or_sphere_integral=True,
            scaling_uses_positive_upper_U_not_unknown_B11=True,
            no_positive_B11_assumption=True,no_division_by_partition_function=True,
            not_an_actual_instrument_or_macroscopic_source_certificate=True,
            no_finite_tau_threshold_or_joint_continuum_claim=True,
            previous_699_and_701_results_preserved=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry=702,formal=701,strict_negative_upper='-67/10000',
                         minimum_PSD_entry_error='67/20000',all_checks_passed=True)))
