"""751: source retention and exact internal-completion audit on original modes."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round751_drafts'))
import record_source_sufficiency as retained
import erased_source_balance as erased
TARGET=HERE/'joint_record_source_retention_results.json'

def run():
    a=retained.run();b=erased.run()
    assert a==json.loads((HERE/'round751_drafts/record_source_sufficiency_results.json').read_text('utf8'))
    assert b==json.loads((HERE/'round751_drafts/erased_source_balance_results.json').read_text('utf8'))
    deps=('research_note_252.md','research_note_350.md','research_note_593.md','research_note_594.md',
          'research_note_634.md','research_note_649.md','research_note_699.md','research_note_708.md',
          'research_note_730.md','research_note_732.md','research_note_750.md',
          'joint_dynamic_continuum_reference.py','joint_relative_source_development.py',
          'round750_drafts/marked_source_entry.py','round751_drafts/record_source_sufficiency.py',
          'round751_drafts/erased_source_balance.py')
    return dict(round=751,tests_run=2,failures=0,errors=0,
        checks=['original_gauge_even_same_record_states_and_24_retained_sources',
                'same_full_cq_output_opposite_original_source_cost_and_dilation_audit'],
        retained_source_audit=a,erased_source_audit=b,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original128 finite local-symbol audit with unknown even sterile inputs: one common record is not sufficient for all sources; identical complete cq outputs can have opposite original Majorana energy costs. A familiar WAY obstruction applies to the exact original Luders task with independent fixed apparatus, conserved additive endpoint energy and finite expectation domains. No universal GR no-go, no original fixed-reference contradiction, no continuum quantum-gravity state or detector realization. Cognitive working hypotheses and alternate joint architectures are distinguished from proofs.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=751,tests_run=2,failures=0,errors=0)))

