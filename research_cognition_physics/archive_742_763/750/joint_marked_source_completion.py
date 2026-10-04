"""750: actual marked sources and joint probability-completion condition."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round750_drafts'))
import marked_source_entry as marked
import branch_probability_audit as audit
TARGET=HERE/'joint_marked_source_completion_results.json'
def run():
    a=marked.run();assert a==json.loads(marked.TARGET.read_text('utf8'))
    b=audit.run();assert b==json.loads(audit.TARGET.read_text('utf8'))
    deps=('research_note_634.md','research_note_730.md','research_note_731.md',
          'research_note_732.md','research_note_735.md','research_note_741.md',
          'research_note_742.md','research_note_749.md',
          'round750_drafts/marked_source_entry.py','round750_drafts/branch_probability_audit.py')
    return dict(round=750,tests_run=3,failures=0,errors=0,
        checks=['actual_conditional_CAR_and_24_sources_independent_Fock_check',
                'marked_source_covariance_and_probability_weight_derivative',
                'outcome_specific_past_normalization_audit'],
        marked_sources=a,joint_probability_audit=b,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Conditional covariances of the original occupation instrument determine all bilinear marked sources on the same pure quasifree reference. Existing continuum smoothness/Ward and linear response apply under their stated hypotheses. Separate background-dependent conditional constructions need a joint normalization contract; two explicit original-matrix history diagnostics fail it. No outcome-dependent geometry has been proved to arise from one causal internal detector, no nonlinear semiclassical completion or GR derivation is claimed.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=750,tests_run=3,joint_probability_contract_checked=True)))
