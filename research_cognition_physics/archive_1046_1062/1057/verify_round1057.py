"""Author freeze/replay; never counts an absent independent review as passed."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import composition_stability as science
import record_bound_check as records

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1057_checks.json'
OWN=['../research_note_1057.md','proof.md','sources.md','dependency_update.md','NEXT.md',
     'composition_stability.py','results.json','record_bound_check.py','record_bound_results.json',
     'verify_round1057.py']
HISTORY=['archive_231_258/research_note_231.md','archive_231_258/research_note_234.md',
         'archive_231_258/research_note_242.md','archive_370_428/research_note_392.md',
         'archive_585_628/research_note_625.md','archive_702_741/research_note_723.md',
         'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
         'archive_1009_1043/research_note_1041.md',
         'archive_1046_/_admission/scalable_organization_after1055/selection.md',
         'archive_1046_/1056/mainline_acceptance.json',
         'archive_1046_/_shared/empirical_interfaces_after1056_checks.json']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(freeze):
    for module,filename in [(science,'results.json'),(records,'record_bound_results.json')]:
        actual=module.calculate()
        science.compare(json.loads((HERE/filename).read_text(encoding='utf-8')),actual)
        assert actual['all_checks_passed']
    links=0
    for name in OWN:
        path=HERE/name
        assert path.is_file(),name
        if path.suffix!='.md': continue
        body=path.read_text(encoding='utf-8')
        body=re.sub(r'\\\[[\s\S]*?\\\]','',body)
        body=re.sub(r'\$\$[\s\S]*?\$\$','',body)
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',body):
            target=target.strip().strip('<>')
            if target.startswith(('http://','https://','#')): continue
            destination=(path.parent/target.split('#',1)[0]).resolve()
            assert destination.is_file(),(name,target)
            links+=1
    record={'round':1057,'date':'2026-10-09','author_checks_passed':True,
            'owned_sha256':{name:digest(HERE/name) for name in OWN},
            'historical_sha256':{name:digest(ROOT/name) for name in HISTORY},
            'local_links_checked':links,'independent_final_review_complete_at_author_freeze':False,
            'accepted_round_at_author_freeze':1056,'accepted_cumulative_at_author_freeze':3832,
            'new_cognitive_axioms':0,'whole_roadmap_complete':False,
            'scope':'finite declared process composition contract; no cognitive derivation of joint permissions'}
    if freeze:
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(record,out,ensure_ascii=False,indent=2)
            out.write('\n')
    else:
        assert record==json.loads(RECEIPT.read_text(encoding='utf-8'))
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    r=run(args.freeze)
    print(json.dumps({'round':1057,'author_checks_passed':True,'assets':len(OWN),
                      'history':len(HISTORY),'links':r['local_links_checked'],
                      'independent_review':'pending','accepted_round':1056,
                      'mode':'exclusive_freeze' if args.freeze else 'read_only'}))
