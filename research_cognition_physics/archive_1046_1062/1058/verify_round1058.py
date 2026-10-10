"""Read-only author reproduction with exclusive first-freeze mode."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import shared_preparation as science

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1058_checks.json'
OWN=['../research_note_1058.md','proof.md','sources.md','dependency_update.md',
     'shared_preparation.py','results.json','verify_round1058.py']
HISTORY=['archive_1046_/1057/proof.md','archive_1046_/1057/composition_stability.py',
         'archive_1046_/1057/research_round_1057_checks.json',
         'archive_231_258/research_note_234.md',
         'archive_1009_1043/research_note_1041.md',
         'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
         'archive_1046_/_admission/scalable_organization_after1055/selection.md']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(freeze):
    actual=science.calculate()
    science.base.compare(json.loads((HERE/'results.json').read_text(encoding='utf-8')),actual)
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
            assert (path.parent/target.split('#',1)[0]).is_file(),(name,target)
            links+=1
    record={'round':1058,'date':'2026-10-09','author_checks_passed':True,
            'owned_sha256':{name:digest(HERE/name) for name in OWN},
            'historical_sha256':{name:digest(ROOT/name) for name in HISTORY},
            'local_links_checked':links,'independent_review_complete_at_author_freeze':False,
            'new_cognitive_axioms':0,'whole_roadmap_complete':False,
            'scope':'shared conditional-iid preparation; finite four-copy grouped CPTP contract',
            'finite_symmetry_does_not_establish_shared_representation':True,
            'quarter_power_bound_not_claimed_optimal':True}
    if freeze:
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(record,out,ensure_ascii=False,indent=2); out.write('\n')
    else:
        assert record==json.loads(RECEIPT.read_text(encoding='utf-8'))
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--freeze',action='store_true')
    args=parser.parse_args()
    r=run(args.freeze)
    print(json.dumps({'round':1058,'author_checks_passed':True,'assets':len(OWN),
                      'history':len(HISTORY),'links':r['local_links_checked'],
                      'independent_review':'pending_at_author_freeze',
                      'mode':'exclusive_freeze' if args.freeze else 'read_only'}))
