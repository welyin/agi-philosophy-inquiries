"""Read-only asset, formula and local-link verification for accepted round 1066."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
from urllib.parse import unquote

HERE=Path(__file__).resolve().parent

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def main():
    receipt=json.loads((HERE/'mainline_acceptance.json').read_text(encoding='utf8'))
    assert receipt['accepted'] and not receipt['whole_roadmap_completed']
    assets=0
    for cat in ('author_assets_sha256','review_assets_sha256','baseline_assets_sha256'):
        for rel,expected in receipt[cat].items():
            assert hashlib.sha256((HERE/rel).read_bytes()).hexdigest()==expected,rel
            assets+=1
    author=module('author1066',HERE/'check.py')
    author.compare(author.run(),json.loads((HERE/'results.json').read_text(encoding='utf8')))
    independent=module('independent1066',HERE/'independent_check.py')
    assert independent.run()==json.loads((HERE/'independent_results.json').read_text(encoding='utf8'))
    links=0
    for rel in ('../research_note_1066.md','plan.md','proof.md','source_scope.md',
                'dependency_update.md','independent_review.md','mathematical_review.md'):
        p=HERE/rel
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf8')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'):continue
            target=unquote(target.split('#')[0].strip('<>'))
            assert (p.parent/target).exists(),(rel,target)
            links+=1
    out={'round':1066,'passed':True,'frozen_assets':assets,'local_links':links,
         'author_recomputed':True,'independent_recomputed':True,'whole_roadmap_completed':False}
    print(json.dumps(out,ensure_ascii=False));return out

if __name__=='__main__':main()
