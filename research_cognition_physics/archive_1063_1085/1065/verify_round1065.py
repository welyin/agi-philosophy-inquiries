"""Read-only verification of accepted round 1065 assets and both computations."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
from urllib.parse import unquote

HERE=Path(__file__).resolve().parent

def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def main():
    receipt=json.loads((HERE/'mainline_acceptance.json').read_text(encoding='utf8'))
    assert receipt['accepted'] and not receipt['whole_roadmap_completed']
    count=0
    for category in ('author_assets_sha256','review_assets_sha256','baseline_assets_sha256'):
        for rel, expected in receipt[category].items():
            path=HERE/rel
            assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, str(path)
            count+=1
    author=load('author1065',HERE/'check.py')
    author.compare(author.run(),json.loads((HERE/'results.json').read_text(encoding='utf8')))
    independent=load('independent1065',HERE/'independent_check.py')
    result=independent.run()
    author.compare(result,json.loads((HERE/'independent_results.json').read_text(encoding='utf8')))
    links=0
    for rel in ('../research_note_1065.md','plan.md','proof.md','dependency_update.md',
                'independent_review.md','mathematical_review.md'):
        p=HERE/rel
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf8')):
            if re.match(r'^[a-zA-Z]+://',target) or target.startswith('#'): continue
            target=unquote(target.split('#')[0].strip('<>'))
            assert (p.parent/target).exists(), (str(p),target)
            links+=1
    output={'round':1065,'passed':True,'frozen_assets_checked':count,
            'local_links_checked':links,'author_recomputed':True,
            'independent_recomputed':True,'whole_roadmap_completed':False}
    print(json.dumps(output,ensure_ascii=False))
    return output

if __name__=='__main__': main()
