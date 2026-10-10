"""Read-only verification of round1091 joint-premise, scoped non-implication."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    a=json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert a['status']=='PASS'
    assert a['joint_FUCP_CO_qualification'] is True
    assert a['scoped_nonimplication_certified'] is True
    assert a['Lorentz_derived'] is False
    assert a['empirical_exclusion_claimed'] is False
    assert a['new_universal_cognitive_axioms']==0
    for path,h in a['assets_sha256'].items():
        assert sha(STAGE/path)==h,path
    for path,h in a['prior_evidence_sha256'].items():
        assert sha(ROOT/path)==h,path
    reviewed=0
    for stem in ['math_review','scope_review']:
        r=json.loads((HERE/(stem+'_checks.json')).read_text(encoding='utf8'))
        assert r['status'].startswith('PASS'),stem
        assert len(r['source_sha256'])>=8,stem
        for path,h in r['source_sha256'].items():
            assert sha(HERE/path)==h,(stem,path)
            reviewed+=1
        assert sha(HERE/(stem+'.md'))==r['review_sha256']
    links=0
    for path in a['assets_sha256']:
        p=STAGE/path
        if p.suffix!='.md':
            continue
        text=p.read_text(encoding='utf8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in text),(path,'control character')
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http:','https:','#')):
                continue
            assert (p.parent/target.split('#',1)[0]).exists(),(path,target)
            links+=1
        sans_code=re.sub(r'\x60\x60\x60[^\n]*\n.*?\x60\x60\x60','',text,flags=re.S)
        assert not re.search(r'\\[\[\]]',sans_code),(path,'legacy display delimiter')
        assert len(re.findall(r'(?m)^\$\$\s*$',sans_code))%2==0,(path,'unpaired display')
    recomputed=[]
    if '--recompute' in sys.argv:
        for p in [HERE/'check.py',HERE/'finite_certificate.py',STAGE/'1090/verify.py']:
            args=[sys.executable,'-B','-X','utf8',str(p)]
            if p.name=='verify.py':
                args+=['--recompute']
            r=subprocess.run(args,capture_output=True,text=True,encoding='utf8')
            assert r.returncode==0,r.stdout+r.stderr
            out=json.loads(r.stdout)
            assert out.get('status')=='PASS' or out.get('passed') is True,p
            recomputed.append({'file':str(p.relative_to(STAGE)),'passed':True})
        for path,h in a['assets_sha256'].items():
            assert sha(STAGE/path)==h,('recompute wrote asset',path)
    completion=HERE/'goal_completion.json'
    tool_complete=False
    if completion.exists():
        actual=json.loads(completion.read_text(encoding='utf8'))
        goal=actual['goal']
        assert goal['status']=='complete'
        source=json.loads((STAGE/'_shared/goal_scope_update_20261010.json').read_text(encoding='utf8'))
        assert goal['objective']==source['goal']['goal']['objective']
        tool_complete=True
    print(json.dumps({'round':1091,'passed':True,'assets_verified':len(a['assets_sha256']),
       'prior_evidence_verified':len(a['prior_evidence_sha256']),'review_source_hashes_verified':reviewed,
       'local_links':links,'recomputed':recomputed,'scoped_nonimplication_certified':True,
       'Lorentz_derived':False,'empirical_exclusion_claimed':False,
       'goal_completion_receipt_verified':tool_complete},ensure_ascii=False))

if __name__=='__main__':
    main()
