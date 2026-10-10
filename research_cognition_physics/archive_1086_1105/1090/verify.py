"""Read-only verification of the declared round1090 model and frozen assets."""
from pathlib import Path
import hashlib,json,re,subprocess,sys
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def main():
    record=json.loads((HERE/'acceptance.json').read_text(encoding='utf8'))
    assert record['status']=='PASS'
    assert record['exact_FUCP_plus_CO_model_certified'] is True
    assert record['joint_SR_goal_proved_or_disproved'] is False
    review=json.loads((HERE/'independent_review_checks.json').read_text(encoding='utf8'))
    assert review['status']=='PASS_WITH_DECLARED_MODEL_INPUTS'
    links=0
    for relative,expected in record['assets_sha256'].items():
        path=STAGE/relative; raw=path.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==expected,relative
        if path.suffix=='.md':
            for target in re.findall(r'\]\(([^)]+)\)',raw.decode('utf8')):
                if target.startswith(('https:','http:','#')): continue
                assert (path.parent/target.split('#',1)[0]).exists(),(relative,target)
                links+=1
    for relative,expected in review['source_sha256'].items():
        assert hashlib.sha256((HERE/relative).read_bytes()).hexdigest()==expected,relative
    assert hashlib.sha256((HERE/'independent_review.md').read_bytes()).hexdigest()==review['review_sha256']
    if '--recompute' in sys.argv:
        run=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'check.py')],capture_output=True,text=True,encoding='utf8')
        assert run.returncode==0,run.stdout+run.stderr
        assert json.loads(run.stdout)['status']=='PASS'
    print(json.dumps({'round':1090,'passed':True,'assets':len(record['assets_sha256']),'local_links':links,
        'review_sources_verified':len(review['source_sha256']),'recomputed':'--recompute' in sys.argv,
        'exact_FUCP_plus_CO_model_certified_under_declared_M1_M2':True,
        'joint_SR_goal_proved_or_disproved':False,'goal_status':'active'},ensure_ascii=False))

if __name__=='__main__':main()
