"""1046 author verification. Default read-only; --write first receipt only.

Independent algebraic algorithm is author verification, not another reviewer.
Future external reviews and root navigation are outside the frozen scope.
"""
from pathlib import Path
from fractions import Fraction as F
from urllib.parse import unquote
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1046_checks.json'
OWN=('../research_note_1046.md','proof.md','photon_direction_receiver.py',
     'photon_direction_receiver_results.json','dependency_update.md','NEXT.md','review.md',
     'verify_round1046.py','../_admission/spacetime/selection_audit_1046.md')


def rational_checks(result):
    # Integrate in x=2cos(theta)-1, without using the author's constants/functions.
    def transverse(j):return (F(1,2*j+3)+F(2,2*j+2)+F(5,2*j+1))/32
    def scalar(j):return F(1,4*(2*j+1))
    a1,a2=transverse(1),transverse(2)
    b1,b2=2*(scalar(1)-a1),2*(scalar(2)-a2)
    j2=scalar(2)
    assert a1==F(71,960) and a2==F(31,672)
    assert b1==F(3,160) and b2==F(13,1680) and j2==F(1,20)
    assert 2*a2+b2==2*j2
    alpha2=2*a1*a1/a2
    assert alpha2>F(12,25)**2
    # Angular mean of the circular packet: positive actual momentum direction.
    c_weighted=(F(1,8)+F(3,7)+F(7,6)+F(5,5))/64
    zeta=c_weighted/a2
    assert F(1,2)<zeta<1
    remainder=F(1,7)**3/6/(1-F(1,7)**2/20)
    assert remainder<F(1,2000)
    amplitude=F(1,10)*F(12,25)*F(7,8)-F(1,2000)
    gap=amplitude**2
    cert=result['exact_certificates']
    for key,value in [('a1',a1),('a2',a2),('b1',b1),('b2',b2),('j2',j2),
                      ('alpha_squared',alpha2),('dyson_remainder_upper',remainder),
                      ('amplitude_lower',amplitude),('antipodal_gap_lower',gap)]:
        assert F(cert[key])==value,key
    assert result['calibration_p_t_one']>float(gap)>F(17,10000)
    assert result['new_science_groups']==1 and result['new_cognitive_axioms']==0
    assert not result['roadmap_completed_here']
    assert len(result['direction_instruments'])==8
    assert len(result['unpolarized_counterexamples'])==8
    assert all(row['anisotropic_residual']<2e-12 for row in result['unpolarized_counterexamples'])
    return {'angular_constants_exact':True,'circular_mean_cosine':str(zeta),
            'strict_gap_lower':str(gap),'dyson_remainder_upper':str(remainder),
            'finite_calibration_not_continuum_proof':True}


def links(allow_pending_receipt=False):
    count=0
    for name in OWN:
        path=(HERE/name).resolve()
        if path.suffix!='.md':continue
        for _,url in re.findall(r'\[([^\]\n]*)\]\(([^)\n]+)\)',path.read_text(encoding='utf8')):
            url=url.strip().strip('<>')
            if '://' in url or url.startswith('#'):continue
            target=(path.parent/unquote(url.split('#',1)[0])).resolve()
            if allow_pending_receipt and target==RECEIPT:pass
            else:assert target.exists(),(name,url)
            count+=1
    return count


def verify():
    result=json.loads((HERE/'photon_direction_receiver_results.json').read_text(encoding='utf8'))
    env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1'
    replay=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'photon_direction_receiver.py')],
                          check=True,capture_output=True,text=True,encoding='utf8',env=env)
    report=json.loads(replay.stdout)
    assert report['all_checks_passed'] and report['mode']=='read_only_compare'
    for path,sha in result['historical_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    return {'round':1046,'date':'2026-10-08','all_checks_passed':True,
            'scientific_read_only_replay':True,'independent_person_review_certified_here':False,
            'new_science_groups':1,'additional_science_groups_from_replay':0,
            'roadmap_completed_here':False,'author_second_algorithm':rational_checks(result),
            'historical_sha256':result['historical_sha256'],
            'owned_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
            'local_links_checked':links(allow_pending_receipt=True),
            'freeze_scope':'Nine named author assets; excludes future independent reviews, integration and live navigation.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    checked=verify()
    if args.write:
        with RECEIPT.open('x',encoding='utf8') as f:json.dump(checked,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert checked==json.loads(RECEIPT.read_text(encoding='utf8'))
    print(json.dumps({'round':1046,'all_checks_passed':True,'mode':'exclusive_write' if args.write else 'read_only_compare',
                      'owned_files':len(OWN),'historical_files':len(checked['historical_sha256']),
                      'local_links':links(),'algebra':checked['author_second_algorithm']},ensure_ascii=False))
