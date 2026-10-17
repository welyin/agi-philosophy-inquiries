"""Exact finite-tree instantiation of certified native instruments; no new physics simulation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import itertools
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def distribution(weight):
    return {(0,0):F(3,4)*(1-weight),(0,1):F(3,5)*weight,
            (1,0):F(1,4)*(1-weight),(1,1):F(2,5)*weight}

def calculate():
    baseline=json.loads((HERE/'source_baseline.json').read_text(encoding='utf-8'))
    for path,h in baseline['source_sha256'].items():
        assert sha(ROOT/path)==h,path
    stage=ROOT/'research_cognition_physics/archive_1234_'
    assert sha(stage/'_shared/goal_prompt.txt')==baseline['goal_sha256']
    old=json.loads((ROOT/'research_cognition_physics/archive_1219_1233/1219/results.json').read_text(encoding='utf-8-sig'))
    old_rows=old['group3_joint_predictions_and_future']['rows']
    matched=0
    for row in old_rows:
        if row['preparation']!='cold_A_hot_B': continue
        q,l=row['unknown_control_basis'],row['program']
        weight=F(1) if q==1 and l in (0,2) else F(1,2)
        d=distribution(weight)
        assert all(abs(float(d[j,r])-row['joint_J_direction'][j][r])<1e-12
                   for j,r in d)
        matched+=1
    assert matched==6
    outputs=[]
    for q in [F(0),F(1,2),F(1)]:
        first=distribution((1+q)/2)
        p=sum(w for y,w in first.items() if y==(1,1))
        assert p==(1+q)/5
        records={}
        for mode in ['feedback','log_only']:
            tree=[]
            for y1,w1 in first.items():
                e1=int(y1==(1,1))
                selected=e1 if mode=='feedback' else 1
                second=distribution(F(1) if selected==0 else F(1,2))
                for y2,w2 in second.items():
                    tree.append((y1,y2,e1,selected,w1*w2))
            assert len(tree)==16 and sum(x[4] for x in tree)==1
            assert all(x[4]>=0 for x in tree)
            event2=sum(w for y1,y2,e1,l,w in tree if y2==(1,1))
            success=sum(w for y1,y2,e1,l,w in tree if e1 or y2==(1,1))
            negative=sum(w for y1,y2,e1,l,w in tree if y2[1]==1)
            mean_x=F(8,5)*(1-2*negative)
            records[mode]={'second_event':event2,'any_event':success,
                           'negative_port':negative,'mean_x_offset':mean_x}
            assert event2==(F(2,5)-p/5 if mode=='feedback' else F(1,5))
            assert success==(F(2,5)+3*p/5 if mode=='feedback' else F(1,5)+4*p/5)
            # All path records remain present, including both failure outcomes.
            assert sum(w for y1,y2,e1,l,w in tree if not e1 and y2!=(1,1))==1-success
        gap=records['feedback']['second_event']-records['log_only']['second_event']
        assert gap>0 and gap>=F(3,25)
        outputs.append({'q_population':str(q),'first_event':str(p),
                        'programs':{m:{k:str(v) for k,v in d.items()} for m,d in records.items()},
                        'event_gap':str(gap)})
    interval=F(2)**2-F(8,5)**2
    assert interval==F(6,5)**2==F(36,25)
    # Fractional example is a conditional margin, not an actual achieved apparatus error.
    eps=F(1,100); margin=F(3,25)-2*eps
    assert margin==F(1,10)>0
    return {'status':'PASS_OLD_NATIVE_PROTOCOL_INSTANTIATION','matched_old_probability_rows':matched,
            'finite_cases':outputs,'every_branch_Minkowski_interval':str(interval),
            'conditional_error_example':{'epsilon':str(eps),'guaranteed_gap':str(margin)},
            'new_scientific_rounds':0,'new_formal_check_groups':0,'new_adopted_axioms':0,
            'scope':'Exact probability composition only; native quantum process and geometry inherited from the cited proofs.',
            'not_certified':['New native instruments','Window-local stress','Nonflat response',
                             'Empirical error achievement','Full strong governance conjecture']}

def main():
    result=calculate(); target=HERE/'results.json'
    if sys.argv[1:]==['--write']:
        with target.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        mode='created_once'
    else:
        assert not sys.argv[1:]
        assert json.loads(target.read_text(encoding='utf-8'))==result
        mode='read_only_match'
    print(json.dumps({'status':result['status'],'mode':mode,
                      'matched_old_rows':result['matched_old_probability_rows'],
                      'minimum_event_gap':'3/25','new_rounds':0},ensure_ascii=False))

if __name__=='__main__': main()
