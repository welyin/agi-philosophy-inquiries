"""Publish the independently accepted conditional stage; preserve navigation originals."""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import relocation_original_bytes
BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1290_'
HERE=STAGE/'_shared'
LAYER=BASE/'_migration/structure_selection_active_1288/layer.json'
OUT=HERE/'conditional_sm_final_navigation.json'
def sha(raw): return hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists(): raise SystemExit('Final navigation already recorded.')
    acceptance=json.loads((HERE/'conditional_sm_stage_acceptance.json').read_text(encoding='utf-8'))
    assert acceptance['status']=='COMPLETE_WITH_SCOPE'
    for rel,expected in acceptance['files_sha256'].items():
        assert sha((ROOT/rel).read_bytes())==expected,rel
    old=LAYER.read_bytes();layer=json.loads(old)
    for rel,expected in layer['unchanged_support_sha256'].items():
        assert sha((ROOT/rel).read_bytes())==expected,rel
    prefix='''
## 最新结项：标准模型结构的条件重建

2026-10-18。按用户新定的“只需条件性成立”目标完成本阶段。[论文13](13_标准模型结构的条件重建与共同来源_阶段论文.md)、[前件与剩余输入](archive_1290_/_shared/conditional_sm_premise_ledger.md)、[逐项验收](archive_1290_/_shared/conditional_sm_stage_acceptance.json)及独立数学／范围审阅已保存。

结论限于公开字段、群／表示、三代、有效阶和保护条件下的标准模型主要结构，以及SM加三个规范单态右手中微子的Dirac共同分支。规范／物质顶点、质量／混合、反常条件和领先来源相容；完整父作用与1305的有限背景量子过程有明确映射。未将标准群或三代称为新推论，也未把背景证书称为全动态规范／Higgs量子实现。

[1307报告](archive_1290_/research_note_1307.md)新增中性Dirac来源定理：指定领先源类中，Yν(R−kI)=0；满列秩锁共同权，核维数d恰保d²个自由中性方向。[本轮验收](archive_1290_/1307/acceptance.json)与[联合复核](archive_1290_/_shared/conditional_sm_joint_checks.json)通过。1290—1307共18轮，作者／独立计算各18组；本次新增1轮，论文整合计0，新公理0。

[新目标](archive_1290_/_shared/goal_after1306.md)的A—D已全部按公开范围验收。后续若研究输入生成、现实参数或更强动态量子匹配，另定范围；本阶段不自行扩成那些更强目标。下方旧目标、启动和“继续／未结项”状态均属历史，以本节为准。
'''
    hashes={}
    for entry in layer['entries']:
        path=ROOT/entry['destination'];raw=path.read_bytes()
        original=relocation_original_bytes(raw,entry)
        enc='utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        assert len(entry['edits'])==1
        edit=entry['edits'][0]
        marker='\n## 当前长程目标：从认知组合到标准模型结构的共同生成'
        assert marker in edit['after']
        edit['after']=prefix+edit['after'][edit['after'].index(marker):]
        txt=original.decode(enc)
        updated=(txt[:edit['start']]+edit['after']+txt[edit['end']:]).encode(enc)
        entry['current_sha256']=sha(updated);entry['current_bytes']=len(updated)
        assert relocation_original_bytes(updated,entry)==original
        path.write_bytes(updated);hashes[entry['destination']]=sha(updated)
    sp=STAGE/'README.md';txt=sp.read_bytes().decode('utf-8')
    marker='\n- [正式目标、五条主线与整体结项标准]'
    assert marker in txt
    tail=txt[txt.index(marker):]
    row='|1307|Dirac中性部门的领先源充要分类；满秩Yukawa将N接到共同应力，核维数d留d²自由；同父背景映射|[报告](research_note_1307.md) · [证明](1307/proof.md) · [数学独审](1307/independent_review.md) · [范围独审](1307/scope_review.md) · [验收](1307/acceptance.json)|'
    assert '|1307|' not in tail
    lines=tail.splitlines()
    idx=next(i for i,line in enumerate(lines) if line.startswith('|1306|'))
    lines.insert(idx+1,row)
    tail='\n'.join(lines)+'\n'
    tail=tail.replace('当前累计正式验收17轮，作者与独立计算各17组，新采用公理0',
                      '当前累计正式验收18轮，作者与独立计算各18组，新采用公理0',1)
    stage_prefix=prefix.replace('archive_1290_/','').replace('](13_标准模型结构的条件重建与共同来源_阶段论文.md)','](../13_标准模型结构的条件重建与共同来源_阶段论文.md)')
    sp.write_bytes(('# 1290—1307：标准模型条件重建阶段\n'+stage_prefix+'\n## 原长程计划与历史接续\n'+tail).encode('utf-8'))
    hashes[str(sp.relative_to(ROOT)).replace('\\','/')]=sha(sp.read_bytes())
    layer['latest_round']=1307;layer['in_progress_round']=None
    layer['scientific_stage_status']='COMPLETE_WITH_SCOPE'
    layer['current_state']='1290-1307 accepted (18 rounds). User-authorized conditional SM reconstruction A-D complete with independent paper reviews. Full interacting parent structure/source identity plus its declared background quantum certificate; no unique group/generations or full dynamic QFT claim. No further round or goal opened.'
    layer['application_goal_status_note']='Science and paper independently accepted; application goal status to be synchronized after final read-only checks.'
    LAYER.write_bytes((json.dumps(layer,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    result=dict(schema='conditional_sm_final_navigation_v1',status='PASS_WITH_SCOPE',date='2026-10-18',
                latest_round=1307,accepted_rounds=18,author_check_groups=18,independent_check_groups=18,
                new_scientific_rounds=[1307],integration_groups=0,new_axioms=0,
                originals_recoverable=len(layer['entries']),unchanged_support_bindings=len(layer['unchanged_support_sha256']),
                old_layer_sha256=sha(old),new_layer_sha256=sha(LAYER.read_bytes()),files_sha256=hashes,
                acceptance_sha256=sha((HERE/'conditional_sm_stage_acceptance.json').read_bytes()))
    with OUT.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':'PASS','originals_recoverable':len(layer['entries']),'rounds':18}))

if __name__=='__main__': main()
