"""Publish 1098 through the existing reversible active-document layer."""
from pathlib import Path
from difflib import SequenceMatcher
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
BASE=ROOT/'research_cognition_physics'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):return hashlib.sha256(raw).hexdigest()


def run(path,*args):
    p=subprocess.run([sys.executable,'-B','-X','utf8',str(path),*args],
        capture_output=True,text=True,encoding='utf8')
    assert p.returncode==0,p.stdout+p.stderr
    return json.loads(p.stdout)


def main():
    acceptance=HERE/'acceptance.json'
    assert not acceptance.exists(),'Already published; refusing overwrite.'
    reviews=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(reviews['reviews'])==2
    hashes=reviews['reviews'][0]['source_sha256'];assert len(hashes)==5
    for r in reviews['reviews']:
        assert r['status']=='PASS_CONDITIONAL_INFLUENCE_BRIDGE'
        assert r['source_sha256']==hashes and r['actual_recomputation']['groups']==6
    for name,value in hashes.items():assert sha((STAGE/name).read_bytes())==value,name
    check=run(HERE/'check.py');assert check['status']=='PASS' and check['groups']==6
    lp=BASE/'_migration/active_documents/layer.json'
    original=lp.read_bytes();layer=json.loads(original)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1098，累计3860。** [本轮报告](archive_1086_/research_note_1098.md)把[1097的条件性任务锥](archive_1086_/research_note_1097.md)接到全部可辨影响：额外要求实际影响运输、各参考共同钟序正向及同一非零纯轴运动的任意有限次实际接续。每个锥外影响只需一份有限参考链检验，无需先放大或传态。[证明](archive_1086_/1098/proof.md) · [6组结果](archive_1086_/1098/results.json) · [验收](archive_1086_/1098/acceptance.json)。新公理0，科学计0。

只保依赖身份不保证所有跨地时标正向；有限次数和方向的实际参考只给更宽的范围。全部影响被限制在任务锥内；包络相等另需实际串接及有限编码／读取接续，未把读取当免费。有限证书同时保非零记录差异与负时间余量，没有现实仪器数据。OC、MC及全影响运输仍是待验证合同，未由六协议生成；有限速度来源、单位和全部物质动力学未得。

[下一步](archive_1086_/1098/NEXT.md)回到共同钟序、实际运动接续及1093跨资源非退化的认知来源，复用294／575／963／1041，停止放大、传态和成熟群分类的技术扩建。整个加强猜想及共同Lorentz物理仍未判定。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091是旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1098](research_note_1098.md)、累计3860，新公理0。四原则、五公理、六协议及H1—H3保留，G未采用。

[1098证明](1098/proof.md) · [结果](1098/results.json) · [验收](1098/acceptance.json) · [下一步](1098/NEXT.md)。在1097条件性共形参考基础上，实际影响运输、全参考钟序正向及任意有限运动接续把全部可辨影响约束于同一任务锥；包络相等另计实际编码／读取接续。6组复算与两份只读独立审阅通过；新增合同未采用，无实际仪器数据，未推出整体Lorentz物理。

有限参考范围、有限方向及含误差任务分别保留余量；只保因果依赖不等于所有同步时标正向。后继回到认知来源，不扩建通信器件、钟或矩阵分类。[1097](research_note_1097.md)与[1093](research_note_1093.md)的前提责任保持。

'''
    addition='''
- [1098报告](research_note_1098.md)：实际参考有限接续、依赖与共同钟序、可辨影响包络及有限误差证书。
- [1098证明](1098/proof.md)、[代码](1098/check.py)、[结果](1098/results.json)、[验收](1098/acceptance.json)、[独立审阅](1098/review_evidence.json)、[核验](1098/verify.py)与[下一步](1098/NEXT.md)。
'''
    roots={f'research_cognition_physics/{n}' for n in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']}
    changes=[]
    for entry in layer['entries']:
        name=entry['destination'];p=ROOT/name;old=p.read_bytes()
        base=relocation_original_bytes(old,entry)
        enc='utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        content=old.decode(enc);nl='\r\n' if '\r\n' in content else '\n';new=content
        if name in roots:
            start=content.index('## 当前目标：');end=content.index('## 狭义相对论连接',start)
            new=content[:start]+common.replace('\n',nl)+content[end:]
        elif name.endswith('archive_1086_/README.md'):
            start=content.index('## 当前目标：');end=content.index('## 阶段目标',start)
            new=content[:start]+stage.replace('\n',nl)+content[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            pos=content.index('\n')+1
            new=content[:pos]+addition.replace('\n',nl)+content[pos:]
        if new==content:continue
        a,b=base.decode(enc).splitlines(keepends=True),new.splitlines(keepends=True)
        offsets=[0]
        for line in a:offsets.append(offsets[-1]+len(line))
        edits=[]
        for tag,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
            if tag!='equal':edits.append(dict(start=offsets[i],end=offsets[j],before=''.join(a[i:j]),
                after=''.join(b[k:l]),kind='round1098_actual_observer_influence'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==base
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(hashes)
    for name in ['verify.py','publish.py','review_evidence.json']:assets['1098/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['archive_1086_/1097/acceptance.json','archive_1086_/1097/proof.md',
            'archive_1086_/_admission/dynamic_contract_audit.md',
            'archive_1086_/1089/proof.md','archive_1086_/1086/proof.md',
            'archive_1009_1043/research_note_1041.md',
            'archive_046_077/research_note_61.md',
            'archive_370_428/research_note_396.md']:
        p=BASE/name;prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1098,date='2026-10-10',status='PASS_CONDITIONAL_INFLUENCE_BRIDGE',
        entire_conjecture_decided=False,Lorentz_derived=False,conditional_influence_bridge=True,
        new_OC_or_MC_adopted=False,actual_observer_resources_certified=False,all_physical_influences_certified=False,physical_clock_scale_fixed=False,
        all_matter_motion_transport_certified=False,
        new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860,
        assets_sha256=assets,prior_sources_sha256=prior,independent_read_only_reviews=2,
        review_delivery='Root records actual read-only agent messages; no reviewer-authored files.',
        goal_status='active')
    assert lp.read_bytes()==original
    for p,old,_ in changes:assert p.read_bytes()==old,str(p)
    acceptance.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for p,_,raw in changes:p.write_bytes(raw)
    lp.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    verification=run(HERE/'verify.py')
    layout=run(ROOT/'scripts/run_research_active.py','--verify-layout')
    previous=run(STAGE/'1097/verify.py')
    receipt=dict(round=1098,date='2026-10-10',status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,previous_round_verification=previous,
        navigation_documents_updated=6,actual_new_review_groups=6,new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json';maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',scope='Round1098 actual observer iteration and influence order bridge',
        active_documents=7,layer_sha256=sha(lp.read_bytes()),layout_verified=True,goal_status='active',
        latest_formal_round=1098,new_formal_rounds=1,receipt='../../archive_1086_/1098/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1098,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=6,layout_passed=True,previous_round_passed=True,
        goal_status='active'),ensure_ascii=False))


if __name__=='__main__':main()
