"""Publish 1097 through the existing reversible active-document layer."""
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
        assert r['status']=='PASS_CONDITIONAL_CONE_BRIDGE'
        assert r['source_sha256']==hashes and r['actual_recomputation']['groups']==7
    for name,value in hashes.items():assert sha((STAGE/name).read_bytes())==value,name
    check=run(HERE/'check.py');assert check['status']=='PASS' and check['groups']==7
    lp=BASE/'_migration/active_documents/layer.json'
    original=lp.read_bytes();layer=json.loads(original)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、已采用六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1097，累计3860。** [本轮报告](archive_1086_/research_note_1097.md)把同一完整任务的正有限速度确界接到闭凸锥包络，再在额外的双向穷尽运动运输、同标定菜单及仿射事件字典下得到共形Lorentz形式。所需新增合同未自动采用。[证明](archive_1086_/1097/proof.md) · [7组结果](archive_1086_/1097/results.json) · [验收](archive_1086_/1097/acceptance.json)。新公理0，科学计0。

闭凸锥首先是完整未知量子资料任务的方向包络，不是全部物理影响锥；边界信号、实际单位尺度及全部物质运输未得。正反逐任务成本分别有限，已足以保持全部有限资源并集，无需额外统一资源开销；固定预算与规模一致服务仍须另验。背景协变不能代替同标定菜单的运动对称。

[下一步](archive_1086_/1097/NEXT.md)优先审计完整量子交付与所有可辨影响的连接，计入放大、预备资源、编码解码的实际时间与误差；保留非退化／一致承接和实际运动能力互惠的认知来源主线。不重做锥或群分类，不扩建[1096迁移器件](archive_1086_/research_note_1096.md)。整个加强猜想及共同Lorentz物理仍未判定。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091属于旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1097](research_note_1097.md)、累计3860，新公理0。四原则、五公理、六条操作协议及H1—H3共同保留，G未采用。

[1097证明](1097/proof.md) · [结果](1097/results.json) · [验收](1097/acceptance.json) · [下一步](1097/NEXT.md)。完整任务速度确界经闭凸包络，在额外双向穷尽、同标定及仿射字典下约束共形Lorentz变换。7组复算与两份只读审阅通过；未推得全影响锥、实际边界信号、单位尺度或全部物质规律，新合同未采用。

保留[1096](research_note_1096.md)的同资料迁移及[1093](research_note_1093.md)的有限包络来源。后继审计完整量子任务与可辨影响的转换，按实际时间／误差／资源验收；不继续修候选或重做群分类。整个加强猜想未判定。

'''
    addition='''
- [1097报告](research_note_1097.md)：完整任务的双向可兑现性、全资源与固定预算、闭凸锥包络及条件性共形Lorentz连接。
- [1097证明](1097/proof.md)、[代码](1097/check.py)、[结果](1097/results.json)、[验收](1097/acceptance.json)、[独立审阅](1097/review_evidence.json)、[核验](1097/verify.py)与[下一步](1097/NEXT.md)。
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
                after=''.join(b[k:l]),kind='round1097_complete_task_cone'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==base
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(hashes)
    for name in ['verify.py','publish.py','review_evidence.json']:assets['1097/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['archive_1086_/1093/proof.md','archive_1086_/1096/acceptance.json',
            'archive_1086_/_admission/dynamic_contract_audit.md','archive_1086_/1086/proof.md',
            'archive_742_763/755/drafts/cognitive_joint_candidate_working_report.md',
            'archive_1009_1043/research_note_1041.md']:
        p=BASE/name;prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1097,date='2026-10-10',status='PASS_CONDITIONAL_CONE_BRIDGE',
        entire_conjecture_decided=False,Lorentz_derived=False,conditional_conformal_bridge=True,
        adopted_R=False,all_physical_influences_certified=False,physical_clock_scale_fixed=False,
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
    previous=run(STAGE/'1096/verify.py')
    receipt=dict(round=1097,date='2026-10-10',status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,previous_round_verification=previous,
        navigation_documents_updated=6,actual_new_review_groups=7,new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json';maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',scope='Round1097 complete task cone and conditional motion transport',
        active_documents=7,layer_sha256=sha(lp.read_bytes()),layout_verified=True,goal_status='active',
        latest_formal_round=1097,new_formal_rounds=1,receipt='../../archive_1086_/1097/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1097,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=7,layout_passed=True,previous_round_passed=True,
        goal_status='active'),ensure_ascii=False))


if __name__=='__main__':main()
