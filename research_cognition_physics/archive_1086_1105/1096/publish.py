"""Publish 1096 through the existing reversible active-document layer."""
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
        assert r['status']=='PASS_FINITE_TASK_BRIDGE'
        assert r['source_sha256']==hashes and r['actual_recomputation']['groups']==6
    for name,value in hashes.items():assert sha((STAGE/name).read_bytes())==value,name
    check=run(HERE/'check.py');assert check['status']=='PASS' and check['groups']==6
    lp=BASE/'_migration/active_documents/layer.json'
    original=lp.read_bytes();layer=json.loads(original)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、已采用六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1096，累计3860。** [本轮报告](archive_1086_/research_note_1096.md)将929正在运行的共同事件协议接入1091实际接触：双份匹配器件交换完整内部状态，同一未知资料、执行钟、故障、记录及参考关联在接收端继续演化。补上了直接叠加处理与接收会遇到的非交换接口。[证明](archive_1086_/1096/proof.md) · [6组结果](archive_1086_/1096/results.json) · [验收](archive_1086_/1096/acceptance.json)。新增公理0，科学计0。

这是一份有限六协议任务与FUCP＋CO及Galilean传播的共同有效实现。匹配硬件、完整载荷接触和双份资源是模型输入；三位读者仍是组织内逻辑角色，共同记录仍指内部仪器事件。未认证发射／到达时标的三方共识，也未闭合接触的全部机械反作用，不能签整个加强猜想的反证。

[下一步](archive_1086_/1096/NEXT.md)停止迁移候选技术扩建，检验跨资源完整交付非退化与实际运动角色互惠的认知来源。H3保旧任务不限制新增任务的最快完成时间。复用[1095](archive_1086_/research_note_1095.md)的实际继承机制边界、[1094](archive_1086_/research_note_1094.md)两侧同钟分类和[1093](archive_1086_/research_note_1093.md)有限包络桥；共同Lorentz结构及整个猜想仍开放。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091属于旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1096](research_note_1096.md)、累计3860，新公理0。四原则、五公理、六条操作协议及H1—H3共同保留，G未采用。

[1096证明](1096/proof.md) · [结果](1096/results.json) · [验收](1096/acceptance.json) · [下一步](1096/NEXT.md)。匹配双份有限器件的全状态接触，让929同一运行中资料、程序钟、记录及参考真实迁移，关闭与1091相接的非交换接口。6组复算和两份只读审阅通过。六协议仅按原有限事件／攻击／读窗合同验收；未给独立运动观察者的发射到达共识或全部机械来源。

[1095](research_note_1095.md)、[1094](research_note_1094.md)、[1093](research_note_1093.md)、[1092](research_note_1092.md)及[六协议审计](_admission/six_protocols_and_common_motion.md)保留。停止该迁移技术分支，回到跨资源非退化、实际角色互惠的认知来源；整个加强猜想未判定。

'''
    addition='''
- [1096报告](research_note_1096.md)：匹配双份有限器件，迁移运行中的完整共同协议；同一未知资料与参考保留，原有限六协议合同与实际传播共同实现。
- [1096证明](1096/proof.md)、[代码](1096/check.py)、[结果](1096/results.json)、[验收](1096/acceptance.json)、[独立审阅](1096/review_evidence.json)、[核验](1096/verify.py)与[下一步](1096/NEXT.md)。
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
                after=''.join(b[k:l]),kind='round1096_running_protocol_transfer'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==base
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(hashes)
    for name in ['verify.py','publish.py','review_evidence.json']:assets['1096/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['archive_923_934/research_note_929.md','archive_923_934/929/joint_protocol_selection.py',
            'archive_742_763/755/drafts/cognitive_joint_candidate_working_report.md',
            'archive_1086_/1090/proof.md','archive_1086_/1091/proof.md','archive_1086_/1095/acceptance.json']:
        p=BASE/name;prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1096,date='2026-10-10',status='PASS_FINITE_TASK_BRIDGE',
        entire_conjecture_decided=False,Lorentz_derived=False,finite_task_joint_qualification=True,
        three_independent_motion_observers_implemented=False,
        arrival_and_emission_records_jointly_certified=False,all_mechanical_sources_closed=False,
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
    previous=run(STAGE/'1095/verify.py')
    receipt=dict(round=1096,date='2026-10-10',status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,previous_round_verification=previous,
        navigation_documents_updated=6,actual_new_review_groups=6,new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json';maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',scope='Round1096 matched running protocol state transfer',
        active_documents=7,layer_sha256=sha(lp.read_bytes()),layout_verified=True,goal_status='active',
        latest_formal_round=1096,new_formal_rounds=1,receipt='../../archive_1086_/1096/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1096,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=6,layout_passed=True,previous_round_passed=True,
        goal_status='active'),ensure_ascii=False))


if __name__=='__main__':main()
