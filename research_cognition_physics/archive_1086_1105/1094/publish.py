"""One-time 1094 publication through the existing reversible navigation layer."""
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


def sha(raw): return hashlib.sha256(raw).hexdigest()


def command_json(path,*args):
    proc=subprocess.run([sys.executable,'-B','-X','utf8',str(path),*args],
        capture_output=True,text=True,encoding='utf8')
    assert proc.returncode==0,proc.stdout+proc.stderr
    return json.loads(proc.stdout)


def main():
    acceptance=HERE/'acceptance.json'
    if acceptance.exists():
        raise SystemExit('Already published; refusing to overwrite acceptance.')
    reviews=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(reviews['reviews'])==2
    reviewed=reviews['reviews'][0]['source_sha256']
    assert len(reviewed)==5
    for review in reviews['reviews']:
        assert review['status']=='PASS_CONDITIONAL_CLASSIFICATION'
        assert review['source_sha256']==reviewed
    for path,value in reviewed.items():
        assert sha((STAGE/path).read_bytes())==value,path
    check=command_json(HERE/'check.py')
    assert check['status']=='PASS' and check['groups']==7

    layerpath=BASE/'_migration/active_documents/layer.json'
    original_layer=layerpath.read_bytes()
    layer=json.loads(original_layer)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、已采用的六条任务限定协议及H1—H3共同保留，G未采用。[六协议审计](archive_1086_/_admission/six_protocols_and_common_motion.md)恢复准确范围；原稿六条自动通向物理理论的箭头尚非定理。

**最新正式1094，累计3860。** [本轮报告](archive_1086_/research_note_1094.md)检验同一运动钟与左右伙伴、多个信号模式的实际滴答互惠。在公开仿射发射律族中，只余源速度不进入传播的共同c支，或源速度直接带入传播的支；后者仍允许不同模式速度，并通过同一事件的三方描述闭合。实际互惠是公开的加强合同，未从六条自动推出。新公理0，科学计0。[证明](archive_1086_/1094/proof.md) · [7组结果](archive_1086_/1094/results.json) · [验收](archive_1086_/1094/acceptance.json)。

这不是六协议全模型分类或整个猜想的反证，也未推出实际尺、所有物质的共同Lorentz运动学。[下一步](archive_1086_/1094/NEXT.md)审计同一信号接口跨实际发射源的可承接能力，联合复用924权限与1093有限尺度承接；不以重命名预设源独立速度，不扩建候选内部技术。

[1093](archive_1086_/research_note_1093.md)的非退化测距—统一承接—有限包络桥和[1092](archive_1086_/research_note_1092.md)的有限能力／真实扩容结论保留。[第三方共同先后](archive_1086_/_admission/bilateral_consensus_event_order.md)及六协议审计是1094之前的分析记录。研究资料在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[新目标回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091结项属于旧目标的限定反证，不代表当前加强猜想已判定。

'''
    stagebody='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1094](research_note_1094.md)、累计3860，新公理0。四原则、五公理、已采用六条操作协议及H1—H3共同保留，G未采用；[六协议现行范围](_admission/six_protocols_and_common_motion.md)不能被早期未采用状态覆盖。

[1094证明](1094/proof.md) · [结果](1094/results.json) · [验收](1094/acceptance.json) · [下一步](1094/NEXT.md)。同一实际钟的左右多模式互惠排除中间源依赖，却在声明发射律族中留下源携带分支。7组复算含216个实际模型发射／到达事件，两代理只读独立审阅通过；有限误差证书复用963。实际互惠尚非六条推论，候选也未取得全部加强认知要求的共同模型资格。

[1093](research_note_1093.md)保留有限包络条件桥，[1092](research_note_1092.md)保留有限能力与扩容限制。[新猜想](../../猜想/可扩展共识与共同洛伦兹结构猜想.md)整体未判定；后续优先审计实际源更换下的信号接口与完整任务承接，不再重做弱三方共识或替候选建设全部物理学。

'''
    changes=[]
    for entry in layer['entries']:
        name=entry['destination']
        p=ROOT/name
        old=p.read_bytes()
        baseline_bytes=relocation_original_bytes(old,entry)
        enc='utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        current=old.decode(enc)
        new=current
        nl='\r\n' if '\r\n' in current else '\n'
        if name in [f'research_cognition_physics/{n}' for n in
                    ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]:
            start=current.index('## 当前目标：')
            end=current.index('## 狭义相对论连接',start)
            new=current[:start]+common.replace('\n',nl)+current[end:]
        elif name.endswith('archive_1086_/README.md'):
            start=current.index('## 当前目标：')
            end=current.index('## 阶段目标',start)
            new=current[:start]+stagebody.replace('\n',nl)+current[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            addition='''
- [1094报告](research_note_1094.md)：六协议与同一运动钟的实际互惠；公开发射律族中的两支分类及有限误差证书。
- [1094证明](1094/proof.md)、[代码](1094/check.py)、[结果](1094/results.json)、[验收](1094/acceptance.json)、[只读审阅记录](1094/review_evidence.json)、[核验](1094/verify.py)与[下一步](1094/NEXT.md)。
- [六协议现行范围审计](_admission/six_protocols_and_common_motion.md)：已采用操作假说、六条原稿物理箭头与实际运动权限的区别；1094准入前记录。
'''.replace('\n',nl)
            pos=current.index('\n')+1
            new=current[:pos]+addition+current[pos:]
        if new==current: continue
        baseline=baseline_bytes.decode(enc)
        a,b=baseline.splitlines(keepends=True),new.splitlines(keepends=True)
        offsets=[0]
        for line in a: offsets.append(offsets[-1]+len(line))
        edits=[]
        for tag,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
            if tag!='equal':
                edits.append(dict(start=offsets[i],end=offsets[j],before=''.join(a[i:j]),
                    after=''.join(b[k:l]),kind='round1094_same_clock_source_law_classification'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==baseline_bytes
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(reviewed)
    for name in ['verify.py','publish.py','review_evidence.json']:
        assets['1094/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['1086/proof.md','1093/proof.md','1093/acceptance.json',
                 '_admission/six_protocols_and_common_motion.md']:
        p=STAGE/name
        prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1094,date='2026-10-10',status='PASS_CONDITIONAL_CLASSIFICATION',
        entire_conjecture_decided=False,Lorentz_derived=False,
        all_six_protocols_implemented=False,joint_strengthened_counterexample_certified=False,
        new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860,
        assets_sha256=assets,prior_sources_sha256=prior,independent_read_only_reviews=2,
        review_delivery='Root records actual read-only agent messages; no reviewer-authored files.',
        known_remaining=['cognitive source of actual reciprocal response',
            'actual source changes and signal mode identity',
            'one joint process meeting all adopted conditions and resource requirements',
            'actual rods, matter and complete motion transport'],goal_status='active')
    assert layerpath.read_bytes()==original_layer
    for p,old,raw in changes: assert p.read_bytes()==old,str(p)
    acceptance.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for p,old,raw in changes: p.write_bytes(raw)
    layerpath.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    verification=command_json(HERE/'verify.py')
    layout=command_json(ROOT/'scripts/run_research_active.py','--verify-layout')
    previous=command_json(STAGE/'1093/verify.py')
    receipt=dict(date='2026-10-10',round=1094,status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,
        previous_round_verification=previous,actual_new_review_groups=7,new_scientific_count=0,
        navigation_documents_updated=6,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json'
    maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',
        scope='Round1094 same actual clock, opposite partners and source-law classification',
        active_documents=7,layer_sha256=sha(layerpath.read_bytes()),layout_verified=True,
        goal_status='active',latest_formal_round=1094,new_formal_rounds=1,
        old_python_json_unchanged=6949,formula_bodies_verified=8340,
        receipt='../../archive_1086_/1094/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1094,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=7,layout_passed=True,
        previous_round_passed=True,goal_status='active'),ensure_ascii=False))


if __name__=='__main__': main()
