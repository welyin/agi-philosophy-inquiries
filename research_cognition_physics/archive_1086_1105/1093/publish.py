"""One-time root publication. Reuses the single reversible navigation layer."""
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


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def command_json(path,*args):
    proc=subprocess.run([sys.executable,'-B','-X','utf8',str(path),*args],
        capture_output=True,text=True,encoding='utf8')
    assert proc.returncode==0, proc.stdout+proc.stderr
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
        assert review['status']=='PASS_CONDITIONAL_BRIDGE'
        assert review['source_sha256']==reviewed
    for path,value in reviewed.items():
        assert sha((STAGE/path).read_bytes())==value, path
    check=command_json(HERE/'check.py')
    assert check['status']=='PASS' and check['groups']==9

    layerpath=BASE/'_migration/active_documents/layer.json'
    original_layer=layerpath.read_bytes()
    layer=json.loads(original_layer)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验。F＋U＋C＋P及CO1—CO5共同保留，G未采用，目标内容未改。

**最新正式1093，累计3860。** [本轮报告](archive_1086_/research_note_1093.md)补上条件桥：同一实际钟尺下，单基线最快传输时间非退化，加完整传输的齐性复制、真实串接和统一有限尺度承接，推出有限速度包络。承接合同尚无认知来源证明；若另有固定正速理想接力，速度界也可反推该合同，须防止循环。新公理0，科学计0。[证明](archive_1086_/1093/proof.md) · [结果](archive_1086_/1093/results.json) · [验收](archive_1086_/1093/acceptance.json)。

反例菜单具有非退化距离、因果次序和普通局部信号，仍容许无界的长程速度优势；它仅区分任务权限，不是全部九项基础的实体反例。速度随观察者变化不自动颠倒因果顺序。1093未推出运动参考不变性、所有弱影响的前沿或Lorentz。[下一步](archive_1086_/1093/NEXT.md)优先审计承接的认知来源与实际运动运输，复用1041、1086、1089。

[1092](archive_1086_/research_note_1092.md)的有限能力与真实扩容结论继续保留；[用户测距猜想准入](archive_1086_/_admission/nondegenerate_distance_and_order.md)为1093之前的分析记录。用户最新补充已记为[双方共识过程的共同先后](archive_1086_/_admission/bilateral_consensus_event_order.md)：必须检验实际发起、接收、回应和确认的完整运输，复用1041，未另立1094。研究资料在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[新目标回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091结项属于旧目标的限定反证，不代表新猜想已判定。

'''
    stagebody='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1093](research_note_1093.md)、累计3860，新公理0。非退化单基线加实际复制串接及统一有限尺度承接，给有限完成速度包络；尚未给承接的认知来源或运动参考不变性。

[1093证明](1093/proof.md) · [结果](1093/results.json) · [验收](1093/acceptance.json) · [下一步](1093/NEXT.md)。两代理只读独立审阅与9组复算通过，工具审阅消息由主线如实记录。权限反例未冒充自治九项联合模型，精确和有限误差范围分别列明。

[1092](research_note_1092.md)的有界类型稳定与真实扩容限制仍有效；[非退化测距准入](_admission/nondegenerate_distance_and_order.md)为此前记录。[新猜想](../../猜想/可扩展共识与共同洛伦兹结构猜想.md)整体未判定，原四原则、五公理及独立G状态均不变。用户强调的[双方共识形成过程及共同先后](_admission/bilateral_consensus_event_order.md)已记入接续，未另立1094。下一步审计真实接力的认知来源及运动主体间完整任务运输，不继续扩建反例候选。

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
- [1093报告](research_note_1093.md)：非退化测距经统一有限尺度承接给有限速度包络；未推出运动参考不变性。
- [1093证明](1093/proof.md)、[代码](1093/check.py)、[结果](1093/results.json)、[验收](1093/acceptance.json)、[只读审阅记录](1093/review_evidence.json)、[核验](1093/verify.py)与[下一步](1093/NEXT.md)。
- [双方共识过程的共同先后](_admission/bilateral_consensus_event_order.md)：用户最新重心澄清，复用1041；保同一协议依赖仍不等于速度不变。
'''.replace('\n',nl)
            pos=current.index('\n')+1
            new=current[:pos]+addition+current[pos:]
            new=new.replace('用户最新猜想的准入与旧结果复用，未占用1093。',
                '用户猜想的准入与旧结果复用，后续正式条件桥见1093。')
        if new==current: continue
        baseline=baseline_bytes.decode(enc)
        a,b=baseline.splitlines(keepends=True),new.splitlines(keepends=True)
        offsets=[0]
        for line in a: offsets.append(offsets[-1]+len(line))
        edits=[]
        for tag,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
            if tag!='equal':
                edits.append(dict(start=offsets[i],end=offsets[j],before=''.join(a[i:j]),
                    after=''.join(b[k:l]),kind='round1093_nondegenerate_refinement_bridge'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==baseline_bytes
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(reviewed)
    for name in ['verify.py','publish.py','review_evidence.json']:
        assets['1093/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['1086/proof.md','1089/proof.md','1092/proof.md','1092/acceptance.json']:
        p=STAGE/name
        prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1093,date='2026-10-10',status='PASS_CONDITIONAL_BRIDGE',
        entire_conjecture_decided=False,Lorentz_derived=False,
        joint_FUCP_CO_counterexample_certified=False,universal_weak_influence_front_derived=False,
        new_adopted_axioms=0,scientific_count_increment=0,scientific_count_total=3860,
        assets_sha256=assets,prior_sources_sha256=prior,independent_read_only_reviews=2,
        review_delivery='Root records actual read-only agent messages in review_evidence.json; no reviewer-authored files.',
        known_remaining=['cognitive source of uniform finite-scale handoff',
            'common actual clock and spatial metrology', 'actual motion transport and unit reciprocity',
            'precision-time uniformity if approximate tasks are extended to the whole family'],goal_status='active')
    assert layerpath.read_bytes()==original_layer
    for p,old,raw in changes: assert p.read_bytes()==old, str(p)
    acceptance.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for p,old,raw in changes: p.write_bytes(raw)
    layerpath.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    verification=command_json(HERE/'verify.py')
    layout=command_json(ROOT/'scripts/run_research_active.py','--verify-layout')
    previous=command_json(STAGE/'1092/verify.py')
    receipt=dict(date='2026-10-10',round=1093,status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,
        previous_round_verification=previous,
        actual_new_review_groups=9,new_scientific_count=0,
        navigation_documents_updated=6,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json'
    maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',
        scope='Round1093 single-baseline nondegeneracy and actual finite-scale handoff bridge',
        active_documents=7,layer_sha256=sha(layerpath.read_bytes()),layout_verified=True,
        goal_status='active',latest_formal_round=1093,new_formal_rounds=1,
        old_python_json_unchanged=6949,formula_bodies_verified=8340,
        receipt='../../archive_1086_/1093/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1093,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=9,layout_passed=True,
        previous_round_passed=True,goal_status='active'),ensure_ascii=False))


if __name__=='__main__': main()
