"""Publish accepted round1095, retaining the single reversible navigation layer."""
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


def run(path,*args):
    p=subprocess.run([sys.executable,'-B','-X','utf8',str(path),*args],
        capture_output=True,text=True,encoding='utf8')
    assert p.returncode==0,p.stdout+p.stderr
    return json.loads(p.stdout)


def main():
    ap=HERE/'acceptance.json'
    assert not ap.exists(),'Already published; refusing overwrite.'
    reviews=json.loads((HERE/'review_evidence.json').read_text(encoding='utf8'))
    assert len(reviews['reviews'])==2
    hashes=reviews['reviews'][0]['source_sha256']
    assert len(hashes)==5
    for review in reviews['reviews']:
        assert review['status']=='PASS_CONDITIONAL_MECHANISM'
        assert review['source_sha256']==hashes
        assert review['actual_recomputation']['groups']==7
    for name,value in hashes.items(): assert sha((STAGE/name).read_bytes())==value,name
    check=run(HERE/'check.py')
    assert check['status']=='PASS' and check['groups']==7
    lp=BASE/'_migration/active_documents/layer.json'
    original=lp.read_bytes()
    layer=json.loads(original)
    common='''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

目标active，[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)继续检验，目标内容未改。F＋U＋C＋P、CO1—CO5、已采用六条任务限定协议及H1—H3保留，G未采用。[六协议范围](archive_1086_/_admission/six_protocols_and_common_motion.md)与原稿物理推论箭头区分。

**最新正式1095，累计3860。** [本轮报告](archive_1086_/research_note_1095.md)用额外实际分裂原语，把有限预装成员、内部储能和反冲、未知资料继承接到真实发射及完整接收。逐份准备有限、分裂运动账守恒，仍允许随真实增员产生越来越快的源；这份源携带实现与1093的同任务有限包络不相容。[证明](archive_1086_/1095/proof.md) · [7组结果](archive_1086_/1095/results.json) · [验收](archive_1086_/1095/acceptance.json)。新公理0，科学计0。

分裂与Newton运动账是模型输入；本轮未生成接收过程的全部机械反作用，未签为全部认知前提的联合反模型。末次分裂已采用源携带律，机制与包络的矛盾不意味着同样前提存在另一Lorentz分支。固定组成能力上限不等于全部操作总次数有限；当前只限制预装资源支持的分裂次数。

[下一步](archive_1086_/1095/NEXT.md)先审计实际运动能力继承和跨任务交付代价的认知来源，保留比例继承、距离相关处理等替代自由，不继续扩建分裂候选。[1094](archive_1086_/research_note_1094.md)保留同钟左右互惠的条件分类，[1093](archive_1086_/research_note_1093.md)保留有限尺度承接桥，[1092](archive_1086_/research_note_1092.md)保留有界类型与真实扩容边界。共同Lorentz结构及整个加强猜想尚未判定。

研究在[archive_1086_](archive_1086_/README.md)，报告在阶段根目录、资产在编号目录。[目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保持；以下1086—1091属于旧目标的限定反证，不覆盖后来加强的猜想。

'''
    stage='''## 当前目标：可扩展共识与共同洛伦兹结构猜想

目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1095](research_note_1095.md)、累计3860，新公理0。四原则、五公理、六条操作协议及H1—H3共同保留，G未采用。

[1095证明](1095/proof.md) · [结果](1095/results.json) · [验收](1095/acceptance.json) · [下一步](1095/NEXT.md)。额外分裂原语在有限预装成员上保持未知资料及其参考、分裂能量／动量账，生成更快的源；完整交付给与1093包络冲突的有限证书。7组复算与两份只读独立审阅通过，删项对照使用实际时耗。接收机械反作用未生成，没有完整认知反模型或Lorentz推导。

[1094](research_note_1094.md)的两侧同钟分类、[1093](research_note_1093.md)的有限包络桥和[1092](research_note_1092.md)的真实扩容边界保持。[六协议审计](_admission/six_protocols_and_common_motion.md)继续限定采用范围。后继检验角色能力继承及完成代价的来源，不重复弱共识，不建设该候选全部物理学。

'''
    addition='''
- [1095报告](research_note_1095.md)：真实成员、分裂运动账和未知资料继承，经完整交付与有限任务包络比较；额外机制未从认知推导。
- [1095证明](1095/proof.md)、[代码](1095/check.py)、[结果](1095/results.json)、[验收](1095/acceptance.json)、[独立审阅记录](1095/review_evidence.json)、[核验](1095/verify.py)与[下一步](1095/NEXT.md)。
'''
    changes=[]
    root_names={f'research_cognition_physics/{n}' for n in
                ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']}
    for entry in layer['entries']:
        name=entry['destination']
        p=ROOT/name
        old=p.read_bytes()
        base=relocation_original_bytes(old,entry)
        enc='utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        content=old.decode(enc)
        nl='\r\n' if '\r\n' in content else '\n'
        new=content
        if name in root_names:
            start=content.index('## 当前目标：')
            end=content.index('## 狭义相对论连接',start)
            new=content[:start]+common.replace('\n',nl)+content[end:]
        elif name.endswith('archive_1086_/README.md'):
            start=content.index('## 当前目标：')
            end=content.index('## 阶段目标',start)
            new=content[:start]+stage.replace('\n',nl)+content[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            pos=content.index('\n')+1
            new=content[:pos]+addition.replace('\n',nl)+content[pos:]
        if new==content: continue
        a,b=base.decode(enc).splitlines(keepends=True),new.splitlines(keepends=True)
        offsets=[0]
        for line in a: offsets.append(offsets[-1]+len(line))
        edits=[]
        for tag,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
            if tag!='equal':
                edits.append(dict(start=offsets[i],end=offsets[j],before=''.join(a[i:j]),
                    after=''.join(b[k:l]),kind='round1095_finite_actual_role_inheritance'))
        raw=new.encode(enc)
        entry.update(current_sha256=sha(raw),current_bytes=len(raw),edits=edits)
        assert relocation_original_bytes(raw,entry)==base
        changes.append((p,old,raw))
    assert len(changes)==6
    assets=dict(hashes)
    for name in ['verify.py','publish.py','review_evidence.json']:
        assets['1095/'+name]=sha((HERE/name).read_bytes())
    prior={}
    for name in ['1091/proof.md','1093/proof.md','1094/proof.md','1094/acceptance.json']:
        p=STAGE/name
        prior[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    evidence=dict(round=1095,date='2026-10-10',status='PASS_CONDITIONAL_MECHANISM',
        entire_conjecture_decided=False,Lorentz_derived=False,
        all_six_protocols_implemented=False,full_receiver_mechanical_conservation=False,
        new_split_primitive_is_input=True,new_adopted_axioms=0,
        scientific_count_increment=0,scientific_count_total=3860,
        assets_sha256=assets,prior_sources_sha256=prior,independent_read_only_reviews=2,
        review_delivery='Root records actual read-only agent messages; no reviewer-authored files.',
        known_remaining=['cognitive source of motion role permissions and velocity loss budget',
            'same-task all-resource nondegeneracy and handoff',
            'joint qualification of all strengthened hypotheses',
            'actual rods, clocks, matter and complete motion transport'],goal_status='active')
    assert lp.read_bytes()==original
    for p,old,_ in changes: assert p.read_bytes()==old,str(p)
    ap.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for p,_,raw in changes: p.write_bytes(raw)
    lp.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    verification=run(HERE/'verify.py')
    layout=run(ROOT/'scripts/run_research_active.py','--verify-layout')
    previous=run(STAGE/'1094/verify.py')
    receipt=dict(round=1095,date='2026-10-10',status='PASS',goal_status='active',
        goal_content_changed=False,entire_conjecture_decided=False,
        round_verification=verification,layout_verification=layout,
        previous_round_verification=previous,navigation_documents_updated=6,
        actual_new_review_groups=7,new_scientific_count=0,
        stage_paper_sha256=sha((BASE/'认知操作五公理与三维关系空间_阶段论文.md').read_bytes()))
    (HERE/'publication_checks.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    p=BASE/'_migration/active_documents/checks.json'
    maintenance=json.loads(p.read_text(encoding='utf8'))
    maintenance['latest_maintenance']=dict(date='2026-10-10',
        scope='Round1095 finite actual members, split inheritance and complete delivery envelope',
        active_documents=7,layer_sha256=sha(lp.read_bytes()),layout_verified=True,
        goal_status='active',latest_formal_round=1095,new_formal_rounds=1,
        receipt='../../archive_1086_/1095/publication_checks.json')
    p.write_text(json.dumps(maintenance,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(round=1095,status='PASS',navigation_documents_updated=6,
        independent_reviews=2,recomputation_groups=7,layout_passed=True,
        previous_round_passed=True,goal_status='active'),ensure_ascii=False))


if __name__=='__main__': main()
