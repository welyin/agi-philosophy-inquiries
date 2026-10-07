"""Publish 989; preserve frozen evidence and all previous versions."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
old=(STAGE/'988/verify_round988.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"relation_withdrawal_results.json")',1)[0]
prefix=prefix.replace('range(776,988)','range(776,989)')
prefix=prefix.replace('Delivery checks for 988','Delivery checks for 989')
prefix=prefix.replace('research_round_988_checks.json','research_round_989_checks.json')
checks=r'''    result=read(HERE/'shared_record_quantum_action_results.json')
    assert result['round']==989 and result['all_scientific_checks_passed']
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    import numpy as np
    from fractions import Fraction as F
    def load(name,path):
        spec=importlib.util.spec_from_file_location(name,path)
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
    core=load('core989',HERE/'shared_record_quantum_action.py')
    core.compare(core.run(),result)
    # Independent representation: 36D active CAR sector and 24D local relational code.
    old=load('core956',STAGE/'956/native_material_interface.py')
    mat=old.material(1.,.1)
    cs=[old.annihilate(j) for j in range(4)];ns=[c.T@c for c in cs]
    sector=np.eye(16)[:,[j for j in range(16) if j.bit_count()==2]]
    q=sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    H=np.kron(mat['H'],np.eye(6))+np.kron(np.eye(6),mat['H'])+.2*np.kron(q,q)
    ev,vec=np.linalg.eigh(H)
    e=np.eye(16)
    code=np.stack([(e[5]-e[6]-e[9]+e[10])/2,
        (2*e[3]+2*e[12]-e[5]-e[6]-e[9]-e[10])/math.sqrt(12)],axis=1)
    local=np.kron(mat['W'],np.eye(4))@code
    phase=load('phase965',STAGE/'965/material_field_window.py').phases
    T=result['inputs']['time_from_958']
    targetlocal=local@np.diag([phase(np.array([-mat['J']]),T)[0],1])
    target=np.kron(targetlocal,targetlocal)@np.array([-1.,1.,1.,1.])/2
    initial=np.kron(local,local)@np.ones(4)/2
    def apply(u,state):
        x=state.reshape(6,4,6,4).transpose(0,2,1,3)
        x=(u@x.reshape(36,16)).reshape(6,6,4,4)
        return x.transpose(0,2,1,3).reshape(576)
    residuals=[]
    for row in result['rows']:
        u=(vec*phase(ev,row['time']))@vec.conj().T
        actual=apply(u,initial)
        fid=float(abs(np.vdot(target,actual))**2)
        residuals.append(abs(fid-row['native_fidelity']))
    assert max(residuals)<1e-8
    assert 29*97441**2<524737**2
    analytic=result['analytic']
    assert F(analytic['native_witness_lower']['exact'])==F(48147,50000)
    assert F(analytic['finite_prediction_gap_lower']['exact'])==F(23147,50000)
    assert F(analytic['unread_classical_flag_witness_lower']['exact'])==F(60647,100000)
    assert result['flag_commutator_norm']==0
    assert result['decision']['every_relation_label_must_be_coherent'] is False
    assert result['decision']['gravitational_quantization_proved'] is False
    note=STAGE/'research_note_989.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','不是设备无关Bell实验','不升级为宇宙最大化预测收益',
        '没有证明空间、几何、引力中介均必须量子化','整个区间下界来自解析误差'):
        assert term in prose,term
    newfiles=[note,HERE/'shared_record_quantum_action.py',HERE/'shared_record_quantum_action_results.json',
        Path(__file__),HERE/'drafts/mechanism_adoption_decision.md',HERE/'drafts/mechanism_map_v0_5.md',
        HERE/'drafts/publish989.py',STAGE/'990/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至989）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        match=re.fullmatch(r'research_note_(\d+).md',p.name)
        if match and p.parent.name.startswith('archive_'):nums.append(int(match[1]))
    assert sorted(n for n in nums if n<=989)==list(range(1,990))
    assert '001—989轮共989份' in nav[0].read_text('utf-8-sig')
    assert '231—989的759份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '989：共享证据与不可替代的量子作用' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'988/research_round_988_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=989,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=989,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_988=prev['cumulative_numbered_test_groups_from_987']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_36D_CAR_fidelity_residuals=residuals,
        analytic_rational_window_bounds_verified=True,
        mechanism_decision='classical evidence may select a quantum contact but cannot replace its quantum resources',
        universal_prediction_optimization_adopted=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_inputs','new_scientific_and_entry_files'):assert before[key]==out[key],key
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
'''
target=STAGE/'989/verify_round989.py';assert not target.exists()
target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 989：共享证据与不可替代的量子作用'
    previous='## 988：关系撤销与内部资料保存'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[989报告]({pre}research_note_989.md)复用958全时间界和269资源界：'
        +'无共享纠缠、只靠局部操作与经典消息，不能恢复原材料的四项局部统计，有限窗口见证差>.46294；'
        +'经典关系标签保留真实量子接触时，未读标签的见证仍>.60647。'
        +f'[结果]({pre}989/shared_record_quantum_action_results.json) · '
        +f'[核验]({pre}989/research_round_989_checks.json)。正式989／累计3774，整体目标未完成。\n\n'
        +f'[机制图v0.5]({pre}989/drafts/mechanism_map_v0_5.md)保留证据政策为内部主体能力，'
        +'不增宇宙收益最大化或所有标签必须相干的原则。停止经典替代与器件分支，'
        +f'接[990]({pre}990/drafts/STATUS.md)回同一量子过程与宏观报告的共同采用。应用目标及957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—988轮共988份','001—989轮共989份')
    if p==paths[4]:s=s.replace('当前正式988／累计3773，988已结项','当前正式989／累计3774，989已结项')
    if p==paths[5]:
        s=s.replace('# 231—988轮阶段成果总览','# 231—989轮阶段成果总览',1)
        s=s.replace('231—988的758份','231—989的759份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至988）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至989）',1)
        rows={
          'C01':'|C01 状态、概率、仪器、复合|1、3、4；001—230、947、955、958—959、989|989复用旧材料全输入界，排除无纠缠及仅经典消息的替代；共享标签仍可选择真实量子作用|四项局部读取须有共同可信字典；未实现全六协议或完整父物理|',
          'C04':'|C04 同一内部动力学与控制|1、6、H2；929、935、947、961—968、987—989|987政策采用为内部有限能力；988区分停活化与撤销，989区分经典标签与量子交互|未增宇宙收益优化或唯一选律要求；控制与物理交互仍需明示|',
          'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—989|989原材料的整个有限窗口见证>.96294，可分替代≤.5；未读经典控制标签仍>.60647|是固定仪器及资源划分下的有效模型检验，非设备无关实验或引力量子化证明|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[989机制图v0.5](../../989/drafts/mechanism_map_v0_5.md)保留共享证据与实际量子作用的分工，不增普遍收益优化或每个关系标签相干要求。停止控制器与经典替代分支，回同一过程与宏观报告的共同采用，复用旧审计而不循环重报。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published 989 mechanism decision to eight live navigation files.')
