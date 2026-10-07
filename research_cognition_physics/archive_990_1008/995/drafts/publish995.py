"""Publish the internal CP bridge, with explicit limits on generation claims."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent;HERE=STAGE/'995'
prefix=(STAGE/'994/verify_round994.py').read_text('utf-8').split('    import importlib.util',1)[0]
prefix=prefix.replace('charge history 994','internal CP bridge 995').replace(
    'research_round_994_checks.json','research_round_995_checks.json').replace('range(776,994)','range(776,995)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core995',HERE/'internal_cp_bridge.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'internal_cp_bridge_results.json');core.compare(core.run(),saved)
    assert saved['round']==995 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('actual_thermal_kernel_computed','net_lepton_charge_generated',
                'expanding_cosmological_history_verified','full_goal_completed'):
        assert not saved[key],key
    assert saved['added_fundamental_species']==0 and saved['added_effective_operator_dimension']==7
    # Independent exact diagonal contraction, rather than reuse of the trace function.
    invariant=F(1,10)*F(1,5)+F(1,5)*F(1,10)+F(3,10)*F(3,20)
    assert invariant==F(saved['invariant']['exact'])==F(17,200)
    assert abs(saved['CP_conjugate_invariant']+float(invariant))<1e-15
    # Exact squeezed-state moments and quartic correction; no Fock cutoff.
    state=saved['parity_even_initial_state']
    q2=F(state['canonical_Q2']['exact']);p2=F(state['canonical_P2']['exact'])
    assert q2*p2==F(1,4) and F(state['canonical_Q4']['exact'])==3*q2*q2
    assert F(state['energy_including_quartic']['exact'])==F(3400003,6400000)
    assert F(state['f_second_derivative_including_quartic']['exact'])==F(1499997,80000000000)>0
    assert state['mean_sigma']==0
    rows=saved['quadratic_control']['rows']
    assert rows[0]['f']==rows[-1]['f'] and rows[2]['f']>rows[0]['f']
    assert rows[2]['two_time_cp_factor']>0
    # Independent multi-harmonic periodic-shift check of the cycle identity.
    angles=2*np.pi*np.arange(2048)/2048
    def periodic(phi):return np.sin(phi)+.37*np.cos(3*phi)-.21*np.sin(5*phi)
    periodic_shift_residual=max(abs(float(np.mean(periodic(angles-tau)-periodic(angles))))
                                for tau in (.137,1.23,7.11))
    assert periodic_shift_residual<1e-14
    assert abs(saved['causal_periodic_control']['cycle_integral'])<1e-16
    assert saved['causal_periodic_control']['pointwise_max']>1e-6
    note=STAGE/'research_note_995.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==14
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,8)]
    for term in ('整体目标未完成','不是实际Weinberg热碰撞核','不含四次项','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'internal_weinberg_adoption.md',HERE/'internal_cp_bridge.py',
        HERE/'internal_cp_bridge_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish995.py',STAGE/'996/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至995）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [pth for pth in newfiles if pth.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for pth in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',pth.name)
        if m and pth.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=995)==list(range(1,996))
    assert '001—995轮共995份' in nav[0].read_text('utf-8-sig')
    assert '231—995的765份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '995：内部变化与CP偏置' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'994/research_round_994_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=995,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=995,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_994=prev['cumulative_numbered_test_groups_from_993']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,periodic_shift_max_residual=periodic_shift_residual,
        adopted_interface='Z2-even internal Weinberg modulation with explicit dimension-seven input',
        actual_thermal_kernel_computed=False,net_lepton_charge_generated=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(pth.relative_to(ROOT)):sha(pth) for pth in oldfiles},
        new_scientific_and_entry_files={str(pth.relative_to(ROOT)):sha(pth) for pth in newfiles})
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
target=HERE/'verify_round995.py';assert not target.exists();target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 995：内部变化与CP偏置';previous='## 994：共同热史与净荷保存'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'+f'[995报告]({pre}research_note_995.md)以明确新增的Z2偶维七项，'
        +'把已有暗模式与中微子耦合变化相接：零均值的偶态仍可改变二阶矩，'
        +'两份耦合具有不可由常量味换基消去的CP不变量，反作用来自同一作用。'
        +f'[结果]({pre}995/internal_cp_bridge_results.json) · [核验]({pre}995/research_round_995_checks.json)。'
        +'正式995／累计3779，整体未完成。\n\n'
        +f'[W995机制补充]({pre}995/internal_weinberg_adoption.md)只采用内部变化与偏置接口；'
        +'平稳线性记忆近似下，完整周期的平均源抵消。实际热核、净荷与宇宙丰度未核。'
        +f'接[996]({pre}996/drafts/STATUS.md)回共同几何造成的包络及保存窗口，'
        +'先作整体判断，停止相位和器件细化。应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1).replace('001—994轮共994份','001—995轮共995份')
    if p==paths[4]:s=s.replace('当前正式994／累计3778，994已结项','当前正式995／累计3779，995已结项')
    if p==paths[5]:
        s=s.replace('# 231—994轮阶段成果总览','# 231—995轮阶段成果总览',1)
        s=s.replace('231—994的764份','231—995的765份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至994）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至995）',1)
        rows={
          'C18':('|C18 共同作用、来源与反作用|1、4、5、H2；940、952、971、981—986、993、995|'
                 'B993共同源保持；W995新增同一σ²Weinberg作用及σ反作用，质量／干涉／洗出共用耦合|'
                 '新增维七系数为输入；完整相互作用态与热复合源未核，实际材料匹配仍开放|'),
          'C24':('|C24 宇宙初态、Λ与暗部门|1、3、6；629、709—710、964—973、983、991—995|'
                 'W995以已有Z2暗模式连接中微子CP干涉，零均值偶态二阶矩可变，994转换／保存保持|'
                 '平稳线性周期源平均抵消；实际热核、净荷、暗丰度、宇宙边界及Λ解释仍开放|')}
        for cid,row in rows.items():
            s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[W995机制补充](../../995/internal_weinberg_adoption.md)'
            '采用内部变化与物理CP干涉接口，不宣称净荷已生成。停止周期振幅／相位优化；'
            '下一项复用332、964／992核共同几何中的非周期包络及994保存窗口，先判整体价值。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 995 internal CP mechanism and limits to eight live navigation files.')
