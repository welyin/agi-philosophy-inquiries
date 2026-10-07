"""Publish thermal-history adoption; no claim of new baryogenesis or completed goal."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent;HERE=STAGE/'994'
prefix=(STAGE/'993/verify_round993.py').read_text('utf-8').split('    import importlib.util',1)[0]
prefix=prefix.replace('joint candidate 993','charge history 994').replace(
    'research_round_993_checks.json','research_round_994_checks.json').replace('range(776,993)','range(776,994)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core994',HERE/'charge_history.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'charge_history_results.json');core.compare(core.run(),saved)
    assert saved['round']==994 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('actual_cosmological_rates_computed','initial_asymmetry_generated',
                'real_singlet_carries_B_minus_L','thermal_flavor_cases_are_actual_history'):
        assert not saved[key],key
    # Independent reduced scalar equations from Yukawa and charge definitions.
    q=F(7,237);h=-12*q/7;ls=-9*q;b=12*q;l=3*ls-3*h
    assert b-l==1 and b==F(saved['conversion_B_over_B_minus_L'])
    # Check integer/rational nullspaces with an independent numerical rank.
    rr=np.array([[float(F(x)) for x in row] for row in saved['reaction_rows']])
    wy=np.array([float(F(w)*F(y)) for w,y in zip(saved['susceptibility_weights'],saved['hypercharges'])])
    for branch in saved['neutral_equilibrium_branches']:
        rows=[r for r in rr]+[wy]
        for i,j in branch['fast_weinberg_pairs']:
            r=np.zeros(10);r[3+i]+=1;r[3+j]+=1;r[9]+=2;rows.append(r)
        matrix=np.array(rows)
        assert 10-int(np.linalg.matrix_rank(matrix))==branch['neutral_equilibrium_dimension']
        for vec in branch['basis_chemical_potentials']:
            assert np.max(abs(matrix@np.array([float(F(x)) for x in vec])))<1e-12
    # Partial-flavor result without importing the exact solver.
    q=F(7,257);h=-12*q/7;le=-h;lm=(4*q+h-1)/3;lt=(4*q+h)/3
    assert 9*q+le+lm+lt==0
    assert 12*q==F(saved['one_fast_flavor_example']['B']['exact'])
    assert saved['relaxation'][0]['rows'][-1]['B']>3.54e-4
    assert abs(saved['relaxation'][1]['rows'][-1]['B'])<1e-12
    note=STAGE/'research_note_994.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','转换关系','不是宇宙时间','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'boundary_adoption_v1.md',HERE/'charge_history.py',
        HERE/'charge_history_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish994.py',STAGE/'995/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至994）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=994)==list(range(1,995))
    assert '001—994轮共994份' in nav[0].read_text('utf-8-sig')
    assert '231—994的764份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '994：共同热史与净荷保存' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'993/research_round_993_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=994,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=994,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_993=prev['cumulative_numbered_test_groups_from_992']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,neutral_equilibrium_dimensions=[3,2,0],
        adopted_boundary='explicit flavor charges and protected reaction window; generation remains open',
        cosmological_rates_verified=False,initial_asymmetry_generated=False,
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
target=HERE/'verify_round994.py';assert not target.exists();target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 994：共同热史与净荷保存';previous='## 993：共同候选与跨部门匹配'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'+f'[994报告]({pre}research_note_994.md)将成熟电弱转换接到B993：'
        +'声明热窗口内B=28(B−L)/79；不同快速Weinberg反应留下3、2或0个中性平衡方向。'
        +'非平衡不替代偏置的产生与保存，实Z2单态不充当B−L仓库。'
        +f'[结果]({pre}994/charge_history_results.json) · [核验]({pre}994/research_round_994_checks.json)。'
        +'正式994／累计3778，整体未完成。\n\n'
        +f'[共同边界补充]({pre}994/boundary_adoption_v1.md)保留明示初始净荷与保护窗口；'
        +'实际速率、生成源及宇宙丰度未核。停止洗出参数细化，'
        +f'接[995]({pre}995/drafts/STATUS.md)筛选现有内部变量改变中微子耦合的净荷生成机制。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1).replace('001—993轮共993份','001—994轮共994份')
    if p==paths[4]:s=s.replace('当前正式993／累计3777，993已结项','当前正式994／累计3778，994已结项')
    if p==paths[5]:
        s=s.replace('# 231—993轮阶段成果总览','# 231—994轮阶段成果总览',1)
        s=s.replace('231—993的763份','231—994的764份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至993）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至994）',1)
        rows={
          'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—980、994|'
                 '980有限热资源不复用；994将初始味净荷、快反应及保护窗口加入共同边界|'
                 '初始净荷来源未生成；归一化化学松弛不是宇宙实际速率或新的量子动力学|'),
          'C24':('|C24 宇宙初态、Λ与暗部门|1、3、6；629、709—710、964—973、983、991—994|'
                 'B993采用Z2暗部门；994给三代SM热窗口的转换／味选择／洗出条件，区分产生和保存|'
                 '实单态不携带B−L；初始不对称、实际速率、暗丰度及Λ生成解释仍开放|')}
        for cid,row in rows.items():
            s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[994边界补充](../../994/boundary_adoption_v1.md)'
            '分开净荷产生、转换与保存。停止洗出速率及丰度扫描；下一项筛选现有内部变量的'
            '动态Weinberg机制，先验物理CP不变量及同源反作用，不自动加新粒子。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 994 boundary and charge adoption to eight live navigation files.')
