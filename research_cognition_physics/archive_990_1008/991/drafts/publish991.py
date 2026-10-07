"""Publish bounded dark-sector mechanism screen; preserve prior frozen work."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'991'
prefix=(STAGE/'990/verify_round990.py').read_text('utf-8').split('    inventory=read(',1)[0]
prefix=prefix.replace('mechanism synthesis 990','mechanism screen 991')
prefix=prefix.replace('research_round_990_checks.json','research_round_991_checks.json')
prefix=prefix.replace('range(776,990)','range(776,991)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core991',HERE/'dark_mode_screen.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'dark_mode_screen_results.json');core.compare(core.run(),saved)
    assert saved['round']==991 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert saved['scalar_is_new_explicit_physical_input']
    assert not saved['scalar_model_is_original']
    assert not saved['current_experimental_viability_claimed']
    assert not saved['relic_density_computed']
    assert F(saved['masses_squared']['dark']['exact'])==F(1,4)
    assert F(saved['analytic_source_energy_ratio_upper']['exact'])==F(1,500)
    # Independent larger representation verifies that moments do not use the edge.
    a=np.diag(np.sqrt(np.arange(1,12,dtype=float)),1)
    q=(a+a.T)/np.sqrt(1000.)
    residuals=[]
    for row in saved['rows']:
        n=row['n']
        residuals.append(abs((q@q)[n,n]-row['s2']['value']))
        residuals.append(abs(np.linalg.matrix_power(q,4)[n,n]-row['s4']['value']))
    assert max(residuals)<1e-13
    note=STAGE/'research_note_991.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,7)]
    for term in ('整体目标未完成','不是两种测量概率之比','没有识别现实暗物质',
                 '不是用8维截断认证长期量子场演化','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'dark_mode_screen.py',HERE/'dark_mode_screen_results.json',
        Path(__file__),HERE/'drafts/selection.md',HERE/'drafts/mechanism_map_v0_6.md',
        HERE/'drafts/publish991.py',STAGE/'992/drafts/STATUS.md']
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至991）',1)[1].split('### 当前取舍',1)[0]
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
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=991)==list(range(1,992))
    assert '001—991轮共991份' in nav[0].read_text('utf-8-sig')
    assert '231—991的761份' in nav[5].read_text('utf-8-sig')
    for p in nav:assert '991：共同几何与弱可见模式' in p.read_text('utf-8-sig')
    oldfiles=[STAGE/'990/research_round_990_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=991,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=991,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_990=prev['cumulative_numbered_test_groups_from_989']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_fock_moment_max_residual=max(residuals),
        mechanism_decision='retain weak visible access with a shared full geometric source as optional D991',
        original_dark_matter_model_claimed=False,physical_dark_matter_identified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
target=HERE/'verify_round991.py';assert not target.exists()
target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 991：共同几何与弱可见模式';previous='## 990：整体假说与联合解释优先'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[991报告]({pre}research_note_991.md)复用旧共同依赖，明确普适几何反馈不要求各内部通道同样可见。'
        +'以成熟Z2实单态门户给可选暗部门D991一个稳定、弱直接接触、完整应力的机制；'
        +'瞬时矩与同源交换完成有界复算。'
        +f'[结果]({pre}991/dark_mode_screen_results.json) · [核验]({pre}991/research_round_991_checks.json)。'
        +'正式991／累计3775，整体目标未完成。\n\n'
        +f'[机制图v0.6]({pre}991/drafts/mechanism_map_v0_6.md)保留新增字段／对称／初态输入；'
        +'未认定现实暗物质，未验丰度或实验允许范围。停止门户扩建，'
        +f'接[992]({pre}992/drafts/STATUS.md)补宇宙状态、记录方向与几何演化的整体机制。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—990轮共990份','001—991轮共991份')
    if p==paths[4]:
        s=s.replace('当前正式990／累计3774，990机制整合已交付','当前正式991／累计3775，991已结项')
    if p==paths[5]:
        s=s.replace('# 231—990轮阶段成果总览','# 231—991轮阶段成果总览',1)
        s=s.replace('231—990的760份','231—991的761份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至990）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至991）',1)
        row=('|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、332—333、553、601、964—973、983、991|'
             '991用明示Z2实单态扩展连接稳定、弱普通门户与共同几何能源；保留为暗部门机制候选D991|'
             '未认定现实暗物质；初态生成、丰度、聚集、现实实验、暗能量与Λ解释仍开放，停止门户细化|')
        s,count=re.subn(r'^\|C24 [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[991机制图v0.6](../../991/drafts/mechanism_map_v0_6.md)'
            '在990共同比较机制上加入可选D991：普通通道弱可见与共同几何能源并存。'
            '不是暗物质认定；停止门户、寿命及丰度扩建，回宇宙边界、热箭头和几何演化的整体解释。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 991 mechanism decision to eight live navigation files.')
