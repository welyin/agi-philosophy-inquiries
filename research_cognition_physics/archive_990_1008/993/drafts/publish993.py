"""Publish the shared candidate decision; no claim of completed unification."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'993'
prefix=(STAGE/'992/verify_round992.py').read_text('utf-8').split('    import importlib.util',1)[0]
prefix=prefix.replace('cooling mechanism 992','joint candidate 993')
prefix=prefix.replace('research_round_992_checks.json','research_round_993_checks.json')
prefix=prefix.replace('range(776,992)','range(776,993)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    spec=importlib.util.spec_from_file_location('core993',HERE/'joint_matching.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'joint_matching_results.json');core.compare(core.run(),saved)
    assert saved['round']==993 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('native_higgs_matching_computed','prospective_budgets_physically_certified',
                'interacting_resonance_claimed','full_common_process_verified'):
        assert not saved[key],key
    old=read(STAGE/'991/dark_mode_screen_results.json')
    mh2=F(old['masses_squared']['higgs']['exact'])
    md2=F(old['masses_squared']['dark']['exact'])
    assert mh2==4*md2 and saved['free_pair_frequency_equals_higgs_mass']
    # Independent finite-difference mixed derivative of the explicit potential.
    mixed=[]
    for row in saved['static_source_rows']:
        j=F(row['ordinary_calibration_source']['exact'])
        d=F(row['dark_source_contrast']['exact'])
        z1=-j/mh2;z2=-(j+d)/mh2;z3=-d/mh2
        E=lambda z,source: mh2*z*z/2+source*z
        delta=E(z2,j+d)-E(z1,j)-E(z3,d)
        assert delta==F(row['cross_energy_density']['exact'])
        mixed.append(str(delta))
    oldthermal=read(STAGE/'992/cooling_availability_results.json')
    a2=next(x for x in oldthermal['rows'] if x['a']==2)
    th=saved['thermal_interface'];delta=th['prospective_state_trace_error']
    ee=read(STAGE/'980/finite_thermal_records_results.json')['parameters']['energies']
    fdelta=-delta*math.log(delta)-(1-delta)*math.log(1-delta)+delta*math.log(3)
    expected=(2*th['prospective_operator_error']+(max(ee)-min(ee))*delta
        +(a2['temperature']+th['prospective_temperature_error'])*fdelta
        +th['prospective_temperature_error']*math.log(4))
    assert abs(expected-th['error_bound'])<1e-12
    assert a2['availability']-expected>.03525724
    note=STAGE/'research_note_993.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==8
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,5)]
    for term in ('整体目标未完成','选择预算不是证明达标','复用','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'common_candidate_v1.md',HERE/'joint_matching.py',
        HERE/'joint_matching_results.json',Path(__file__),HERE/'drafts/selection.md',
        HERE/'drafts/publish993.py',STAGE/'994/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至993）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=993)==list(range(1,994))
    assert '001—993轮共993份' in nav[0].read_text('utf-8-sig')
    assert '231—993的763份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '993：共同候选与跨部门匹配' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'992/research_round_992_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=993,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=993,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_992=prev['cumulative_numbered_test_groups_from_991']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        conditional_availability_lower_bound=a2['availability']-expected,
        candidate_version='B993 v1: P981 + D991; retain common Higgs dynamics',
        matching_budgets_certified=False,full_common_process_verified=False,
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
target=HERE/'verify_round993.py';assert not target.exists()
target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 993：共同候选与跨部门匹配';previous='## 992：差异冷却与非平衡资源'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[993报告]({pre}research_note_993.md)将P981＋D991及热资源机制写成'
        +f'[共同候选B993]({pre}993/common_candidate_v1.md)：普通与暗源必须共同匹配，'
        +'保交叉作用及来源；991自由频率在成对阈值，不能自动消去Higgs。'
        +f'[结果]({pre}993/joint_matching_results.json) · [核验]({pre}993/research_round_993_checks.json)。'
        +'正式993／累计3777，整体目标未完成。\n\n'
        +'992资源机制有非零的条件性误差窗口，但所选谱／态／温度预算尚未物理认证。'
        +'停止门户和热设备细化，'
        +f'接[994]({pre}994/drafts/STATUS.md)回共同宇宙边界、非平衡与物质不对称的机制。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—992轮共992份','001—993轮共993份')
    if p==paths[4]:s=s.replace('当前正式992／累计3776，992已结项','当前正式993／累计3777，993已结项')
    if p==paths[5]:
        s=s.replace('# 231—992轮阶段成果总览','# 231—993轮阶段成果总览',1)
        s=s.replace('231—992的762份','231—993的763份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至992）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至993）',1)
        rows={
          'C18':('|C18 共同作用、来源与反作用|1、4、5、H2；940、952、971、981—986、993|'
                 'B993采用P981＋D991，共同Higgs源先合并再匹配，普通／暗交叉作用与来源共同保留|'
                 '源校准不是原生材料Higgs匹配；显式传播与已消去接触不重计，完整共同过程未核|'),
          'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—980、992—993|'
                 '992差异冷却形成有界非平衡资源；993给谱／态／温度误差下的条件性稳定窗口|'
                 '选定误差预算尚未物理认证；资源使用与宇宙边界未共同实现，不自动产生箭头|'),
          'C24':('|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、332—333、553、601、964—973、983、991—993|'
                 'B993明确采用Z2暗部门，保普通材料反馈和共同Higgs动力学；991自由阈值限制重场消去|'
                 '并非相互作用共振或现实暗物质认证；初态、物质不对称、丰度、Λ生成解释继续开放|')}
        for cid,row in rows.items():
            s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[B993共同候选](../../993/common_candidate_v1.md)'
            '明确采用P981＋D991，谱、源及热资源共同匹配。'
            '停止门户／共振／热设备细化；下一项回共同边界、非平衡和物质不对称，先查机制与独立输入。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published B993 and its adoption scope to eight live navigation files.')
