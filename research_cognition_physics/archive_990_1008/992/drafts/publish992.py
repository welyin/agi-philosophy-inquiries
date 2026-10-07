"""Publish differential-cooling mechanism, not a proof of the cosmological arrow."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'992'
prefix=(STAGE/'991/verify_round991.py').read_text('utf-8').split('    import importlib.util',1)[0]
prefix=prefix.replace('mechanism screen 991','cooling mechanism 992')
prefix=prefix.replace('research_round_991_checks.json','research_round_992_checks.json')
prefix=prefix.replace('range(776,991)','range(776,992)')
checks=r'''
    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location('core992',HERE/'cooling_availability.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'cooling_availability_results.json');core.compare(core.run(),saved)
    assert saved['round']==992 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('physical_reset_implemented','work_extraction_implemented',
                'initial_cosmological_boundary_explained','full_quantum_geometry_verified'):
        assert not saved[key],key
    assert saved['additional_blank_registers_created']==0
    # Original rational/radical four-level matrix, independent of the material importer.
    h=np.array([[0.,-.2,0.,-math.sqrt(3)/80],[-.2,1.,0.,0.],
                [0.,0.,1.,0.],[-math.sqrt(3)/80,0.,0.,-.025]])
    ev,u=np.linalg.eigh(h);beta0=2*math.log(2)
    p=np.exp(-beta0*ev);p/=sum(p)
    actual_entropy=float(-p@np.log(p));e0=float(p@ev)
    errors=[]
    for row in saved['rows']:
        beta=beta0*row['a'];q=np.exp(-beta*ev);z=sum(q);q/=z
        # F(tau) = -T log Z, with unshifted physical h.
        A=e0-actual_entropy/beta+math.log(z)/beta
        errors.append(abs(A-row['availability']))
        assert abs(float(np.sum(abs(p-q))/2)-row['equilibrium_state_trace_distance'])<1e-12
    assert max(errors)<1e-12
    assert abs(saved['expansion_availability_limit']-(e0-min(ev)))<1e-12
    note=STAGE/'research_note_992.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','不是材料已经达到的状态','不自动提取工作',
                 '两者还不是已经共同实现的一次宇宙记录过程','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'cooling_availability.py',HERE/'cooling_availability_results.json',
        Path(__file__),HERE/'drafts/selection.md',HERE/'drafts/mechanism_map_v0_7.md',
        HERE/'drafts/publish992.py',STAGE/'993/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至992）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=992)==list(range(1,993))
    assert '001—992轮共992份' in nav[0].read_text('utf-8-sig')
    assert '231—992的762份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '992：差异冷却与非平衡资源' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'991/research_round_991_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=992,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=992,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_991=prev['cumulative_numbered_test_groups_from_990']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_free_energy_max_residual=max(errors),
        mechanism_decision='retain differential scaling as athermality formation; keep actual use and boundary selection separate',
        cosmological_arrow_derived=False,full_quantum_geometry_verified=False,
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
target=HERE/'verify_round992.py';assert not target.exists()
target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 992：差异冷却与非平衡资源';previous='## 991：共同几何与弱可见模式'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[992报告]({pre}research_note_992.md)将旧固定材料谱与自由光红移接到305相对熵身份：'
        +'初始同温后，差异尺度响应形成有界非平衡自由能；原980材料在a=2给A≈.04299333，'
        +'熵仍保持，膨胀侧上限≈.22019126。'
        +f'[结果]({pre}992/cooling_availability_results.json) · [核验]({pre}992/research_round_992_checks.json)。'
        +'正式992／累计3776，整体目标未完成。\n\n'
        +f'[机制图v0.7]({pre}992/drafts/mechanism_map_v0_7.md)区分资源形成、实际使用与初始边界；'
        +'未完成工作提取、记忆复位或宇宙箭头。停止热装置／温度扩建，'
        +f'接[993]({pre}993/drafts/STATUS.md)把近期机制并入一份候选，先核共同身份和解释范围。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—991轮共991份','001—992轮共992份')
    if p==paths[4]:s=s.replace('当前正式991／累计3775，991已结项','当前正式992／累计3776，992已结项')
    if p==paths[5]:
        s=s.replace('# 231—991轮阶段成果总览','# 231—992轮阶段成果总览',1)
        s=s.replace('231—991的761份','231—992的762份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至991）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至992）',1)
        row=('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—980、992|'
             '992原材料的差异冷却可形成正且有界的非平衡自由能，实际态熵保持；资源形成与980的使用合同分开|'
             '隔离、稳定谱、初始同温与几何分支输入；未共同实现宇宙记录、实际复位或初始箭头来源|')
        s,count=re.subn(r'^\|C23 [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[992机制图v0.7](../../992/drafts/mechanism_map_v0_7.md)'
            '加入差异尺度响应形成非平衡资源，区分资源使用与初始边界。'
            '停止局部热设备与门户扩建；把990—992与P981并入同一明确候选，先核重复计数和共同采用。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 992 mechanism decision to eight live navigation files.')
