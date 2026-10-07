"""Publish the actual native dipole-metric bridge without changing frozen files."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'985/verify_round985.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"two_active_sources_results.json")',1)[0]
for a,b in [('range(776,985)','range(776,986)'),('Delivery checks for 985','Delivery checks for 986'),
            ('research_round_985_checks.json','research_round_986_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"metric_dipole_bridge_results.json")
    assert result["round"]==986 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    import numpy as np
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location("core986",HERE/"metric_dipole_bridge.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    bound=result["infinite_domain_bound"]
    assert Fraction(bound["rational_Duhamel_upper"]["exact"])<Fraction(134,10**9)
    assert Fraction(bound["record_contrast_lower"]["exact"])>Fraction(4859478,10**7)
    assert bound["all_times_in_original_interval"] and bound["unknown_input_and_passive_reference"]
    # Independent second source moment check: retain cross terms explicitly.
    material=mod.load("native_check986",STAGE/"968/internal_relay.py")
    m,h,q,_,w,_=material.material();n=5
    a=np.diag(np.sqrt(np.arange(1,n)),1);I=np.eye(4)
    S=mod.kron(q,I)+mod.kron(I,q)
    Cbare=.5*mod.kron(np.eye(16),a@a+a.T@a.T)
    mixed=.003*mod.kron(S,a+a.T);C=Cbare+mixed
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    plus=w@np.ones(2)/math.sqrt(2)
    W=np.column_stack([np.kron(np.kron(w[:,i],plus),photon) for i in (0,1)])
    difference=(C@W).T@(C@W)-(Cbare@W).T@(Cbare@W)
    d=(1-1/math.sqrt(1.16))/2
    expected=.003**2*d*np.diag([3.,1.])
    assert np.linalg.norm(difference-expected)<1e-14
    assert np.linalg.norm(W.T@(Cbare@mixed+mixed@Cbare)@W)<1e-14
    for key in ("independent_quadrupole_matching_derived","full_U1_completed",
        "spatial_mode_matching_error_certified","full_GR_error_certified",
        "Newton_motion_985_and_TT_glued","all_SM_and_cosmology_completed",
        "physical_minimum_scale_assumed","full_goal_completed"):
        assert result["scope"][key] is False
    note=STAGE/"research_note_986.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","停止本模式、耦合和精度优化",
        "不声称该延拓恢复强场GR","独立固有四极在所选TT偏振上的投影为零",
        "不等于所有无界高阶物理余项都有同样误差界"):
        assert term in prose,term
    newfiles=[note,HERE/"metric_dipole_bridge.py",HERE/"metric_dipole_bridge_results.json",
        Path(__file__),HERE/"drafts/metric_dipole_decision.md",HERE/"drafts/metric_dipole_adoption.md",
        HERE/"drafts/publish986.py",STAGE/"987/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [('截至985','截至986'),('n<=985','n<=986'),('list(range(1,986))','list(range(1,987))'),
            ('001—985轮共985份','001—986轮共986份'),('231—985的755份','231—986的756份'),
            ('985：两个活动来源与有限域共同过程','986：原电磁交互与传播几何的共同生成元'),
            ('STAGE/"984/research_round_984_checks.json"','STAGE/"985/research_round_985_checks.json"')]:tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=986,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=986,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_985=prev["cumulative_numbered_test_groups_from_984"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,native_EM_metric_source_joined=True,
        independent_mixed_source_second_moment_verified=True,
        full_Fock_finite_time_record_bound_verified=True,
        entire_U1_completed=False,full_physical_error_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'986/verify_round986.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 986：原电磁交互与传播几何的共同生成元'
    previous='## 985：两个活动来源与有限域共同过程'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[986报告]({pre}research_note_986.md)从同一模式作用将976原电磁交互接入动态TT过程；'
        +'几何来源含原偶极固定的混合项。原记录全区间的联合态误差<1.34×10⁻⁷，'
        +'记录对比>.4859478。'
        +f'[结果]({pre}986/metric_dipole_bridge_results.json) · '
        +f'[核验]({pre}986/research_round_986_checks.json)。正式986／累计3771，整体目标未完成。\n\n'
        +f'[采用范围]({pre}986/drafts/metric_dipole_adoption.md)保模式／偏振／支撑输入；'
        +'未将985运动、982非零固有四极及完整父理论自动合并。停止本模式和精度优化，'
        +f'接[987]({pre}987/drafts/STATUS.md)回共同恢复对象及采用范围。应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—985轮共985份','001—986轮共986份')
    if p==paths[4]:s=s.replace('当前正式985／累计3770，985已结项','当前正式986／累计3771，986已结项')
    if p==paths[5]:
        s=s.replace('# 231—985轮阶段成果总览','# 231—986轮阶段成果总览',1)
        s=s.replace('231—985的755份','231—986的756份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至985）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至986）',1)
        rows={
        'C18':'|C18 共同作用、来源与反作用|1、4、5、H2；940、971、981—986|986从实际偶极—模式作用给同一生成元、lapse及TT混合来源；系数由旧g固定|只核声明模式及零固有四极投影部门，未将985运动或完整约束自动合入|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—986|986来源二阶矩有实际漏项对照；原976记录在动态TT过程有无限占据全时间界|误差只属所列模式，未认证父物理空间／支撑匹配，非实际引力探测证书|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[986采用范围](../../986/drafts/metric_dipole_adoption.md)加入原电磁交互的几何来源与共同正过程，整体U1—U3仍未结项。停止模式／耦合／精度优化，接[987](../../987/drafts/STATUS.md)回共同恢复对象、参数与有效域，不把不同合同的数值拼成总证书。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 986 verifier.')
