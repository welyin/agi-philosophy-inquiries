"""Publish shared vacuum-source and causal-accessibility scope."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent;HERE=STAGE/'997'
old=(STAGE/'996/verify_round996.py').read_text('utf-8')
prefix=old.split('    import importlib.util',1)[0].replace('cosmic CP window 996','vacuum accessibility 997')
prefix=prefix.replace('research_round_996_checks.json','research_round_997_checks.json').replace('range(776,996)','range(776,997)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core997',HERE/'vacuum_accessibility.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'vacuum_accessibility_results.json');core.compare(core.run(),saved)
    assert saved['round']==997 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('physical_vacuum_value_derived','early_boundary_generated',
                'actual_cosmological_fit','full_quantum_geometry_verified','full_goal_completed'):
        assert not saved[key],key
    assert not saved['ideal_causal_tasks']['actual_signal_device_verified']
    assert saved['conditional_future_window']['depends_on_future_continuation']
    # Independent derivative/source calculation from the common cell energy.
    ip=saved['inputs'];M=F(ip['material_mass']['exact']);R=F(ip['radiation_R']['exact'])
    V=F(ip['V0']['exact']);L=F(ip['vacuum_cell_L']['exact']);r=R/L
    assert M==L and F(ip['q']['exact'])==1 and F(ip['r']['exact'])==r
    for row in saved['source_rows']:
        a=F(row['a']);rho=(M/a**3+R/a**4+L)/V;p=(R/(3*a**4)-L)/V
        accel=-(rho+3*p)/(2*L/V)
        assert abs(float(accel)-row['acceleration_over_HLambda_squared'])<1e-15
        assert accel>0 and row['physical_cell_volume_ratio']==int(a**3)
    win=saved['finite_conformal_window'];lo=F(win['lower']['exact']);hi=F(win['upper']['exact'])
    assert F(6577,10000)<lo<hi<F(6582,10000)
    # Dense independent trapezoid is not the certificate; exact monotone sums are.
    x=np.linspace(.25,1,32769)
    estimate=float(np.trapezoid(1/np.sqrt(1+x**3+float(r)*x**4),x))
    assert float(lo)<estimate<float(hi)
    assert 2*F(1,5)<lo and F(11,20)<lo and 2*F(11,20)>1
    fut=saved['conditional_future_window']
    assert F(fut['upper']['exact'])==hi+F(1,4)<1
    assert F(21,10)>2*F(fut['general_upper'])
    # Source matching does not turn a material rest mass into a vacuum density.
    assert F(2)-F(2)*2**3==-14
    for row in saved['availability_same_scale']:assert row['vacuum_shift_residual']<8e-13
'''
tail='    note=STAGE/'+old.split('    note=STAGE/',1)[1]
tail=re.sub(r'99[4-7]',lambda m:str(int(m[0])+1),tail)
for before,after in {
    'mechanism_map_v0_8.md':'cosmological_adoption_v1.md',
    'cosmic_cp_window.py':'vacuum_accessibility.py',
    'cosmic_cp_window_results.json':'vacuum_accessibility_results.json',
    "('整体目标未完成','一维筛选方程','同步取C_ref','386或425')":
        "('整体目标未完成','未来假设','固定几何','386或425')",
    '231—997的766份':'231—997的767份',
    '997：膨胀包络与生成保存约束':'997：有效Λ与实际可达范围',
    'response_verified_with_shared_coupling=True':'finite_causal_window_certified_by_rational_bounds=True',
    "adopted_interface='conditional nonperiodic envelope with shared source-washout parameters and H/T suppression'":
        "adopted_interface='effective positive vacuum source with explicitly limited causal access and future assumptions'",
    'thermal_source_coefficient_evaluated=False,finite_physical_net_charge_certified=False':
        'physical_vacuum_value_derived=False,actual_cosmological_future_predicted=False'
}.items():
    assert before in tail,before
    tail=tail.replace(before,after)
target=HERE/'verify_round997.py';assert not target.exists();target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 997：有效Λ与实际可达范围';previous='## 996：膨胀包络与生成保存约束'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'+f'[997报告]({pre}research_note_997.md)将正有效Λ作为明示的加速机制，'
        +'区分参数未定与机制缺失。同一材料／辐射／真空来源给有限类光窗口；'
        +'体积增长不自动增加可完成的回读权限。'
        +f'[结果]({pre}997/vacuum_accessibility_results.json) · [核验]({pre}997/research_round_997_checks.json)。'
        +'正式997／累计3781，整体未完成。\n\n'
        +f'[B997宇宙补充]({pre}997/cosmological_adoption_v1.md)保共同真空账、有限任务与未来条件；'
        +'未预测Λ值、实际宇宙未来或全量子几何。停止真空／视界细化，'
        +f'接[998]({pre}998/drafts/STATUS.md)先对齐同一宇宙状态、约束和记录资源的边界合同。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1).replace('001—996轮共996份','001—997轮共997份')
    if p==paths[4]:s=s.replace('当前正式996／累计3780，996已结项','当前正式997／累计3781，997已结项')
    if p==paths[5]:
        s=s.replace('# 231—996轮阶段成果总览','# 231—997轮阶段成果总览',1)
        s=s.replace('231—996的766份','231—997的767份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至996）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至997）',1)
        rows={
          'C05':('|C05 因果序、局域性、传播界|1、4；465—466、933、957、970、997|'
                 '旧传播界保持；997在给定FLRW中给有限回读域及带未来假设的可达限制|'
                 '不是从认知生成光锥；理想因果许可不等于实际仪器，局部加速不指定全部未来|'),
          'C24':('|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、332、553、601、613、964—973、983、991—997|'
                 '997采用正有效Λ的加速机制及共同真空账；D991与W995—996按原条件保留|'
                 'Λ取值／自然性、早期边界、暗丰度及实际净荷未推得；体积增长不等于可访问容量|')}
        for cid,row in rows.items():
            s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[B997宇宙补充](../../997/cosmological_adoption_v1.md)'
            '采用明示正Λ的有效机制，区分参数、状态边界与实际可达量词。停止真空／视界／参数优化；'
            '下一项复用既有约束态和准备工具，对齐同一宇宙边界及记录资源，不要求制造全部过去。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 997 common vacuum and accessibility scope to eight navigation files.')
