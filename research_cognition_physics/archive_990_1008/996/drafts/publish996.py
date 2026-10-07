"""Publish the finite cosmological envelope screen, not a physical abundance claim."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent;HERE=STAGE/'996'
old=(STAGE/'995/verify_round995.py').read_text('utf-8')
prefix=old.split('    import importlib.util',1)[0].replace('internal CP bridge 995','cosmic CP window 996')
prefix=prefix.replace('research_round_995_checks.json','research_round_996_checks.json').replace('range(776,995)','range(776,996)')
checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core996',HERE/'cosmic_cp_window.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'cosmic_cp_window_results.json');core.compare(core.run(),saved)
    assert saved['round']==996 and saved['all_scientific_checks_passed']
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('thermal_source_coefficient_evaluated','finite_physical_net_charge_certified',
                'actual_thermal_hierarchy_verified','cosmological_abundance_fitted','full_goal_completed'):
        assert not saved[key],key
    state=saved['phase_averaged_quadratic_preparation']
    assert state['constant_coefficient_shifted_together']
    assert not state['quartic_and_thermal_corrections_certified']
    assert F(state['occupation']['exact'])==F(9,16)
    assert F(state['f_excitation_initial']['exact'])==F(9,800000)
    assert F(state['initial_vacuum_reference_shift']['exact'])==F(1,100000)
    # Exact initial full-reference identity and trace polynomial.
    fi=F(9,800000);fv=F(1,100000);r2=F(233,2500)
    assert fi+fv==F(17,800000)
    shared=saved['shared_coupling'];A=F(shared['norm_constant']['exact'])
    B=F(shared['norm_a_minus3']['exact']);D=F(shared['norm_a_minus6']['exact'])
    gi=F(shared['norm_initial']['exact'])
    assert A==F(7,50)+fv/F(100)+fv*fv*r2
    assert B==fi/F(100)+2*fi*fv*r2 and D==fi*fi*r2 and gi==A+B+D
    window=saved['finite_window'];W=F(window['optical_depth']['exact'])
    # Independent coefficient integration on [1,2], exact arithmetic.
    assert W==(A/2+B*F(15,64)+D*F(127,896))/(10*gi)
    assert F(window['unwashed_shape']['exact'])==F(31,160)
    lower=math.exp(-float(W))*float(F(31,160));y=window['normalized_response']
    assert 0<lower<=y<=float(F(31,160))
    expected=3*float(F(17,200)*fi/gi/F(10))*y
    assert abs(window['yield_per_abs_alpha_over_beta_times_H_over_T']-expected)<1e-18
    assert window['independent_ODE_residual']<1e-12
    assert window['quadrature_comparison']<1e-13
    # Derivative of cumulative washout equals minus the rate/H/a, independently.
    for a in (1.1,1.4,1.8):
        h=1e-5
        def integral(x):return float((A*(F(1)/F(str(x))-F(1,2))
            +B*(F(str(x))**-4-F(1,16))/4+D*(F(str(x))**-7-F(1,128))/7)/(10*gi))
        deriv=(integral(a+h)-integral(a-h))/(2*h)
        rhs=-.1*(float(A)+float(B)*a**-3+float(D)*a**-6)/(float(gi)*a*a)
        assert abs(deriv-rhs)<1e-10
    assert F(saved['counter_kernel_first_moment']['exact'])==0
'''
tail='    note=STAGE/'+old.split('    note=STAGE/',1)[1]
tail=re.sub(r'99[3-6]',lambda m:str(int(m[0])+1),tail)
for before,after in {
    'internal_weinberg_adoption.md':'mechanism_map_v0_8.md',
    'internal_cp_bridge.py':'cosmic_cp_window.py',
    'internal_cp_bridge_results.json':'cosmic_cp_window_results.json',
    "prose.count('$$')==14":"prose.count('$$')==16",
    'range(1,8)':'range(1,9)',
    "('整体目标未完成','不是实际Weinberg热碰撞核','不含四次项','386或425')":
        "('整体目标未完成','一维筛选方程','同步取C_ref','386或425')",
    '231—996的765份':'231—996的766份',
    '996：内部变化与CP偏置':'996：膨胀包络与生成保存约束',
    'periodic_shift_max_residual=periodic_shift_residual':'response_verified_with_shared_coupling=True',
    "adopted_interface='Z2-even internal Weinberg modulation with explicit dimension-seven input'":
        "adopted_interface='conditional nonperiodic envelope with shared source-washout parameters and H/T suppression'",
    'actual_thermal_kernel_computed=False,net_lepton_charge_generated=False':
        'thermal_source_coefficient_evaluated=False,finite_physical_net_charge_certified=False'
}.items():
    assert before in tail,before
    tail=tail.replace(before,after)
target=HERE/'verify_round996.py';assert not target.exists();target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
    STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                     '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 996：膨胀包络与生成保存约束';previous='## 995：内部变化与CP偏置'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'+f'[996报告]({pre}research_note_996.md)把内部占据稀释接到W995：'
        +'非周期包络与有限保存可在声明归约内相容；同一耦合的生成／洗出关系带H/T慢变限制。'
        +'归一化响应非零不等于现实净荷或足够产额。'
        +f'[结果]({pre}996/cosmic_cp_window_results.json) · [核验]({pre}996/research_round_996_checks.json)。'
        +'正式996／累计3780，整体未完成。\n\n'
        +f'[整体机制图v0.8]({pre}996/mechanism_map_v0_8.md)保留条件性宇宙链及新增输入；'
        +'实际热核、味输运、共同热窗口与丰度未核。停止该候选的产额细化，'
        +f'接[997]({pre}997/drafts/STATUS.md)先回查宇宙边界、Λ与共同真空来源的整体缺口。'
        +'应用目标及957保持。\n\n')
    s=s.replace(previous,block+previous,1).replace('001—995轮共995份','001—996轮共996份')
    if p==paths[4]:s=s.replace('当前正式995／累计3779，995已结项','当前正式996／累计3780，996已结项')
    if p==paths[5]:
        s=s.replace('# 231—995轮阶段成果总览','# 231—996轮阶段成果总览',1)
        s=s.replace('231—995的765份','231—996的766份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至995）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至996）',1)
        rows={
          'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—980、994—996|'
                 '996把占据差与参考系数共同匹配；初始因子化的领先C₀C₂干涉只需一次σ²均值|'
                 '冷占据、相位混合及热层级为输入；高阶关联、真实复合基线与共同准备未核|'),
          'C24':('|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、332、553、601、964—973、983、991—996|'
                 '996条件性非周期包络与有限保存相容；生成／洗出共用参数并显式受H/T限制|'
                 '真实净荷、热核、味窗口及暗丰度未核；停止候选细化，回早期边界与Λ整体机制|')}
        for cid,row in rows.items():
            s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement=('4. 957数学假说v0.2保持；[996整体机制图](../../996/mechanism_map_v0_8.md)'
            '保留非周期CP包络的条件桥及慢变幅度限制，不认定现实净荷生成。停止热核／丰度／相位细化；'
            '下一项回查宇宙边界、Λ与共同真空来源，先分辨数值未定、机制缺失和域外延伸。')
        s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items()),'Concurrent navigation edit'
for p,data in updates.items():p.write_bytes(data)
print('Published 996 finite envelope screen and whole-model adoption to eight navigation files.')
