"""Prepare672 without modifying frozen evidence or the active goal."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_671.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：671静态规范曲率与原联合反射泛函',
    '# 联合条件总账：672独立规范历史与共同电动能拼接')
ledger=ledger.replace('接[670全账](unified_physics_condition_ledger_670.md)，回填[671报告](research_note_671.md)。[结果](joint_static_gauge_reflection_results.json)、[核验](research_round_671_checks.json)。',
    '接[671全账](unified_physics_condition_ledger_671.md)，回填[672报告](research_note_672.md)。[结果](joint_gauge_history_kernel_results.json)、[核验](research_round_672_checks.json)。')
ledger=ledger.replace('671静态空间曲率下的原联合未归一反射泛函|',
    '671静态联合反射；672独立历史与电热核的共同必要限制|')
ledger=ledger.replace('群选择与动态规范过程仍缺|',
    '672原商群电热核与手征核须联合验收；群选择、共享边界与动态Gauss过程仍缺|')
ledger=ledger.replace('动态规范正性、指定归一与全局测度仍缺|',
    '672固定边界的联合裸核出现数值负方向，正电热核也不自动修复；动态Gauss正性、指定归一与全局测度仍缺|')
ledger=ledger.replace('手征辅助过程与原全图物理时间的识别；实际无限动力学存在|',
    '672固定边界二步候选与643原Gauss有序过程仍须匹配；实际无限动力学存在|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 672不改变671静态定理；新负方向位于不同规范半区历史之间。固定共享时间链路的裸核及二步电热核候选不等于完整Gauss核，数值候选失败不能扩大为原H_F失败。

## 本轮合并与下一项

672连接C04、C14、C16、C17、C19、C22：保留原完整相位、16通道及全部质量，检查不同半区规范历史之间的联合核。其固定背景物理Gram均正，裸联合标量核却有最小特征值约−.43191；帧重相位和原大矩阵Pfaffian身份已核。数值见证不称区间认证证明。

共同K=B∘F的解析二阶主子式条件要求电因子相关强度不超过原手征核所允许的值。使用原Z6商群完整表示热核，指定对角电度量及二步分割下，τ=.04的有限样本通过、τ=.25失败。这是共同条件限制，不引入自由U1替代物，也不优化τ作为现实参数。

独立正玻色核和静态物理正性不足以签收共同动态理论。原共享边界、Gauss端点／费米闭合、完整有序影响与来源仍须同步；643原过程已存在，不重复重建一般路径积分。原非对角电度量、完整势及极限不由本项分支覆盖。

接[673](round673_drafts/STATUS.md)：核固定边界见证与规范不变观测的实际关系，恢复共同Gauss闭合，再比较原有序影响。不要继续扫描热时或增加对角样本。空间旧定理及已消去条件完整继承；四分支未统一，目标不改。
'''
write('unified_physics_condition_ledger_672.md',ledger)
mapping={'joint_static_gauge_reflection':'joint_gauge_history_kernel',
         '670':'671','671':'672','672':'673','3247':'3249','3249':'3251',
         '1322':'1325','1325':'1328','2376':'2389','2389':'2403'}
verify=replace((HERE/'verify_round671.py').read_text('utf8'),mapping)
verify=verify.replace('Verify static noncommuting spectral density and original joint reflection.',
    'Verify independent signed half histories and original-group electric sewing.')
verify=verify.replace('import joint_nonflat_mass_measure as prior_model','import joint_static_gauge_reflection as prior_model')
verify=verify.replace('magnetic_reflection_probe','dynamic_commutator_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_672.md','round672_drafts/research_note_672_draft.md',
           'round672_drafts/final_review.txt','round673_drafts/STATUS.md',
           'round672_drafts/dynamic_commutator_entry.md','round672_drafts/dynamic_commutator_probe.py',
           'round672_drafts/dynamic_commutator_probe_results.json','round672_drafts/half_history_probe.py',
           'round672_drafts/half_history_probe_results.json','round672_drafts/first_heat_attempt.py',
           'round672_drafts/first_heat_results.json')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==18","assert text['display_formulas']==16")
write('verify_round672.py',verify)
publish=replace((HERE/'publish_round671.py').read_text('utf8'),mapping|{'## 317.':'## 318.','## 222.':'## 223.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第672轮完成：** [独立规范历史、原手征权重与电动能的共同拼接条件]({p}research_note_672.md)'
         '原带相位权重在独立半区历史间出现数值负方向；解析主子式条件限制共同电热核，原商群实际热核复现通过与失败分支。'
         '固定边界候选不等于完整Gauss过程。'
         '两组、十六式通过，最新672／3251，1328份编号科学文件、2403份保护证据。'
         '[核验]({p}research_round_672_checks.json)、[全条件账]({p}unified_physics_condition_ledger_672.md)。'
         '共享边界、原完整过程、连续及量子GR仍开放。')
order=('**当前执行顺序（672后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[673共享边界与原Gauss闭合]({p}round673_drafts/STATUS.md)，'
       '核负方向的实际物理代数与完整有序过程，停止热时扫描；目标不改。')
'''+publish[b:]
publish=publish.replace('静态规范曲率、原手征物质与共同反射泛函','独立规范历史、原手征权重与电动能的共同拼接条件')
publish=publish.replace('静态曲率与原物理辅助质量的共同反射条件','独立规范历史与原电动能的共同限制')
publish=publish.replace('静态反射泛函与真实动态规范过程的边界','固定边界核与完整Gauss过程的边界')
write('publish_round672.py',publish)
write('postcheck_round672.py',replace((HERE/'postcheck_round671.py').read_text('utf8'),mapping))
write('round672_drafts/research_note_672_draft.md',(HERE/'research_note_672.md').read_text('utf8'))
review='''672 primary-agent mathematical/code/scope review; no independent agent.
Goal unchanged; no tasks, automations, image checks, installs, or old evidence edits.
Reused382-386/425/522-523 via667 full mapping. No restored Lipschitz, standard
antipode, half-homomorphism, global contraction or repeated joint-existence gap.
671 static RP not invalidated:672 compares independent gauge halves, not merely
another reflection-compatible diagonal. Actual sixteen channels, full original
mass, auxiliary E, frame determinant and Pfaffian phase are retained.
Original670 projector dictionary used; direct large Pfaffian, source covariance,
triangular identity and frame rephasing verified. No entrywise absolute weights.
Diagonal congruence is inertia preserving. Bare scalar kernel negative direction
is a numerical finite witness with stored vector/residual, not interval proof.
Exact two-by-two Cauchy-Schwarz condition stated for declared K=B*F factorization.
Original quotient center selection 2(p+2q)+3ell+n=0 mod6; all color/weak reps
included even on pure hypercharge argument. Heat coefficient normalization is
declared, not inferred from cognition. Two links and two interfaces give power4.
Gaussian upper tail for d3^2 and C3 checked; union bound gives t4*s2+s4*t2.
SU2/U1 bounds and quotient restriction handled conservatively. Initial fixed
cutoffs underflowed a tail; attempt preserved. Adaptive exponent~150 gives
positive finite tail bound. Floating rounding not interval certified.
Two samples distinguish positive electric factor from positive combined kernel;
small tau sample positivity is not all-polynomial RP or physical parameter fit.
Shared temporal boundary fixed, Gauss closure absent. Therefore no original H_F,
strong-coupling or continuum no-go.643 endpoint and Fock closure both remain
required; positive Gram of original process already known, not rediscovered.
New report16 formulas;2 groups. Next673 tests actual boundary/Gauss closure,
not further thermal tuning. Original four-branch unified goal remains open.
'''
for name in ('research_note_672.md','joint_gauge_history_kernel.py','joint_gauge_history_kernel_results.json','unified_physics_condition_ledger_672.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round672_drafts/final_review.txt',review)
print('672 preparation completed; inherited evidence untouched')
