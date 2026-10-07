"""Prepare 767 scoped evidence and publication helpers with exclusive writes."""
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def write(name, content):
    p=HERE/name
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:
        f.write(content)


lines=(HERE/'unified_physics_condition_ledger_766.md').read_text('utf8').splitlines()
lines[0]='# 联合条件总账：767同一完整线性物理商的正Hadamard态'
lines[2]='2026-10-04。接[766全账](unified_physics_condition_ledger_766.md)，回填[767报告](research_note_767.md)。[结果](joint_physical_hadamard_positivity_results.json)、[核验](research_round_767_checks.json)。固定原背景上自由物理正性已接通；下方各旧轮次的未完成语句保留历史语义，以本条及末节为当前状态。'
updates={
    'C01':'767原完整线性物理商存在正Hadamard准自由态；不替代原Q共同过程或完整相互作用',
    'C04':'767只加平滑核，保765完整因果传播与766精确规范/CCR',
    'C09':'767同一无稳定子背景上实现正物理态；原参考的实际量子读取未完成',
    'C11':'767物理约束商正性成立，完整非线性量子约束和重整化来源仍独立',
    'C19':'767正确频率椭圆块与无限秩平滑余项分别处理，共同协方差修正保规范/实性/短距离',
    'C22':'767自由Hadamard准备存在；735只一费米圈，全部门局部Ward、实际来源及原记录映射仍开放'}
seen=set()
for i,line in enumerate(lines):
    for key,value in updates.items():
        if line.startswith('|'+key+' '):
            parts=line.split('|');parts[2]+='；'+value
            lines[i]='|'.join(parts);seen.add(key)
assert seen==set(updates)
write('unified_physics_condition_ledger_767.md','\n'.join(lines)+"""

## 767同一完整线性物理商的Hadamard正性

在固定原753无连续稳定子、紧初片、正F及完整二导数作用下，E微扰量子化自由部门存在正的物理Hadamard准自由态。新构造以两个正椭圆Gram逆实现精确物理代表空间，用766全部耦合频率投影在该空间分频；物理主符号的正性来自原正动能和横向/横向无迹规范商。

正确频率块有严格正的椭圆主部，因此只有有限低谱需有限秩修正。剩余交叉和错误频率项是一般平滑核，可无限秩；通过椭圆本征基的矩阵元行和给一份正平滑共同上界。把同一个协方差加到正负二点核，精确保完整CCR、规范和实结构，且不改变短距离波前。这补767入口，而不是把另一份任意Gaussian态与766拼接。

态一般混合且不唯一，修正会改变真实平滑关联及来源值，未证明原指定态或实际仪器可有限成本制备。原Q到E连续映射、共同真实记录、全部门Ward/二阶来源和相互作用仍未完成。没有推导原输入的维数、群、作用、正F或普遍认知原则。

## 下一项：共同自由态到同阶全部门来源

接768，先明确玻色、费米、规范固定/ghost的二次来源与同一物理代数的映射，复用734—735完整参考及接触项，不以Hadamard状态存在直接宣布局部Ward。现有平方根尺度下二次玻色与费米来源同阶；按原完整作用和二阶约束核验。旧空间、604、649/699及统一目标保持。
""")
write('round767_drafts/research_note_767_draft.md',(HERE/'research_note_767.md').read_text('utf8'))
write('round768_drafts/STATUS.md',"""# 第768轮入口：共同自由态与同阶完整来源

接[767](../research_note_767.md)、[当前全账](../unified_physics_condition_ledger_767.md)。

1. 原753无连续稳定子背景上的完整线性玻色场已有同一正物理Hadamard准自由态；正性存在任务关闭，不再优化辅助谱阈值或噪声常数。
2. 764明确B=B*+sqrt(epsilon)v+epsilon w，二阶约束来源包含1/2 D²C[v,v]及费米来源；Hadamard只解决局部乘积可减除的入口，不自动解决约束/Ward。
3. 回查734—735的完整过去参考、局部接触项、一费米圈处方；735没有量子玻色场或引力圈，不能直接扩大其签收范围。
4. 新增对象要有同一原作用的玻色Hessian、费米算符、ghost/规范固定及BRST身份；先核准确成熟定理前提与当前完整非零物质背景的映射，再作联合来源结论。
5. 必须区分物理规范商上的态、完整未约束二点扩展和局部复合算符。单独“引力涨落的应力张量”未自动成为规范不变可观测。
6. 仍须共同连接原Q、实际原记录/装置及受控反作用。不要把自由态存在或预置Einstein作用的验证称为认知生成全部物理。
7. 旧空间382—386、425、522—523与604、649/699边界保持；目标不改，无新应用任务/调度/图像。
""")
write('round767_drafts/scope_and_dedup_review.md',"""# 767范围与去重复核

- 先读取当前七项导航中的主文件、765—766报告与结果、753无稳定子、734—735和764阶次；无运行中的Python任务。
- 766已核验发布，3679等旧计数不重新使用；本轮只新增三组当前公式校准，历史内容按哈希保留。
- R的两个像H正交依赖K*QK=0。T、A、M、其逆均有明确有限阶；K核零继承753，不另做秩扫描。
- R与c全阶交换至平滑，不只是主符号对易。C是正收缩且幂等至平滑；Riesz谱投影存在于真实投影Hilbert空间，平滑差由椭圆逆与缺陷因子给出。
- 物理正主符号来自既定正作用及规范商；不能以辅助H正替代真实Q正。
- 只在严格正椭圆频率块上使用有限低谱修复；错误频率和平滑交叉可无限秩，另外使用共同对角行和上界。
- 平滑上界证明逐项给快速衰减、谱可求和与Hermitian行和不等式；不假定任意紧扰动可有限秩处理。
- 噪声同加λ±，湮灭规范且实对称化，保完整CCR及交织。不是指定原态不变、纯态或物理免费操作。
- 三组有限矩阵只核对象和代数；真实全变系数逆、频率投影及Hadamard存在由解析论证承担。
- 不新增认知公理；原维数、群、作用、正F、紧初片、在壳背景和微扰量子化明确继承。总目标未完成。
- 主代理逐式复核；未做独立代理审查。外部独立数学审查仍有价值，不能把机器精度当证明可靠性的替代。
""")
write('round767_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/2311.11043',
             checked='Sections3.1-3.2: general indefinite normally hyperbolic Hadamard projectors. Exact projection, charge adjoint and wavefront properties; not physical positivity.',
             mapping='Inherited full coupled D1 and normalization from765-766. No transfer of the later vacuum Einstein positivity theorem.'),
        dict(url='https://arxiv.org/html/1403.7153',
             checked='Theorem3.17 and Remark3.18 distinguish exact gauge/CCR from physical positivity.',
             mapping='Uses766 exact candidates; this round adds the same positive smoothing covariance on the physical quotient.')],
    project_increment='Explicit physical orthogonal representative, exact frequency splitting modulo smoothing, finite correct-branch spectral repair, and infinite-rank common smooth majorant preserving gauge/CCR/reality.',
    no_claimed_new_discovery_of_standard_pseudodifferential_calculus=True,
    inherited_fixed_background_free_state_existence_proven=True,
    full_interacting_quantum_gravity_or_original_Q_mapping_proven=False),ensure_ascii=False,indent=2)+'\n')
files=('research_note_767.md','joint_physical_hadamard_positivity.py',
       'joint_physical_hadamard_positivity_results.json','unified_physics_condition_ledger_767.md')
write('round767_drafts/final_review.txt',
      'Primary review only. Checked fixed compact irreducible background, full coupled symbols and '
      'finite-order operators. T commutes exactly with c; K and TK have orthogonal closed ranges. '
      'Physical R is a Psi0 orthogonal projection, and R|N=Pi|N. Positive physical principal symbols '
      'hold after the actual gravity/YM/scalar gauge quotient, not merely because H is positive. '
      'Projected frequency Riesz construction is exact modulo a smoothing defect. Only strictly '
      'positive elliptic branch blocks use finite low-spectrum repair; all remaining smoothing '
      'tails use a proved common infinite-rank matrix majorant. Same noise preserves CCR, '
      'annihilates gauge, and real symmetrization preserves both signs of positivity. This yields '
      'existence of a free physical Hadamard state, not a unique state, fixed preparation, '
      'finite-resource instrument, all-sector Ward normalization or the original Q continuum limit.\n'
      +'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n')

mapping={'766':'767','765':'766','3691':'3705','3705':'3714',
         '3482':'3485','3485':'3488','1610':'1613',
         'joint_microlocal_gauge_projection':'joint_physical_hadamard_positivity'}
pattern=re.compile('|'.join(map(re.escape,sorted(mapping,key=len,reverse=True))))
def convert(s):
    return pattern.sub(lambda m:mapping[m[0]],s)

verify=convert((HERE/'verify_round766.py').read_text('utf8'))
a,b=verify.index('    extra = ('),verify.index('    new = ')
verify=verify[:a]+"""    extra = ('unified_physics_condition_ledger_767.md',
             'round767_drafts/research_note_767_draft.md',
             'round767_drafts/final_review.txt',
             'round767_drafts/literature_scope_audit.json',
             'round767_drafts/scope_and_dedup_review.md',
             'round768_drafts/STATUS.md')
"""+verify[b:]
verify=verify.replace("checks['display_formulas'] == 16","checks['display_formulas'] == 18")
a,b=verify.index('                actual_coupled_Cauchy_gauge_ellipticity_and_injectivity_proven='),verify.index('                primary_code_and_note_review_completed=')
verify=verify[:a]+"""                exact_physical_representative_constructed=True,
                common_positive_smoothing_majorant_proven=True,
                physical_positive_Hadamard_state_proven=True,
                physical_state_scope='Fixed inherited background; complete linear bosonic physical quotient only.',
                numerical_checks_only_finite_algebraic_calibrations=True,
                preparation_uniqueness_or_finite_resource_implementation_proven=False,
                all_sector_quantum_Ward_completed=False,
                original_Q_to_E_process_equivalence_proven=False,
"""+verify[b:]
write('verify_round767.py',verify)

publish=convert((HERE/'publish_round766.py').read_text('utf8'))
a,b=publish.index('summary = '),publish.index('planned = ')
publish=publish[:a]+"""summary = '**第767轮完成（同一自由物理态的正性）：** [共同平滑修正与Hadamard正态]({p}research_note_767.md)在原紧初片无稳定子背景上，对完整耦合线性玻色场共同满足正性、规范、CCR和短距离条件。态选择不唯一，未证原实际制备、全来源或相互作用。三组、十八式通过，最新767／3488，1613份编号科学文件、3714份保护证据。[核验]({p}research_round_767_checks.json)、[条件账]({p}unified_physics_condition_ledger_767.md)。'
order = '**当前执行顺序（767后，优先于下方历史安排）：** 接[768同阶全部门来源与共同Ward]({p}round768_drafts/STATUS.md)，复用734—735，核同一原作用、物理商及玻色/费米/ghost来源；不继续态的辅助常数。[范围审计]({p}round767_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
"""+publish[b:]
publish=publish.replace('第767轮完成（精确规范与短距离共同候选）','第767轮完成（同一自由物理态的正性）')
publish=publish.replace('## 412.','## 413.').replace('## 317.','## 318.')
publish=publish.replace('频率相容的规范投影与短距离共同候选','同一完整自由物理商的Hadamard正性')
publish=publish.replace('旧空间合同保持，精确规范与物理正性','旧空间合同保持，自由物理态与完整来源')
publish=publish.replace('频率相容的规范投影','共同平滑修正与Hadamard正态')
publish=publish.replace('next_round=767','next_round=768')
write('publish_round767.py',publish)
post=convert((HERE/'postcheck_round766.py').read_text('utf8'))
post=post.replace('range(584, 767)','range(584, 768)')
post=post.replace('第767轮完成（精确规范与短距离共同候选）','第767轮完成（同一自由物理态的正性）')
write('postcheck_round767.py',post)
print('767 report, ledger, audits and publication helpers prepared.')
