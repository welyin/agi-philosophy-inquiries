"""Prepare 768 evidence without overwriting frozen work."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:
        f.write(text)


lines=(HERE/'unified_physics_condition_ledger_767.md').read_text('utf8').splitlines()
lines[0]='# 联合条件总账：768共同BRST扩展与完整玻色相对来源'
lines[2]='2026-10-04。接[767全账](unified_physics_condition_ledger_767.md)，回填[768报告](research_note_768.md)。[结果](joint_brst_relative_source_results.json)、[核验](research_round_768_checks.json)。同一自由态的辅助扩展及指定平滑变化的玻色来源差已接通；绝对全部门来源仍开放。'
updates={
    'C01':'768原物理Hadamard态有保同一物理商的精确线性BRST二点扩展，未约束辅助代数不要求全正',
    'C04':'768完整扩展Green和二点核保同一因果结构及精确分次对易关系',
    'C09':'768源与原共同参考仍属于753同一模型，未完成实际量子读取与内部装置',
    'C11':'768平滑物理二点差为P双解，来源满足联合Noether并接原全来源约束右逆；非线性量子约束未完成',
    'C19':'768物理协方差差可只改字段块，辅助/ghost准备保持；不把一般态差当独立正噪声',
    'C22':'768完整玻色相对二次来源有限并有形式首阶响应；绝对全部门Ward、有限ε自洽和原Q映射未完成'}
seen=set()
for i,line in enumerate(lines):
    for key,value in updates.items():
        if line.startswith('|'+key+' '):
            parts=line.split('|');parts[2]+='；'+value
            lines[i]='|'.join(parts);seen.add(key)
assert seen==set(updates)
write('unified_physics_condition_ledger_768.md','\n'.join(lines)+"""

## 768同一自由物理态的BRST扩展与相对来源

翻转辅助规范参数束的表示配对后，765原P、K、D1、D0精确适配成熟线性BRST扩展。767的完整Cauchy身份输送为时空二点交织；完整111字段块含全部16规范参数、辅助乘子、ghost和anti-ghost，gh=0物理态与767相同。物理正性不推广到未约束ghost代数。

对同一准备族中实相容、平滑、湮灭规范的物理协方差变化，两端均须为正态；辅助二点块保持不变。变化核满足原P双解和subsidiary条件，故原完整玻色作用三阶变分给有限来源差，联合Noether由同一核严格成立。来源代表不单独规范不变，二阶背景修正随之共同变换。

原731/754全来源右逆与732短时相对发展可直接接该同态来源，给形式epsilon阶的平均背景响应。这里没有新增任意诊断源、补偿物种或参考场；也没有证明两个绝对量子分支存在、有限epsilon非线性闭合或原Q连续映射。费米来源按735准确范围独立保留，不删除绝对高阶有效项或有限反项。

## 下一项：绝对全部门局部来源与量子Ward

接769。成熟BRST线性扩展已映射，不再重复分块代数。先核同一原作用的完整一圈/局部有效理论规范化，区分背景Ward与量子BV主方程、绝对处方与平滑态差；回用735一费米圈、580/600高阶项及628反常范围。参考/记录、原Q映射和所有独立输入保留，总目标不变。
""")
write('round768_drafts/research_note_768_draft.md',(HERE/'research_note_768.md').read_text('utf8'))
write('round769_drafts/STATUS.md',"""# 第769轮入口：同阶绝对来源及全部门量子Ward

接[768](../research_note_768.md)、[全账](../unified_physics_condition_ledger_768.md)。

1. 原完整自由玻色物理态已经有精确BRST二点扩展；指定平滑物理态差的完整来源及形式首阶响应已接通，不继续投影/ghost分块扫描。
2. 下一项是同一原作用下绝对的玻色、引力、ghost与费米来源。735只给一费米圈，不能直接扩为全部门或全圈。
3. 先核成熟局部EFT/BV规范化结果及真实前提。背景变换下协变，不自动等同于量子规范主方程；无局部反常也不自动固定全部有限反项。
4. 原H5非线性目标和Einstein作用应按明确的有效理论阶次处理；580/600已有高阶局部项不应重新宣布发现。不得把二导数原作用当未经证明可重整化的闭合集合。
5. 767固定在壳物理商不是任意离壳背景族；绝对背景变分需明确离壳BV或等价映射。准备和实际记录保持同一历史，不随每次变分重选真空。
6. 继续核来源、二阶约束、原Q连续匹配及实际装置的共同验收。旧空间382—386、425、522—523与604、649/699保持。目标未改，无新任务/调度/图像。
""")
write('round768_drafts/scope_and_dedup_review.md',"""# 768范围与去重复核

- 上一目标轮766—767及768入口改变了正式状态，为有效进展。本轮先读导航、最新笔记/结果及进程，无活跃Python实验。
- 回用732相对发展、734接触项、735一费米圈、754全来源右逆及756平滑差方法；不将成熟定理复述或普通Noether恒等式计为独立新轮。
- 项目增量是原全部线性几何/规范/H5物理态共同接精确BRST二点扩展，再给该同态族的完整玻色来源差及原受约束响应映射。
- 时空dagger、W翻号后的star与Cauchy电荷ddagger区分；形式与算符二点核的识别同步翻号。
- ghost负频块具有额外负号；完整分次CCR、场方程成立，不能只模规范成立。正性只在gh=0物理商。
- 平滑物理差湮灭K，故辅助二点差零、P双解。该强结论不施加到整份未约束Hadamard核。
- 来源从原作用第三变分定义，保体积、目标度量及全部混合；共用局部反项在差中消失，不声称绝对Ward。
- 二次来源依赖代表，与二阶背景修正一同变换；未把引力应力单独说成规范不变观测。
- 相对首阶Cauchy问题使用真实来源和旧约束右逆，不对非紧时间支撑源随意写负无穷退迟积分或乘开关。
- 原完整势的离壳局部Taylor jet只校准Noether接触项，不是连续量子来源数值。其他两组也仅有限代数校准。
- 原输入维数、群、作用、参考与量子化声明保持；自主装置、原Q、绝对全部门来源及相互作用不报完成。
- 主代理复核，无独立代理或图像审查。
""")
write('round768_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/1407.8079',
             checked='Section2.7, Lemma2.24, Prop3.4 and Remark3.3. General auxiliary BRST Green/two-point extension and exact versus modulo-gauge identities.',
             mapping='Original full P,K with W pairing sign reversal; D1=P+KKstar,D0=KstarK.767 covariance addition preserves the strong Cauchy intertwiner, which transports to spacetime.',
             not_imported='Specific Maxwell/YM physical-state existence proofs, interacting BRST master equation or all-sector source renormalization.')],
    inherited_scope=['732/754 constrained relative linear response','734 same-past contact terms',
                     '735 one-fermion-loop normalization','756 smooth fermionic source differences'],
    project_increment='Complete coupled bosonic physical state, auxiliary extension, smooth original-action third-variation source and inherited response refer to the same background and preparation family.',
    full_BRST_linear_extension_proven=True,
    smooth_relative_bosonic_source_joint_Noether_proven=True,
    absolute_all_sector_Ward_proven=False,
    finite_epsilon_semiclassical_closure_proven=False),ensure_ascii=False,indent=2)+'\n')
files=('research_note_768.md','joint_brst_relative_source.py',
       'joint_brst_relative_source_results.json','unified_physics_condition_ledger_768.md')
write('round768_drafts/final_review.txt',
      'Primary review only. Verified auxiliary sign flip, full coupled operator blocks, '
      'Green support reasoning and transport of exact Cauchy covariance intertwining. '
      'Negative ghost covariance signs and physical-only positivity are retained. '
      'The admissible smooth covariance difference changes only the original field block '
      'and is a genuine P bisolution with both subsidiary conditions. '
      'The original action third variation defines the complete relative bosonic source; '
      'Noether follows by contraction with this smooth bisolution. It does not establish '
      'absolute quantum Ward or make the isolated source representative gauge independent. '
      'The inherited constraint inverse and source Cauchy problem give only a formal '
      'first-order relative response, not existence of absolute quantum branches or finite-epsilon '
      'nonlinear control. Matrix calibrations cover signs and algebra, not continuum kernels.\n'
      +'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n')
mapping={'767':'768','766':'767','3705':'3714','3714':'3726',
         '3485':'3488','3488':'3491','1613':'1616',
         'joint_physical_hadamard_positivity':'joint_brst_relative_source'}
pat=re.compile('|'.join(map(re.escape,sorted(mapping,key=len,reverse=True))))
def conv(s):
    return pat.sub(lambda m:mapping[m[0]],s)

verify=conv((HERE/'verify_round767.py').read_text('utf8'))
a,b=verify.index('    extra = ('),verify.index('    new = ')
verify=verify[:a]+"""    extra = ('unified_physics_condition_ledger_768.md',
             'round768_drafts/research_note_768_draft.md',
             'round768_drafts/final_review.txt',
             'round768_drafts/literature_scope_audit.json',
             'round768_drafts/scope_and_dedup_review.md',
             'round769_drafts/STATUS.md',
             'round768_drafts/joint_source_entry.md',
             'round768_drafts/joint_source_entry_checks.json',
             'round768_drafts/entry_postpublication_checks.json')
"""+verify[b:]
a,b=verify.index('                exact_physical_representative_constructed='),verify.index('                primary_code_and_note_review_completed=')
verify=verify[:a]+"""                complete_linear_BRST_extension_constructed=True,
                original_physical_Hadamard_state_preserved=True,
                relative_full_bosonic_source_Noether_proven=True,
                relative_response_scope='Formal first-order coefficient using inherited constrained Cauchy response.',
                numerical_checks_only_finite_algebraic_calibrations=True,
                absolute_all_sector_Ward_completed=False,
                finite_epsilon_nonlinear_semiclassical_closure_proven=False,
                original_Q_to_E_process_equivalence_proven=False,
"""+verify[b:]
write('verify_round768.py',verify)
publish=conv((HERE/'publish_round767.py').read_text('utf8'))
a,b=publish.index('summary = '),publish.index('planned = ')
publish=publish[:a]+"""summary = '**第768轮完成（共同态与相对量子来源）：** [BRST扩展与完整玻色来源差]({p}research_note_768.md)保767同一物理态，平滑物理变化可保持辅助部门不变；原作用给有限联合守恒来源差，并接原形式首阶约束响应。绝对全部门Ward和有限强度闭合未证。三组、十八式通过，最新768／3491，1616份编号科学文件、3726份保护证据。[核验]({p}research_round_768_checks.json)、[条件账]({p}unified_physics_condition_ledger_768.md)。'
order = '**当前执行顺序（768后，优先于下方历史安排）：** 接[769绝对全部门来源与量子Ward]({p}round769_drafts/STATUS.md)，核同一原作用的局部EFT/BV规范化、离壳背景及有限项；不重复自由辅助分块。[范围审计]({p}round768_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
"""+publish[b:]
publish=publish.replace('第768轮完成（同一自由物理态的正性）','第768轮完成（共同态与相对量子来源）')
publish=publish.replace('## 413.','## 414.').replace('## 318.','## 319.')
publish=publish.replace('同一完整自由物理商的Hadamard正性','同一物理态的BRST扩展与完整相对来源')
publish=publish.replace('旧空间合同保持，自由物理态与完整来源','旧空间合同保持，同态来源与绝对规范化')
publish=publish.replace('共同平滑修正与Hadamard正态','BRST扩展与完整玻色来源差')
publish=publish.replace('next_round=768','next_round=769')
write('publish_round768.py',publish)
post=conv((HERE/'postcheck_round767.py').read_text('utf8'))
post=post.replace('range(584, 768)','range(584, 769)')
post=post.replace('第768轮完成（同一自由物理态的正性）','第768轮完成（共同态与相对量子来源）')
write('postcheck_round768.py',post)
print('768 evidence and publication helpers prepared.')
