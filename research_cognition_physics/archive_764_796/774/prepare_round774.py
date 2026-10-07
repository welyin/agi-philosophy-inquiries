"""Freeze a conditional result and prepare guarded navigation publication."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent


def write(name, value):
    p = HERE/name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as stream:
        stream.write(value)


note = (HERE/'research_note_774.md').read_text('utf8')
write('round774_drafts/research_note_774_draft.md', note)
ledger = (HERE/'unified_physics_condition_ledger_773.md').read_text('utf8').splitlines()
ledger[0] = '# 联合条件总账：774来源首项与共同量子作用的条件性连接'
ledger[2] = '2026-10-04。接[773全账](unified_physics_condition_ledger_773.md)，回填[774报告](research_note_774.md)。[结果](joint_source_action_lift_results.json)、[核验](research_round_774_checks.json)。在明确N1—N3下，可同解一圈作用、771首阶来源和关系观测插入；N1—N3的原对象共同实现尚未签收。'
updates = {
    'C01': '774只给共同规范化N1—N3下的局部一圈连接；原相互作用物理正态与实际仪器仍待构造',
    'C11': '同一in-in适当顶点身份若实现，则反常首项为J_raw(Rc)，771局部接触可作为作用反项的一次jet；非局部/初态边界差不可略',
    'C19': '原关系观测可加形式源并与作用反常同解；理论插入不是实际测量或新物理参考',
    'C22': '773收缩在N1—N3下同时修复作用和观测源接触；原共同正规化前提、旧UV类与全局物理实现保持开放',
}
found = set()
for i, line in enumerate(ledger):
    for key, value in updates.items():
        if line.startswith('|'+key+' '):
            parts = line.split('|')
            parts[2] += '；'+value
            ledger[i] = '|'.join(parts)
            found.add(key)
assert found == set(updates)
write('unified_physics_condition_ledger_774.md', '\n'.join(ledger)+'''

## 774当前结论（优先于保留的773执行安排）

首圈共同规范化前提分为N1（同一Wick/因果/字段独立体系及场方程）、N2（局部量子作用原理与一致性）、N3（该体系的实际in-in一次导数等于原来源，保所有接触和边界）。原完整对象的N1—N3共同实现尚未签收，不把成熟框架的引用当作对象映射。

在这些前提下，在壳Slavnov身份把反常首项与原J_raw绑定，771来源Ward恰给773保首项递推的相容条件。显式加入ℏ阶局部作用后，可保771来源与原首阶约束响应；给关系观测加形式源，再以同一局部收缩联立多插入接触。于是不用为作用、来源、观测各加一项独立“可以修复”的假设。

该结果是局部正则参考片、背景依赖扩大系数类、模ℏ²及形式涨落/插入源展开的条件连接。未计算实际圈反常，未构造共同CTP生成泛函；N3不是已完成事实。未证明严格UV类别、背景独立、相互作用正态、实际仪器或Q连续映射。新增有限量子作用可改变更高关联，不称无物理影响。

下一项：[775](round775_drafts/STATUS.md)核原共同因果对象和初态/接触的真实映射，优先验证N3并核N1/N2函数类别。不要重复773同调、774制造反常的代数校准或旧空间结果。
''')
write('round775_drafts/STATUS.md', '''# 第775轮入口：共同规范化前提的原对象映射

接[774条件连接](../research_note_774.md)与[774全账](../unified_physics_condition_ledger_774.md)。

1. 774只证明N1—N3下的共同作用/来源/观测修复，不代表原相互作用QME已经完成。
2. 优先核N3：保772原准备的因果/闭时路径一次适当顶点，与770—771来源的局部接触和初态边界，是否真属同一规则。不可用in-out真空行列式替代。
3. 同时核N1/N2：原111辅助混合系统、ghost/反场/frame与费米，能否在773扩大背景jet类中满足因果、字段独立、基本方程及局部QAP；明确是否牺牲了旧幂计数或背景独立。
4. 若对象映射成立，可直接复用774递推，不需逐项计算全部反常才能证明局部存在；若映射失败，记录具体冲突、边界项及替代路线。不要把未证明当作反例。
5. 首阶来源与实际初值保持共同；更高有限量子项不是自动无物理影响。实际相互作用正态/仪器、记录资源及Q连续仍开放。
6. 382—386、425、522—523及604、649/699按原范围保留。无目标编辑、新任务、定时或图像检查。
''')
write('round774_drafts/scope_and_dedup_review.md', '''# 774范围与去重审查

- 770外腿计数、771来源Ward、773保首项引理均直接复用，未计作新证明。
- 新增连接是在同一规范化N1—N3下，将真实来源首项、作用修复和关系观测源的修复绑定；显式量子作用与普遍时间序重定义分开。
- N1—N3的原对象共同实现仍未签收。报告没有把条件定理改写为实际完整QME、正态或仪器。
- 原圈反常未计算。三组代码为BV投影、非线性含反场的精确诊断、有限Gaussian阶次诊断；不冒充原连续实验。
- 新增量子有限作用是明示分支；保原一圈来源并不保证全部高阶关联不变。
- 775核共同因果对象和初态/局部接触，不重复本文递推。原认知动机、旧空间成果和限定反例不扩张。
- 主代理审查；未使用独立代理或图像检查。
''')
audit = dict(
    primary_agent_review=True, original_loop_anomaly_computed=False,
    conditional_one_loop_lift_proven=True,
    original_N1_N2_N3_joint_realization_proven=False,
    original_interacting_QME_unconditionally_proven=False,
    source_projection_uses_actual_in_in_mapping_premise=True,
    initial_boundary_identity_not_assumed_globally=True,
    explicit_quantum_action_branch=True, higher_finite_terms_physical=True,
    sources=[
        dict(url='https://arxiv.org/pdf/1803.10235', locations=['Section 2, equation (25)', 'Theorems 3, 10-12', 'Qualification after equation (219)'],
             use='Local anomaly framework and distinction between common field-independent normalization and observable contact repairs.',
             limitation='Does not by itself identify the original source/initial state or validate all enlarged coefficient-class axioms.'),
        dict(url='https://arxiv.org/abs/0705.3160', locations=['Abstract; local quantum action principle'],
             use='Primary statement of the local anomaly/finite-renormalization framework.',
             limitation='Full PDF unavailable in this fetch; no uninspected theorem is claimed to apply to the full original gauge-gravity model.'),
        dict(url='https://arxiv.org/html/1306.1058v5', locations=['Relational observables', 'Section 3.4, equations (64)-(65)'],
             use='Comparison for relational insertions and local master identities.',
             limitation='Pure-gravity cohomology is not substituted for the original material theory.'),
        dict(url='https://arxiv.org/html/2604.26941v1', locations=['Sections 3.4-3.5, equations (3.46)-(3.49)'],
             use='Initial-state/boundary terms in closed-time-path BRST identities.',
             limitation='Not a proof of the original dynamical gravity/material process or interacting prepared state.'),
    ])
write('round774_drafts/literature_scope_audit.json', json.dumps(audit, ensure_ascii=False, indent=2)+'\n')
main = ['research_note_774.md', 'joint_source_action_lift.py',
        'joint_source_action_lift_results.json', 'unified_physics_condition_ledger_774.md']
review = '774主代理最终范围审查：条件性连接成立；原N1—N3尚未共同签收。无独立代理、无图像检查。\n'
review += '\n'.join(p+' '+hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in main)+'\n'
write('round774_drafts/final_review.txt', review)

mapping = {'773': '774', '772': '773', '774': '775',
           '3505': '3508', '3502': '3505', '1631': '1634',
           '3789': '3806', '3806': '3815',
           'joint_material_local_brst': 'joint_source_action_lift'}
regex = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda txt: regex.sub(lambda match: mapping[match[0]], txt)
verify = convert((HERE/'verify_round773.py').read_text('utf8'))
start, stop = verify.index('    extra = ('), verify.index('    new = ')
verify = verify[:start]+'''    extra = ('unified_physics_condition_ledger_774.md',
             'round774_drafts/research_note_774_draft.md',
             'round774_drafts/final_review.txt',
             'round774_drafts/literature_scope_audit.json',
             'round774_drafts/scope_and_dedup_review.md',
             'round775_drafts/STATUS.md')
'''+verify[stop:]
verify = verify.replace("checks['display_formulas'] == 14", "checks['display_formulas'] == 15")
start, stop = verify.index('                original_physical_free_state_preserved='), verify.index('                primary_code_and_note_review_completed=')
verify = verify[:start]+'''                original_physical_free_state_preserved=True,
                conditional_source_preserving_action_lift_proven=True,
                conditional_joint_observable_source_lift_proven=True,
                N1_N2_N3_joint_original_realization_proven=False,
                original_loop_anomaly_computed=False,
                original_interacting_QME_unconditionally_proven=False,
                strict_global_polynomial_counterterm_result=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
'''+verify[stop:]
write('verify_round774.py', verify)
publish = convert((HERE/'publish_round773.py').read_text('utf8'))
summary = '**第774轮完成（来源与量子作用的条件性连接）：** [保首项及共同插入]({p}research_note_774.md)在明确N1—N3下，同一局部一圈作用修复可保771来源并共同提升关系观测；原N1—N3共同实现、相互作用正态及仪器仍未签收。三组、十五式通过，最新774／3508，1634份编号科学文件、3815份保护证据。[核验]({p}research_round_774_checks.json)、[条件账]({p}unified_physics_condition_ledger_774.md)。'
order = '**当前执行顺序（774后，优先于下方历史安排）：** 接[775共同规范化的原对象映射]({p}round775_drafts/STATUS.md)，优先核同一实际来源、因果处方及初态/接触边界；不重复本轮条件递推或773同调。[范围审计]({p}round774_drafts/scope_and_dedup_review.md)。原物理准备、旧空间与目标保持。'
start, stop = publish.index('summary = '), publish.index('planned = ')
publish = publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish = publish.replace('第774轮完成（原物质参考与局部BRST）', '第774轮完成（来源与量子作用的条件性连接）')
publish = publish.replace('## 419.', '## 420.').replace('## 324.', '## 325.')
publish = publish.replace('原物质参考与局部规范左逆', '来源与共同量子作用的条件性连接')
publish = publish.replace('旧空间合同保持，局部参考与反项类别', '旧空间合同保持，条件性一圈共同插入')
publish = publish.replace('局部规范左逆与形式修复', '保来源首项与共同插入')
write('publish_round774.py', publish)
post = convert((HERE/'postcheck_round773.py').read_text('utf8'))
post = post.replace('第774轮完成（原物质参考与局部BRST）', '第774轮完成（来源与量子作用的条件性连接）')
write('postcheck_round774.py', post)
print('Prepared conditional round 774; no old evidence modified.')
