"""Prepare 770 evidence and navigation scripts without replacing history."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(name, content):
    path = HERE/name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)


lines = (HERE/'unified_physics_condition_ledger_769.md').read_text('utf8').splitlines()
lines[0] = '# 联合条件总账：770同态BRST减除、完整来源候选及局部Ward余项'
lines[2] = '2026-10-04。接[769全账](unified_physics_condition_ledger_769.md)，回填[770报告](research_note_770.md)。[结果](joint_local_brst_subtraction_results.json)、[核验](research_round_770_checks.json)。自由BRST减除及背景插入同一性已接通；来源有限不等于守恒，局部Ward修复仍开放。'
updates = {
    'C01': '770局部BRST减除保767/768同一实际物理态，减除核不是额外量子准备',
    'C11': '770绝对来源的剩余Ward写为局部方程余项超迹；未证明允许局部修复，不能先接绝对约束反馈',
    'C19': '770原111字段扩展允许精确自由BRST减除；全相互作用Ward不由此自动成立',
    'C22': '770完整有限tadpole候选与背景二次插入相同，自由规范固定插入抵消；联合守恒须另核局部余项',
}
seen = set()
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith('|'+key+' '):
            parts = line.split('|'); parts[2] += '；'+value
            lines[i] = '|'.join(parts); seen.add(key)
assert seen == set(updates)
write('unified_physics_condition_ledger_770.md', '\n'.join(lines)+'''

## 770同一BRST减除及来源候选

对原D1/D0取同尺度局部Hadamard奇异核，按768辅助块结构组装。其自由BRST身份精确成立，不要求两个减除核之间精确交织；该交织差是光滑的。完整方程余项由r1、r0及q三块给出，并不因为自由BRST成立而为零。

同一原作用的实际三阶tadpole含完整玻色混合、ghost顶点及原费米来源。背景二次作用变分还含规范固定插入；在这份减除下其b-v、ghost项共同抵消，b-b块为零。因此当前自由二次插入接口已关闭，不必为此先求全阶背景独立。

已构造有限来源候选，而其在壳联合Ward仍是原减除核的局部方程余项超迹。735已守恒费米来源沿原范围加入。需要同一允许局部有限修复抵消剩余项；没有计算出非零不可消去反常，也没有证明全来源守恒。不能以非局部投影或约束来源右逆代替局部重整化。

## 下一项：771局部余项与实际规范化

直接计算完整D1/D0/K所需有限重合jet，核成熟变分Hadamard/参数子处方与当前完整来源的映射，再处理联合Ward修复。不继续有限矩阵精度或独立纯引力真空模型。保持原Q、共同实际记录、内部资源及所有冻结空间条件。
''')
write('round770_drafts/research_note_770_draft.md', (HERE/'research_note_770.md').read_text('utf8'))
write('round771_drafts/STATUS.md', '''# 第771轮入口：局部方程余项与同一绝对来源的Ward修复

接[770](../research_note_770.md)、[全账](../unified_physics_condition_ledger_770.md)。

1. 同一767/768自由态可取精确BRST相容的局部减除；三阶tadpole来源与背景二次插入的自由规范固定差已经相消，不再重做矩阵块或G变分。
2. 尚缺的是770式(13)的实际局部Ward余项以及式(14)的共同局部修复。r1、r0及q来自原完整耦合D1/D0/K的Hadamard递推；不能替换成真空引力加独立标量。
3. 先核可用成熟变分减除/参数子构造，明确密度Hessian、真实in-in核、ghost与原F>0系数。背景协变不自动等于未规范固定的量子Ward；但770已经消去了当前自由规范固定插入，不能把它再次列为缺口。
4. 735费米处方直接保留，不重算628荷表。固定背景首阶来源不需要先证明全阶背景独立，也不需要无限背景幅度Taylor级数收敛。局部反项函数类别及其原始元仍须交代。
5. 证明剩余局部项可修复后再接754绝对首阶约束；未证明时保有限候选的准确范围。禁止用状态依赖正常序置零或非局部投影伪装局部修复。
6. 原Q连续过程、自主实际记录、资源和有限量子强度反馈仍开放；旧空间382—386/425/522—523、604及649/699范围保持。目标不改，无新任务、调度或图像。
''')
write('round770_drafts/scope_and_dedup_review.md', '''# 770范围与去重复核

- 上一目标轮完成769及770工作入口，是有效进展。本轮读导航、最新报告/结果，未见活跃Python进程。
- 直接使用734实际历史来源、735费米局部规范化、765完整双曲复形、768强BRST核及769变量字典；一般Hadamard和BRST方法不是新理论。
- 本轮增量是原同态扩展的局部减除可以精确保自由BRST，即便局部核不精确交织；具体r1/r0/q给其完整方程余项。
- H为Wightman型奇异核；方程余项无Feynman逆的delta。辅助块精确BRST是代数事实；物理全态、正性和因果性沿旧成果保持。
- 密度Hessian与原纤维算符分开。普通规范参数固定后求背景导数，协变参数字典的背景项不丢。
- 实际三阶来源含ghost。背景二次插入额外的规范固定项为自由BRST精确插入，在本处方中抵消；没有升级成全相互作用量子主方程。
- 余项超迹使用光滑W-H与紧支测试进行分部积分，不对发散迹随意循环。有限来源不自动守恒，局部修复仍是实际待验收问题。
- 数值仅为完整111块及原H5势/SU2的零维校准；声明辅助G的小混合，不当原753实际时空规范或连续圈结果。早期开发校准G过近奇异，最终显式取1/4混合并检查谱；不是调整物理耦合或隐藏理论反例。
- 770工作入口的有限价数枚举保持冻结，不计入本轮两个科学测试组。
- 主代理逐式/代码审查，没有独立审稿、图像检查或目标改动。
''')
write('round770_drafts/literature_scope_audit.json', json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/1708.00444', sections='I, III and IV',
             use='Local Hadamard kernels, free Ward identities and parameter dependence.',
             boundary='Its vacuum tensor background is not used as a theorem for the original mixed Hessian.'),
        dict(url='https://arxiv.org/html/1804.07640', sections='3.1.2-3.1.3 and 4',
             use='Background/field split and BRST-exact gauge-fixing insertion.',
             boundary='Only the current free quadratic insertion is closed here; interacting anomaly insertion remains unproved.'),
        dict(url='https://arxiv.org/html/2403.16806', sections='I.1-II',
             use='One-loop radiation currents, ghost inclusion, in-in/in-out distinction.',
             boundary='No transfer of a pure tensor/FLRW calculation to the full original mixed model; no absorption of arbitrary state dependence into local counterterms.')],
    exact_free_BRST_subtraction_proven=True,
    local_equation_remainders_identified=True,
    free_background_gauge_fixing_insertion_cancelled=True,
    all_sector_finite_source_candidate=True,
    all_sector_conserved_source_proven=False,
    interacting_QME_proven=False,
    local_repair_primitive_constructed=False,
    original_state_unchanged=True,
), ensure_ascii=False, indent=2)+'\n')

main = ['research_note_770.md', 'joint_local_brst_subtraction.py',
        'joint_local_brst_subtraction_results.json', 'unified_physics_condition_ledger_770.md']
review = '770主代理终审：十六式、两组复算及代码核验；不是独立审稿。\n'
review += '核：H块只用自由BRST；不要求H1K=KH0；r1/r0/q完整；真实三阶来源含ghost；密度转换与两腿方程；有限候选不当守恒来源。\n'
review += '\n'.join(p+' '+hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in main)+'\n'
write('round770_drafts/final_review.txt', review)

mapping = {'769':'770','768':'769','3494':'3496','3491':'3494','1619':'1622',
           '3742':'3757','3726':'3742','joint_quantum_frame_transport':'joint_local_brst_subtraction'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda s: pattern.sub(lambda m: mapping[m[0]], s)
verify = convert((HERE/'verify_round769.py').read_text('utf8'))
start, stop = verify.index('    extra = ('), verify.index('    new = ')
verify = verify[:start]+'''    extra = ('unified_physics_condition_ledger_770.md',
             'round770_drafts/research_note_770_draft.md',
             'round770_drafts/final_review.txt',
             'round770_drafts/literature_scope_audit.json',
             'round770_drafts/scope_and_dedup_review.md',
             'round771_drafts/STATUS.md',
             'round770_drafts/research_note_770_working.md',
             'round770_drafts/finite_valence_entry.py',
             'round770_drafts/finite_valence_entry_results.json',
             'round770_drafts/entry_scope_review.json',
             'round770_drafts/finite_valence_entry_checks.json',
             'round770_drafts/entry_postpublication_checks.json')
    entry = core.read(HERE/'round770_drafts/finite_valence_entry_checks.json')
    for p, digest in entry['working_artifact_hashes'].items():
        assert core.digest(HERE/p) == digest, p
'''+verify[stop:]
verify = verify.replace('== (3, 0, 0)', '== (2, 0, 0)').replace('fresh_tests=dict(run=3,', 'fresh_tests=dict(run=2,')
start, stop = verify.index('                original_free_state_frame_transport_proven='), verify.index('                primary_code_and_note_review_completed=')
verify = verify[:start]+'''                exact_free_BRST_subtraction_proven=True,
                finite_same_state_complete_source_candidate=True,
                free_gauge_fixing_insertion_cancelled=True,
                Ward_defect_local_equation_remainder_expression=True,
                all_sector_Ward_repair_proven=False,
                absolute_conserved_backreaction_proven=False,
                interacting_QME_proven=False,
                original_Q_to_E_process_equivalence_proven=False,
'''+verify[stop:]
write('verify_round770.py', verify)
publish = convert((HERE/'publish_round769.py').read_text('utf8'))
start, stop = publish.index('summary = '), publish.index('planned = ')
summary = '**第770轮完成（同态局部减除与来源候选）：** [BRST减除及局部Ward余项]({p}research_note_770.md)构造保持自由BRST的共同减除，原完整来源与背景二次插入的规范固定差抵消；剩余守恒缺陷明确为局部方程余项，尚未证明可修复。两组、十六式通过，最新770／3496，1622份编号科学文件、3757份保护证据。[核验]({p}research_round_770_checks.json)、[条件账]({p}unified_physics_condition_ledger_770.md)。'
order = '**当前执行顺序（770后，优先于下方历史安排）：** 接[771局部余项与共同规范化]({p}round771_drafts/STATUS.md)，计算原完整算符的局部Ward并核有限修复；不重复自由规范固定插入或矩阵精度。[范围审计]({p}round770_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
publish = publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish = publish.replace('第770轮完成（同态换帧与来源接触）', '第770轮完成（同态局部减除与来源候选）')
publish = publish.replace('## 415.', '## 416.').replace('## 320.', '## 321.')
publish = publish.replace('同态换帧、平均场接触与相对来源', '同态BRST减除与局部Ward来源')
publish = publish.replace('旧空间合同保持，同态换帧与绝对来源', '旧空间合同保持，同态减除与局部规范化')
publish = publish.replace('原两帧的共同量子字典', 'BRST减除及局部Ward余项')
publish = publish.replace('next_round=770', 'next_round=771')
write('publish_round770.py', publish)
post = convert((HERE/'postcheck_round769.py').read_text('utf8'))
post = post.replace('range(584, 770)', 'range(584, 771)')
post = post.replace('第770轮完成（同态换帧与来源接触）', '第770轮完成（同态局部减除与来源候选）')
write('postcheck_round770.py', post)
print('Prepared 770 evidence and publication helpers.')
