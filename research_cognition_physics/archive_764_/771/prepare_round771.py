"""Freeze 771 evidence and prepare its guarded navigation publication."""
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


lines = (HERE/'unified_physics_condition_ledger_770.md').read_text('utf8').splitlines()
lines[0] = '# 联合条件总账：771同一完整来源的局部一圈Ward修复'
lines[2] = '2026-10-04。接[770全账](unified_physics_condition_ledger_770.md)，回填[771报告](research_note_771.md)。[结果](joint_local_ward_repair_results.json)、[核验](research_round_771_checks.json)。原同态完整一圈来源可作共同局部Ward修复；这不等于全相互作用QME、规范固定独立或有限强度自洽。'
updates = {
    'C01': '771局部来源修复不改变767/768实际态与735费米准备，所有同背景状态差保持',
    'C11': '771在原局部Wick来源类别内显式抵消完整一圈联合Ward缺陷，可接原形式首阶约束；全自洽响应尚待验收',
    'C19': '771完整辅助消元保局部接触，两个原耦合双曲算符的余项jet给所需修复；全部插入QME未证',
    'C22': '771完整有限首阶来源可共同守恒；同态记录、原Q连续映射和有限量子强度仍开放',
}
seen = set()
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith('|'+key+' '):
            parts = line.split('|'); parts[2] += '；'+value
            lines[i] = '|'.join(parts); seen.add(key)
assert seen == set(updates)
write('unified_physics_condition_ledger_771.md', '\n'.join(lines)+'''

## 771最新结论（优先于上述保留的770历史状态）

在原753在壳背景及同一767/768实际自由态上，完成辅助平方会留下局部接触；原密度配对变分也有局部接触。两项共同保留后，来源的Ward归结为原D1/D0两个完整耦合算符的交换余项。

原算符是保持非退化纤维配对的正常双曲算符，可用矩阵Hadamard输运与双点伴随对称。原余项的对角值为零，一阶右jet为−2ν∇[v1]；内部规范零阶插入消失，微分同胚缺陷由一个联合局部体积接触抵消。ghost权−1与玻色半权均保留。加入735同背景费米处方后，首阶完整来源联合守恒。

新增选择属于实、状态无关、局部协变的有限jet来源Wick重定义；没有要求它单独是某个不变作用的Euler导数。是否可扩展为所有相互作用时间序插入的一份规范化，仍待证明。不能把本轮来源守恒升级为QME、全部背景独立或规范固定独立。

原四维、群、物质谱、作用、正F、量子化分支及准备仍是明示输入；不是认知原则独立生成GR或标准模型。实际绝对来源的数值和物理有限常数尚未确定。

## 下一项：772绝对首阶来源、同态响应与实际过程

复用754的完整约束右逆，核771来源下的形式共同响应、初始修正及有限局部常数。把来源Ward、形式首阶解、量子态随背景改变、有限强度闭合与实际记录分别验收。优先推进共同过程接口，不继续局部系数或矩阵精度。旧空间、604、649/699、原Q和总目标保持。
''')
write('round771_drafts/research_note_771_draft.md', (HERE/'research_note_771.md').read_text('utf8'))
write('round772_drafts/STATUS.md', '''# 第772轮入口：同一绝对一圈来源与共同首阶响应

接[771](../research_note_771.md)、[最新总账](../unified_physics_condition_ledger_771.md)。

1. 原同态玻色/ghost来源已在明确局部Wick插入类别内完成Ward修复；735费米来源相加，得到原在壳背景的完整守恒首阶来源。不再把局部Ward修复列为缺失，也不升级为全相互作用QME。
2. 复用754全颜色/引力约束右逆与732形式响应，核绝对来源、背景初始修正及同一历史准备。不能仅重述“守恒源可响应”便计新轮次。
3. 优先核771处方的剩余有限常数、规范固定/变量依赖及原实际记录对应。区分自由态在原背景、态随背景的响应、形式首阶均值和有限强度自洽；新困难不重新否定已证窄结果。
4. 回用734同一过去的退迟比较、769字段字典与756同态记录/噪声。一般耦合玻色态的跨背景输送不得直接套费米定理；明确新增连接。
5. 原Q连续映射、自主实际仪器、内部资源与全相互作用共同规范化仍开放。四维、群/表示/参数、原作用及F>0仍为输入。空间382—386/425/522—523和604/649/699范围保持。
6. 先对接成熟响应/背景比较结果；无新任务、调度或图像，目标不改。若本轮只能做旧结果的直接应用，不虚增正式轮次，保存工作报告并推进真正新增连接。
''')
write('round771_drafts/scope_and_dedup_review.md', '''# 771范围与去重复复核

- 上一目标轮完成770并保存771入口，判为有效进展；本次先读导航、报告、结果与当前进程，未见运行Python。
- 734、735、765、768—770直接继承；工作入口密度接触不再计作独立科学轮次。
- 本轮增量为完整辅助消元接触＋原两个耦合算符的真实余项jet＋联合体积修复。不是原111块上多做精度扫描。
- Schur变换是可逆局部微分变换，导数接触Ω保留；q不要求精确为零，只用其光滑性。字段配对变化另计，不用发散迹循环。
- 原D1/D0保全部混合，唯一compatible connection允许非正定配对。矩阵系数的双点伴随对称有准确文献映射。
- 正常双曲递推按D=+Box+E、sigma=半间隔、nu=1/(8 pi^2)重推；左右jet及scalar Moretti 1/3回检。
- 体积接触是来源一形式，不把delta integral b误写为仅体积项；全相互作用插入可积性/QME仍未证明。
- 仅签收固定背景一圈一次来源的局部Wick规范化，包含原ghost和735费米。未声称规范固定独立、全阶背景独立、有限强度闭合或全量子引力。
- 校准W=0不是物理态；原实际物理态在解析构造中保留。非交换矩阵势的输运检验是真实有限Hadamard jet，但不是原753背景的数值模拟。
- 开发阶段delta K曾误用实矩阵，与原Fourier K=i*K_real不一致，使实接触校准为零；改为同一Fourier约定后非零，未改物理耦合或掩盖理论失败。结果首次冻结前完成修正。
- 旧空间382—386/425/522—523与604/649/699不重证或扩大。主代理逐式和代码复核；无独立代理审稿，无图像，无目标改动。
''')
write('round771_drafts/literature_scope_audit.json', json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/1904.03708', sections='I and Theorem 1',
             use='Sesqui-symmetry for the original formally self-adjoint matrix wave operators with nondegenerate, possibly indefinite pairing.',
             boundary='No positive Euclidean quantum-gravity measure or interacting state follows.'),
        dict(url='https://arxiv.org/html/1407.1994', sections='4, Propositions 4.3 and 4.8',
             use='Normally hyperbolic matrix transport remainder method, rederived in the original plus-Box and half-distance conventions.',
             boundary='No chiral-fermion anomaly result is transplanted; fermion scope stays at 735.'),
        dict(url='https://arxiv.org/html/hep-th/0306138', sections='2.1, 2.3 and equation 4.28',
             use='Unique connection/endomorphism and local heat coefficient with all original mixing.',
             boundary='Local coefficient identity only, no substitution of a Euclidean vacuum for the actual state.'),
        dict(url='https://arxiv.org/html/gr-qc/0109048', sections='Theorem 2.1 and Appendix B',
             use='Scalar 1/3 normalization cross-check after the matrix derivation.',
             boundary='Scalar theorem alone is not the joint construction.')],
    local_one_insertion_Ward_repair_constructed=True,
    complete_original_mixed_operators_retained=True,
    original_state_unchanged=True,
    interacting_QME_proven=False,
    gauge_fixing_independence_proven=False,
    all_insertion_common_normalization_proven=False,
    finite_strength_self_consistency_proven=False,
), ensure_ascii=False, indent=2)+'\n')

main = ['research_note_771.md', 'joint_local_ward_repair.py',
        'joint_local_ward_repair_results.json', 'unified_physics_condition_ledger_771.md']
review = '771主代理终审：十九式、三组复算及原对象映射；不是独立审稿。\n'
review += '核：Schur与Omega；密度/纤维配对接触；Hadamard矩阵左右jet；ghost权−1；局部体积一形式；真实态差保持；一圈来源Ward不升级成QME或有限强度闭合。\n'
review += '\n'.join(p+' '+hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in main)+'\n'
write('round771_drafts/final_review.txt', review)

mapping = {'770':'771', '769':'770', '3496':'3499', '3494':'3496',
           '1622':'1625', '3757':'3773', '3742':'3757',
           'joint_local_brst_subtraction':'joint_local_ward_repair'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda s: pattern.sub(lambda m: mapping[m[0]], s)
verify = convert((HERE/'verify_round770.py').read_text('utf8'))
start, stop = verify.index('    extra = ('), verify.index('    new = ')
verify = verify[:start]+'''    extra = ('unified_physics_condition_ledger_771.md',
             'round771_drafts/research_note_771_draft.md',
             'round771_drafts/final_review.txt',
             'round771_drafts/literature_scope_audit.json',
             'round771_drafts/scope_and_dedup_review.md',
             'round772_drafts/STATUS.md',
             'round771_drafts/research_note_771_working.md',
             'round771_drafts/local_contact_entry.py',
             'round771_drafts/local_contact_entry_results.json',
             'round771_drafts/entry_scope_review.json',
             'round771_drafts/entry_condition_ledger.md',
             'round771_drafts/local_contact_entry_checks.json',
             'round771_drafts/entry_postpublication_checks.json')
    entry = core.read(HERE/'round771_drafts/local_contact_entry_checks.json')
    for p, digest in entry['working_artifact_hashes'].items():
        assert core.digest(HERE/p) == digest, p
'''+verify[stop:]
verify = verify.replace('== (2, 0, 0)', '== (3, 0, 0)').replace('fresh_tests=dict(run=2,', 'fresh_tests=dict(run=3,')
verify = verify.replace("checks['display_formulas'] == 16", "checks['display_formulas'] == 19")
start, stop = verify.index('                exact_free_BRST_subtraction_proven='), verify.index('                primary_code_and_note_review_completed=')
verify = verify[:start]+'''                original_same_free_state_preserved=True,
                complete_Schur_and_pairing_contacts_retained=True,
                original_normal_operator_Ward_jets_constructed=True,
                one_loop_one_insertion_local_Ward_repair_proven=True,
                all_insertion_common_normalization_proven=False,
                gauge_fixing_independence_proven=False,
                finite_strength_self_consistency_proven=False,
                interacting_QME_proven=False,
                original_Q_to_E_process_equivalence_proven=False,
'''+verify[stop:]
write('verify_round771.py', verify)
publish = convert((HERE/'publish_round770.py').read_text('utf8'))
start, stop = publish.index('summary = '), publish.index('planned = ')
summary = '**第771轮完成（完整一圈来源的局部Ward修复）：** [同态来源与联合守恒]({p}research_note_771.md)保原耦合算符、辅助/配对接触与实际态，以共同局部体积项关闭当前一次来源的Ward。全相互作用QME、规范固定独立和有限强度闭合仍未证。三组、十九式通过，最新771／3499，1625份编号科学文件、3773份保护证据。[核验]({p}research_round_771_checks.json)、[条件账]({p}unified_physics_condition_ledger_771.md)。'
order = '**当前执行顺序（771后，优先于下方历史安排）：** 接[772绝对来源与共同首阶响应]({p}round772_drafts/STATUS.md)，核原完整约束、同一准备与实际过程；不重做局部Ward或矩阵精度。[范围审计]({p}round771_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
publish = publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish = publish.replace('第771轮完成（同态局部减除与来源候选）', '第771轮完成（完整一圈来源的局部Ward修复）')
publish = publish.replace('## 416.', '## 417.').replace('## 321.', '## 322.')
publish = publish.replace('同态BRST减除与局部Ward来源', '同态完整来源的局部Ward修复')
publish = publish.replace('旧空间合同保持，同态减除与局部规范化', '旧空间合同保持，完整来源与共同响应')
publish = publish.replace('BRST减除及局部Ward余项', '同态来源与联合守恒')
publish = publish.replace('next_round=771', 'next_round=772')
write('publish_round771.py', publish)
post = convert((HERE/'postcheck_round770.py').read_text('utf8'))
post = post.replace('range(584, 771)', 'range(584, 772)')
post = post.replace('第771轮完成（同态局部减除与来源候选）', '第771轮完成（完整一圈来源的局部Ward修复）')
write('postcheck_round771.py', post)
print('Prepared 771 evidence and publication helpers.')
