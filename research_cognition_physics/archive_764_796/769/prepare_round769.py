"""Prepare 769 evidence and publishers using exclusive creates."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(name, content):
    p = HERE/name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)


lines = (HERE/'unified_physics_condition_ledger_768.md').read_text('utf8').splitlines()
lines[0] = '# 联合条件总账：769同态换帧、均值接触与共同来源'
lines[2] = '2026-10-04。接[768全账](unified_physics_condition_ledger_768.md)，回填[769报告](research_note_769.md)。[结果](joint_quantum_frame_transport_results.json)、[核验](research_round_769_checks.json)。同一自由态和指定相对来源的换帧已接通；绝对全部门来源仍开放。'
updates = {
    'C01': '769原F>0换帧及实际拉回规范固定输送767同一物理Hadamard态，不另选准备',
    'C09': '769同一报告的二阶字典包含均值接触与协方差两项；尚非实际自主仪器',
    'C11': '769指定平滑物理态差的来源随Hessian接触变换，与均值/初值共同保首阶约束',
    'C19': '769共同自由物理态的换帧成立；独立连续量子测度和绝对Wick处方等价未证',
    'C22': '769相对来源和响应跨原两帧共同输送；BV反常类别只在匹配函数空间中同构，绝对来源仍开放',
}
seen = set()
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith('|'+key+' '):
            parts = line.split('|'); parts[2] += '；'+value
            lines[i] = '|'.join(parts); seen.add(key)
assert seen == set(updates)
write('unified_physics_condition_ledger_769.md', '\n'.join(lines)+'''

## 769：同一量子模型的原换帧字典

在原753在壳F>0条带中，原Jordan/Einstein变量、规范固定、辅助系统和物理态可以一起输送。自由完整Green/BRST/Hadamard及物理商正性由可逆局部变换保持。并未比较两个独立选择的量子测度。

768指定平滑物理协方差差决定非线性平均场接触k；完整相对来源的正确变换是T转置乘(Delta J_E+H_E k)，平均场则变为T Delta w_J+k。原形式首阶方程及约束由此一起对应，既定报告的二阶字典也保同一读数。绝对奇异核不可直接代入此平滑差公式。

同一局部BV函数类别的可逆换帧不消除非平凡反常类。经典Jordan作用多项式不自动使全部EFT局部项满足有限物质多项式合同，曲率平方拉回即有F逆幂。独立测度、反常原始元类别、全部门绝对局部来源及实际in-in连接仍未签收。未新增认知公理、物种或联络，未改变目标。

## 下一项：770固定背景全部门绝对来源

以769工作入口的最小首阶验收为准，回到原完整三阶作用和局部减除，核共同规范化与反常原始元类别，回用628/735。不要先要求全部高阶背景独立，也不要重复场重定义、物理态或本轮均值接触。原Q、实际记录、内部资源和空间旧接口保持。
''')
write('round769_drafts/research_note_769_draft.md', (HERE/'research_note_769.md').read_text('utf8'))
write('round770_drafts/STATUS.md', '''# 第770轮入口：固定背景首阶绝对来源的共同规范化

接[769](../research_note_769.md)、[全账](../unified_physics_condition_ledger_769.md)及[769绝对来源工作入口](../round769_drafts/research_note_769_working.md)。

1. 原753在壳背景、767/768同一自由物理态及辅助扩展保持。769已输送原两帧的相同自由对象和指定相对来源，不再重复T、物理正性或均值接触计算。
2. 下一验收是同一固定背景的全部门绝对局部来源及其在壳联合Ward，先做形式首阶，再接754。不是先解决全阶背景独立性。
3. 成熟BV局部异常定理给局部缺陷，不自动给其可消去性。回用628费米反常抵消和735实际费米处方；实核原U1、重力、H5的同一函数类别及反常原始元。
4. 可选择有限EFT/外线阶的形式合同，但必须说明它与原光滑背景首阶来源的关系；不能将形式任意阶可解直接当作函数收敛或一个实际局部处方。
5. 换帧应同时拉回规范固定、减除、复合字段与来源。绝对均值接触有处方依赖，不能把769平滑差公式用于奇异二点核。也不能把状态依赖正常序当作共同局部减除而令真空来源为零。
6. 反常抵消、相互作用量子主方程、背景分裂比较和有限epsilon闭合分开验收。保持原Q连续字典、实际记录、内部资源及旧空间382—386/425/522—523；604和649/699范围不变。无新目标、任务、调度或图像。
''')
write('round769_drafts/scope_and_dedup_review.md', '''# 769范围与去重复核

- 直接复用583联合Hessian、600共同源/观测重定义、619费米权重、664原联络范围，不重报一般换帧定理。
- 新接口：当前767/768同一完整物理态及其指定平滑变化实际接入原Jordan/Einstein变量，来源差、均值接触和约束响应共用同一协方差。
- 密度Hessian和纤维配对分开；可变系数T的导数保留；实际规范固定也拉回，未比较独立选定的两个规范。
- P双解只施加到768平滑物理差，未施加到整份未约束二点对象。
- 原完整作用三阶链式法则在在壳背景、完整双解条件下化为来源式(10)；数值离壳势块则保留式(16)全部接触，不伪装在壳全波解。
- 平滑差k无UV歧义，绝对k需要共同Wick处方。来源和平均场不是各自不变量，组合及同一报告才保持。
- 物理态的自由代数同构不等于非线性量子主方程、两个独立测度等价或实际自主读取。
- 局部BV同构限定为匹配函数类别，经典Jordan多项式不自动处理全部量子EFT菜单。未宣布已计算异常或已将全异常消去。
- 三组代码是有限代数校准。111字段用原主部；两点分布和坐标平方不当真实Hadamard态/物理记录；密度jet不当圈计算。
- 769入口的旧两组工作校准保持冻结，不重新计入本轮三个测试组。
- 主代理逐式复核，未使用独立代理或图像检查。当前目标保持。
''')
write('round769_drafts/literature_scope_audit.json', json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/1811.09413', sections='2 and 6',
             use='Joint field/source/Jacobian changes and limitations of ignoring them.',
             boundary='No unqualified unit Jacobian for the original finite F-dependent transformation; no independent continuum measure comparison.'),
        dict(url='https://arxiv.org/html/hep-th/9505173', sections='2, local functions and action assumptions',
             use='Matter/jet polynomial category with smooth undifferentiated vielbein dependence.',
             boundary='Pulled curvature-squared counterterms have inverse F. A polynomial classical Jordan action alone does not establish the full admissible anomaly/counterterm category.')],
    reused_project_rounds=[583,600,619,664,765,766,767,768],
    same_free_physical_state_transport=True,
    smooth_relative_source_and_mean_contact=True,
    absolute_source_exists_proven=False,
    measure_equivalence_proven=False,
    all_sector_quantum_Ward_proven=False,
), ensure_ascii=False, indent=2)+'\n')

main = ['research_note_769.md', 'joint_quantum_frame_transport.py',
        'joint_quantum_frame_transport_results.json', 'unified_physics_condition_ledger_769.md']
review = '769主代理终审：十六式、三组复算及源代码已核。非独立审稿。\n'
review += '重点：完整拉回规范固定；物理商正性；仅平滑P双解使用消项；均值与来源接触同号；有限量子强度和绝对Ward未签收。\n'
review += '\n'.join(p+' '+hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in main)+'\n'
write('round769_drafts/final_review.txt', review)

mapping = {'768':'769','767':'768','3491':'3494','3488':'3491','1616':'1619','3714':'3726','3726':'3742',
           'joint_brst_relative_source':'joint_quantum_frame_transport'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda x: pattern.sub(lambda m: mapping[m[0]], x)
verify = convert((HERE/'verify_round768.py').read_text('utf8'))
start, stop = verify.index('    extra = ('), verify.index('    new = ')
verify = verify[:start]+'''    extra = ('unified_physics_condition_ledger_769.md',
             'round769_drafts/research_note_769_draft.md',
             'round769_drafts/final_review.txt',
             'round769_drafts/literature_scope_audit.json',
             'round769_drafts/scope_and_dedup_review.md',
             'round770_drafts/STATUS.md',
             'round769_drafts/research_note_769_working.md',
             'round769_drafts/bv_source_entry_calibration.py',
             'round769_drafts/bv_source_entry_calibration_results.json',
             'round769_drafts/entry_condition_ledger.md',
             'round769_drafts/entry_scope_review.json',
             'round769_drafts/absolute_source_entry_checks.json',
             'round769_drafts/entry_postpublication_checks.json')
    entry = core.read(HERE/'round769_drafts/absolute_source_entry_checks.json')
    for p, digest in entry['working_artifact_hashes'].items():
        assert core.digest(HERE/p) == digest, p
'''+verify[stop:]
verify = verify.replace("checks['display_formulas'] == 18", "checks['display_formulas'] == 16")
start, stop = verify.index('                complete_linear_BRST_extension_constructed='), verify.index('                primary_code_and_note_review_completed=')
verify = verify[:start]+'''                original_free_state_frame_transport_proven=True,
                nonlinear_relative_mean_contact_proven=True,
                relative_source_and_response_transport_proven=True,
                matched_category_BV_isomorphism_only=True,
                absolute_all_sector_Ward_completed=False,
                independent_measure_equivalence_proven=False,
                finite_epsilon_nonlinear_semiclassical_closure_proven=False,
                original_Q_to_E_process_equivalence_proven=False,
'''+verify[stop:]
write('verify_round769.py', verify)

publish = convert((HERE/'publish_round768.py').read_text('utf8'))
start, stop = publish.index('summary = '), publish.index('planned = ')
summary = '**第769轮完成（同态换帧与来源接触）：** [原两帧的共同量子字典]({p}research_note_769.md)输送同一自由物理态；指定平滑变化的平均场接触恰补完整来源及首阶响应。经典作用多项式化未自动关闭反常类别与绝对来源。三组、十六式通过，最新769／3494，1619份编号科学文件、3742份保护证据。[核验]({p}research_round_769_checks.json)、[条件账]({p}unified_physics_condition_ledger_769.md)。'
order = '**当前执行顺序（769后，优先于下方历史安排）：** 接[770固定背景首阶绝对来源]({p}round770_drafts/STATUS.md)，核共同局部规范化及真实反常原始元类别；不重复换帧或均值接触，不要求先完成全阶背景独立。[范围审计]({p}round769_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
publish = publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish = publish.replace('第769轮完成（共同态与相对量子来源）', '第769轮完成（同态换帧与来源接触）')
publish = publish.replace('## 414.', '## 415.').replace('## 319.', '## 320.')
publish = publish.replace('同一物理态的BRST扩展与完整相对来源', '同态换帧、平均场接触与相对来源')
publish = publish.replace('旧空间合同保持，同态来源与绝对规范化', '旧空间合同保持，同态换帧与绝对来源')
publish = publish.replace('BRST扩展与完整玻色来源差', '原两帧的共同量子字典')
publish = publish.replace('next_round=769', 'next_round=770')
write('publish_round769.py', publish)
post = convert((HERE/'postcheck_round768.py').read_text('utf8'))
post = post.replace('range(584, 769)', 'range(584, 770)')
post = post.replace('第769轮完成（共同态与相对量子来源）', '第769轮完成（同态换帧与来源接触）')
write('postcheck_round769.py', post)
print('Prepared 769 note ledger, review, verification and publication helpers.')
