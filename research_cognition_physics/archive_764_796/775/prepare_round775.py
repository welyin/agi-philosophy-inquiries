"""Freeze the direct source/BV result and prepare guarded publication."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent


def write(name, value):
    p = HERE/name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as f:
        f.write(value)


note = (HERE/'research_note_775.md').read_text('utf8')
write('round775_drafts/research_note_775_draft.md', note)
ledger = (HERE/'unified_physics_condition_ledger_774.md').read_text('utf8').splitlines()
ledger[0] = '# 联合条件总账：775同态Wick来源与BV反常首项'
ledger[2] = '2026-10-04。接[774全账](unified_physics_condition_ledger_774.md)，回填[775报告](research_note_775.md)。[结果](joint_wick_source_anomaly_results.json)、[核验](research_round_775_checks.json)。原自由体系直接给实际来源与反常首项匹配；不必先有完整in-in有效作用。高次N1/N2仍开放。'
updates = {
    'C01': '775以同一适配自由态直接核Wick/BV首项；高次相互作用物理正态、荷及仪器仍待构造',
    'C11': '原完整三阶作用的一次Wick重排给770来源，其规范投影正是局部反常首项；771指定首项可用，不再以完整N3作此前置',
    'C19': '774关系插入提升继续有条件；775只关闭实际来源首项接口，未把理论插入实现为记录',
    'C22': '原111体系Schur辅助场的零态协方差与非零Feynman接触已区分；N1/N2分布延拓、字段独立与QAP共同映射仍须核验',
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
write('unified_physics_condition_ledger_775.md', '\n'.join(ledger)+'''

## 775当前结论（优先于保留的774执行安排）

原完整自由Wick/BV身份直接给三阶作用的来源字典：A的一圈线性项为J_H(Rc)。因此771的指定局部接触是773/774保首项递推所需的同一个C₁。它在772适配扩展上使用；与适配前比较须共同运输均值和初值。735费米局部差仍是存在构造，未列全数值系数。

高次逐图规范化显式保持ℏ阶与外腿数；增加的ℏ局部作用另行记账。完整N3初态边界问题对局部首项不再是前置，却未因此得到全局证明。原相互作用态、真实仪器与Q连续仍开放。

下一项接[776](round776_drafts/STATUS.md)：保原辅助代数接触，核N1/N2在原混合系统与773正则背景jet类中的共同实现。若成功，直接复用773/774，不重新做同调、首项或旧空间证明。理论维数、群及经典作用的输入属性保持。
''')
write('round776_drafts/STATUS.md', '''# 第776轮入口：原混合体系的高次因果规范化

接[775直接连接](../research_note_775.md)、[775条件总账](../unified_physics_condition_ledger_775.md)。

1. 775已从同一自由态、原三阶作用、Wick减除证明来源与反常首项相容。不要再次以完整in-in有效作用或全局迹类密度矩阵作为这个局部接口的必要前置。
2. 优先核N1/N2：原111混合场、完整费米、反场和frame，是否可在同一因果/实/字段独立时间序体系中满足基本方程与局部量子作用原理。Schur辅助协方差虽为零，其Feynman接触−iδ必须保留。
3. 773的逆参考Jacobian、逆Gram仅在正则片内是光滑背景有限jet系数。每固定阶局部，不坚持旧UV幂计数，也不擅称全局或背景独立。核成熟延拓定理能否容纳这个类以及实际导数插入。
4. 映射成立后直接复用773/774的实局部一圈递推，保771首项及共同关系插入。无需逐项列出所有反常系数才允许作存在证明；也不能把尚未计算称为系数为零。
5. 完整相互作用正态、实际准备/记录及Q连续映射仍另行核验。不重复775有限Wick诊断、自由主符号或旧空间382—386、425、522—523。
6. 604、649/699的限定失败及总目标保持。无目标改写、新应用任务、定时或图像检查。
''')
write('round775_drafts/scope_and_dedup_review.md', '''# 775范围与去重审查

新增：原完整三阶BV作用经同态Wick重排给770来源，直接识别一圈线性反常；保771首项不再需要先建完整N3。原辅助Feynman代数接触也已明确。

复用：735费米局部修复存在性，768完整自由身份，770一外腿三价菜单，771局部来源Ward，772适配态和均值运输，773局部同调，774条件递推。没有将这些旧结果重算为新定理。

显式范围：在壳正则片、紧支体内测试、原物理输入、逐图保外腿/圈次数、固定背景的一次局部反项。两个数值/精确代数组只校准身份和错误对照，不求原连续圈系数。

未证明：高次N1/N2共同实现、完整N3/初态边界、全部原相互作用QME、正态/算符域/实际仪器、严格UV与全局背景独立、Q连续。原空间及604、649/699结论不改写。

审查方式：主代理读原代码、完整报告及原始文献框架，并复算已存结果；无独立代理或图像检查。
''')
audit = dict(primary_agent_review=True,
    direct_original_cubic_source_anomaly_match=True,
    prescribed_source_firstjet_compatible=True,
    complete_in_in_functional_not_needed_for_this_local_firstjet=True,
    complete_in_in_or_initial_boundary_identity_proven=False,
    original_full_loop_anomaly_coefficients_computed=False,
    higher_N1_N2_joint_realization=False,
    sources=[
        dict(url='https://arxiv.org/pdf/1803.10235',
             locations=['Theorem 3', 'Equations (137)-(140) and Theorem 5', 'Theorems 7-9'],
             use='Quadratic BV derivation, antifield noncontraction, anomalous first-insertion identity, consistency and equation normalization.',
             limitation='Higher original mixed-system normalization remains to be mapped; actual source identification is the present argument.'),
        dict(url='https://arxiv.org/html/0705.3160',
             locations=['Section 4, Proposition 4', 'Sections 5.3 and 5.4.1, Proposition 13', 'Appendix A'],
             use='Proper-vertex and quantum action principle framework; full HTML inspected after the previous abstract-only audit.',
             limitation='Scalar framework does not directly establish the complete original gauge/gravity system or a prepared in-in state.'),
        dict(url='https://arxiv.org/html/2604.26941v1',
             locations=['Sections 3.4-3.5'],
             use='Full prepared-state boundary problem kept separate from local first-jet identification.',
             limitation='Not needed as a premise of the new first-jet proof; no global interacting-state identity is claimed.'),
    ])
write('round775_drafts/literature_scope_audit.json', json.dumps(audit, ensure_ascii=False, indent=2)+'\n')
main = ['research_note_775.md', 'joint_wick_source_anomaly.py',
        'joint_wick_source_anomaly_results.json', 'unified_physics_condition_ledger_775.md']
review = '775主代理最终范围审查：直接首项连接成立；高次共同规范化、完整相互作用态仍未签收。无独立代理或图像检查。\n'
review += '\n'.join(p+' '+hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in main)+'\n'
write('round775_drafts/final_review.txt', review)

mapping = {'774': '775', '773': '774', '775': '776',
           '3508': '3510', '3505': '3508', '1634': '1637',
           '3806': '3815', '3815': '3829',
           'joint_source_action_lift': 'joint_wick_source_anomaly'}
regex = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda txt: regex.sub(lambda match: mapping[match[0]], txt)
verify = convert((HERE/'verify_round774.py').read_text('utf8'))
verify = verify.replace("'round776_drafts/STATUS.md')", """'round776_drafts/STATUS.md',
             'round775_drafts/research_note_775_working.md',
             'round775_drafts/entry_condition_ledger.md',
             'round775_drafts/entry_scope_review.json',
             'round775_drafts/causal_source_entry_checks.json',
             'round775_drafts/entry_postpublication_checks.json')""")
verify = verify.replace('    result = model.run()', """    entry = core.read(HERE/'round775_drafts/causal_source_entry_checks.json')
    for name, digest in entry['working_artifact_hashes'].items():
        assert core.digest(HERE/name) == digest, name
    result = model.run()""")
verify = verify.replace("== (3, 0, 0)", "== (2, 0, 0)").replace('fresh_tests=dict(run=3', 'fresh_tests=dict(run=2')
verify = verify.replace("checks['display_formulas'] == 15", "checks['display_formulas'] == 14")
start, stop = verify.index('                original_physical_free_state_preserved='), verify.index('                primary_code_and_note_review_completed=')
verify = verify[:start]+'''                original_physical_free_state_preserved=True,
                direct_original_cubic_source_anomaly_match=True,
                prescribed_source_firstjet_compatible=True,
                auxiliary_Feynman_contact_retained=True,
                complete_in_in_or_initial_boundary_identity_proven=False,
                original_full_loop_anomaly_coefficients_computed=False,
                higher_N1_N2_joint_realization=False,
                original_interacting_QME_unconditionally_proven=False,
                strict_global_polynomial_counterterm_result=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
'''+verify[stop:]
write('verify_round775.py', verify)
publish = convert((HERE/'publish_round774.py').read_text('utf8'))
summary = '**第775轮完成（同态来源与反常首项）：** [原Wick/BV直接匹配]({p}research_note_775.md)同一三阶作用直接给原来源的反常投影，771指定首项可接774；无需先建完整in-in有效作用。高次N1/N2和相互作用态仍待验收。两组、十四式通过，最新775／3510，1637份编号科学文件、3829份保护证据。[核验]({p}research_round_775_checks.json)、[条件账]({p}unified_physics_condition_ledger_775.md)。'
order = '**当前执行顺序（775后，优先于下方历史安排）：** 接[776高次因果规范化]({p}round776_drafts/STATUS.md)，保原辅助Feynman接触，核混合体系的分布延拓、字段独立及局部QAP；不重复首项匹配或773同调。[范围审计]({p}round775_drafts/scope_and_dedup_review.md)。原物理准备、旧空间与目标保持。'
start, stop = publish.index('summary = '), publish.index('planned = ')
publish = publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish = publish.replace('第775轮完成（来源与量子作用的条件性连接）', '第775轮完成（同态来源与反常首项）')
publish = publish.replace('## 420.', '## 421.').replace('## 325.', '## 326.')
publish = publish.replace('来源与共同量子作用的条件性连接', '同一Wick来源与BV首项直接连接')
publish = publish.replace('旧空间合同保持，条件性一圈共同插入', '旧空间合同保持，来源首项已直接匹配')
publish = publish.replace('保来源首项与共同插入', '同态来源与反常首项')
write('publish_round775.py', publish)
post = convert((HERE/'postcheck_round774.py').read_text('utf8'))
post = post.replace('第775轮完成（来源与量子作用的条件性连接）', '第775轮完成（同态来源与反常首项）')
write('postcheck_round775.py', post)
print('Prepared direct round 775 result; old evidence unchanged.')
