"""Prepare 764 artifacts and derive checked publication tools from 763."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(name, text):
    p = HERE/name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x', encoding='utf8', newline='\n') as f:
        f.write(text)


lines = (HERE/'unified_physics_condition_ledger_763.md').read_text('utf8').splitlines()
lines[0] = '# 联合条件总账：764共同线性物理相空间与内禀正量子态'
lines[2] = '2026-10-04。接[763全账](unified_physics_condition_ledger_763.md)，回填[764报告](research_note_764.md)。[结果](joint_linear_physical_phase_results.json)、[核验](research_round_764_checks.json)。新增E微扰量子化分支，不等同原Q已完成连续匹配。'
updates = {
    'C01': '764原完整经典线性约束可辛约化并构造正则正Weyl态；领先费米CAR可接同一背景',
    'C03': '764仅有物理代数中的抽象仪器可能性，原L_r(s)的关系实现和注能字典仍须证明',
    'C04': '764原短时经典线性发展诱导物理CCR同构；不保证单一Fock表示内全时间酉实现',
    'C09': '764给R到各向同性补空间S及规范截面P的显式连接；截面非物理自主纠错操作',
    'C11': '764同一物理辛空间包含几何内禀方向；全非线性量子约束与排序/反常未完成',
    'C19': '764辅助正协方差仍为新增准备；regular不等于既定Fock的normal，也不保证Hadamard',
    'C22': '764同阶来源含玻色约束Hessian二次项及费米来源；735单费米圈不能覆盖全部部门'}
seen = []
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith('|'+key+' '):
            parts = line.split('|')
            parts[2] += '；'+value
            lines[i] = '|'.join(parts)
            seen.append(key)
assert set(seen) == set(updates)
write('unified_physics_condition_ledger_764.md', '\n'.join(lines)+'''

## 764从完整经典约束到共同领先量子对象

763的诱导响应不含独立齐次几何资料。764明确采用E微扰量子化分支：将原g、pi、规范A/E及五标量phi/p纳入同一canonical切空间，复用731/754的光滑右逆，修成各向同性补空间并构造辛规范截面。完整线性约束核除规范方向的配对非退化；原短时发展保持此配对。

在该物理空间上选辅助正协方差，存在正的正则Weyl态；原零经典费米背景使领先Hessian的费米交叉块消失，可接同一背景CAR正态。该分支给领先自由共同对象，不把原固定图Gauss物质态、590控制几何或763来源派生量混作同一量子引力态。

新增输入是微扰量子化及准备选择。Hadamard短距离、二次来源重整化、完整量子约束、相互作用和原实际记录映射尚未完成。局部31/29自由度对照来自既定维数、群和作用，不是认知推导或全局模式维数。数值只校准导数对角核心和有限辛恒等式，没有模拟全部变系数约束。

## 下一项：物理短距离态、同阶全来源与任务映射

接765，核当前耦合物理CCR准备如何同时有足够的短距离性质、Ward身份和原来源/记录字典。原平方根尺度下，二次玻色来源与费米来源同阶，不能只回填735一费米圈。优先成熟的线性规范场Hadamard/微扰方法及其真正前提，不继续投影矩阵、自由度或读口常数扫描。旧空间、604、649/699及统一目标保持。
''')
write('round764_drafts/research_note_764_draft.md', (HERE/'research_note_764.md').read_text('utf8'))
write('round765_drafts/STATUS.md', '''# 第765轮入口：共同物理态的短距离条件与同阶完整来源

接[764](../research_note_764.md)、[条件账](../unified_physics_condition_ledger_764.md)与[H1—H3](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 764已在原753经典背景上，从完整约束右逆构造辛截面、物理CCR及正则正态；领先费米CAR可接同一背景。该E微扰量子化是新增分支，不是原Q连续等价。
2. 任意辅助Gaussian协方差虽正，并不保证Hadamard或二次来源可重整化。优先核完整耦合几何/规范/五标量的短距离态、规范相容与共同Ward；735仅一费米圈，不能替代。
3. 原缩放的v在sqrt(epsilon)阶，D²C[v,v]/2与费米q_F在epsilon同阶。原R处理剩余初始系数，但不自行提供非线性CCR、量子约束无反常或有限epsilon误差。
4. 复用573发展、731/754右逆、735范围、752已有BV文献审计及764物理空间；不要重做一般投影或全态赋态反例。
5. 仍须给原Q的来源乘积、未知输入及实际记录到E物理代数的同一映射。Weyl代数可定义抽象CP仪器不等于原L_r(s)已经接通。
6. 数值主阶对角核心不是全部变系数主符号，31/29是局部自由度计数，不作空间推导。正则regular亦不等于既定Fock的normal。
7. 目标未改；旧空间382—386、425、522—523及604、649/699保持。不新建应用任务、调度或图像。
''')
write('round764_drafts/scope_and_dedup_review.md', '''# 764范围与去重审查

- 前一目标轮已完成763及764对象回查，为有效进展。本轮核导航、最新报告/结果、运行进程及未冻结工作内容后接续。
- 568规范辛约化和650边界辛通量仅作历史接口；不把它们直接视为当前全引力约化。574/758的gamma外给，590是另加控制分支。
- 完整经典C的在壳第一类恒等式与731/754光滑右逆共同给本轮截面。伴随在同一canonical密度配对上定义；原R由光滑微分/椭圆/有限秩构成，其转置保持光滑。无边界T³是实质前提。
- R不是态通道，P不是自主修复。I-RL保规范方向，本轮另将右逆像修成各向同性并除去规范方向；弱非退化由公式证明，不靠数值秩。
- 正则Gaussian Weyl态是新增准备。没有签收Hadamard、正能量基态、既定Fock中的normal或全时间酉实现；也没有强加kinematical正常零约束本征态。
- 领先费米解耦依赖经典费米背景零及大作用量分阶；残余有限稳定子识别保留。相互作用、原Q匹配及实际记录不能由自由乘积态领取。
- 数值仅通用导数对角核心和有限校准；显式注明省略背景相关导数混合项。初稿代码/结果保留，避免把它称作完整变系数主符号。
- 同阶约束包含玻色二次来源。有限模对称排序只给形式系数恒等式，不给连续Wick构造、量子约束闭合、非线性CCR或全域正度规。
- 旧空间与604、649/699边界、完整统一目标均保持。三组新检查不复计工作对象回查；不再延伸矩阵数值扫描。
''')
write('round764_drafts/literature_scope_audit.json', json.dumps({
    'sources': [
        {'url': 'https://arxiv.org/html/1205.3484',
         'checked': 'Definition 3.4 distinguishes the full gauge complex; Definition 4.3/Theorem 4.5 and Remarks 4.6-4.7 supply Weyl algebra quantization and caution on centres.',
         'use': 'Mature algebraic interface only. Coupled canonical reduction and positive covariance are proved explicitly in note764.',
         'not_imported': 'Separate vacuum GR / Yang-Mills examples are not treated as a ready-made all-matter theorem; no Hadamard or full quantum Einstein conclusion.'}],
    'historical_reference': '752 already audited BFR perturbative BV scope; 735 controls only the declared one-fermion-loop sector.',
    'new_project_connection': 'The existing full canonical constraint right inverse gives an isotropic complement and physical symplectic section, which admits a positive regular CCR preparation.',
    'scope': 'Linear/full-Cauchy-surface algebra and leading free branch; no original-instrument matching or complete interacting quantum gravity.'
}, ensure_ascii=False, indent=2)+'\n')
files = ('research_note_764.md', 'joint_linear_physical_phase.py',
         'joint_linear_physical_phase_results.json', 'unified_physics_condition_ledger_764.md')
write('round764_drafts/final_review.txt',
      'Primary agent review; no independent agent. Checked J signs and adjoints: G=JL*, L R=I, B=R*JR, S=R+GB/2, and P=I-SL+GS*J. Verified weak nondegeneracy without relying on numerical rank. R and its transpose are smooth operators in the inherited compact canonical setting; no blanket Sobolev bound is asserted. Gaussian positivity is for the reduced CCR, not an unreduced normal constraint eigenstate. The preparation is regular, not certified Hadamard/normal in an existing Fock representation. Numerical matrices omit background-dependent derivative mixing and are only diagonal-core calibrations. Leading CAR factorization uses zero classical fermions; finite stabilizers and higher-order coupling remain. Full sources need bosonic quadratic terms at the same order and require new renormalization/task matching.\n'
      +'\n'.join(name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files)+'\n')

# Explicit one-pass substitutions avoid cascaded round/counter replacements.
mapping = {'763': '764', '762': '763', '3653': '3666', '3666': '3679',
           '3473': '3476', '3476': '3479', '1601': '1604',
           'joint_operator_constraint_lift': 'joint_linear_physical_phase'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
def convert(s):
    return pattern.sub(lambda m: mapping[m[0]], s)

verify = convert((HERE/'verify_round763.py').read_text('utf8'))
a, b = verify.index('    extra = ('), verify.index('    new = ')
extra = """    extra = ('unified_physics_condition_ledger_764.md',
             'round764_drafts/research_note_764_draft.md',
             'round764_drafts/final_review.txt',
             'round764_drafts/literature_scope_audit.json',
             'round764_drafts/scope_and_dedup_review.md',
             'round765_drafts/STATUS.md',
             'round764_drafts/research_note_764_working.md',
             'round764_drafts/object_inventory_checks.json',
             'round764_drafts/joint_linear_physical_phase_before_scope_review.py',
             'round764_drafts/results_before_scope_review.json')
"""
verify = verify[:a]+extra+verify[b:]
verify = verify.replace("checks['display_formulas'] == 13", "checks['display_formulas'] == 15")
a = verify.index('                source_and_response_same_algebra=True,')
b = verify.index('                primary_code_and_note_review_completed=True,')
verify = verify[:a]+'''                full_classical_linear_reduction_proven=True,
                positive_regular_CCR_state_constructed=True,
                auxiliary_covariance_is_new_input=True,
                numerical_checks_only_derivative_core_and_finite_symplectic=True,
                complete_background_dependent_constraint_matrix_simulated=False,
                Hadamard_bosonic_state_certified=False,
                nonlinear_quantum_constraints_completed=False,
                original_Q_to_E_process_equivalence_proven=False,
'''+verify[b:]
write('verify_round764.py', verify)

publish = convert((HERE/'publish_round763.py').read_text('utf8'))
summary = '**第764轮完成（领先物理量子分支）：** [共同线性物理相空间与内禀量子态]({p}research_note_764.md)由原完整约束右逆构造辛截面、物理CCR及正则正态，接同背景领先费米部门。微扰量子化与准备是明示输入；Hadamard全来源、实际记录及原Q连续匹配仍开放。三组、十五式通过，最新764／3479，1604份编号科学文件、3679份保护证据。[核验]({p}research_round_764_checks.json)、[条件账]({p}unified_physics_condition_ledger_764.md)。'
order = '**当前执行顺序（764后，优先于下方历史安排）：** 接[765物理短距离态与同阶完整来源]({p}round765_drafts/STATUS.md)，核耦合物理态、玻色二次来源及费米来源的共同Ward和任务字典；不继续投影或自由度扫描。[范围审计]({p}round764_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及总目标保持。'
publish = re.sub(r'^summary = .*$', lambda m: 'summary = '+repr(summary), publish, flags=re.M)
publish = re.sub(r'^order = .*$', lambda m: 'order = '+repr(order), publish, flags=re.M)
publish = publish.replace('**第764轮完成（条件性初始连接）：**', '**第764轮完成（领先物理量子分支）：**')
publish = publish.replace('409.', '410.').replace('314.', '315.')
publish = publish.replace('同代数初始响应与几何因子边界', '共同线性物理相空间与正态')
publish = publish.replace('同一来源与受约束响应', '共同量子相空间与同阶来源')
publish = publish.replace('同代数受约束响应与几何因子边界', '共同线性物理相空间与内禀量子态')
publish = publish.replace('next_round=764', 'next_round=765')
write('publish_round764.py', publish)

post = convert((HERE/'postcheck_round763.py').read_text('utf8'))
post = post.replace('range(584, 764)', 'range(584, 765)')
post = post.replace('**第764轮完成（条件性初始连接）：**', '**第764轮完成（领先物理量子分支）：**')
post = post.replace(" == publication['navigation_links'] == 12249", " == publication['navigation_links']")
write('postcheck_round764.py', post)
print('Prepared 764 metadata and publication tools; 765 entry saved.')
