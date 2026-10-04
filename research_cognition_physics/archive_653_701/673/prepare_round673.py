"""Prepare673; preserve all prior evidence and the active goal."""
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


ledger=(HERE/'unified_physics_condition_ledger_672.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：672独立规范历史与共同电动能拼接',
    '# 联合条件总账：673共同边界、变秩字典与有限动态泛函')
ledger=ledger.replace('接[671全账](unified_physics_condition_ledger_671.md)，回填[672报告](research_note_672.md)。[结果](joint_gauge_history_kernel_results.json)、[核验](research_round_672_checks.json)。',
    '接[672全账](unified_physics_condition_ledger_672.md)，回填[673报告](research_note_673.md)。[结果](joint_gauss_boundary_functional_results.json)、[核验](research_round_673_checks.json)。')
ledger=ledger.replace('672独立历史与电热核的共同必要限制|',
    '672独立历史限制；673原H_b与全边界的有限未归一泛函|')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '673物理来源矩形字典保同一观测；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('672固定边界二步候选与643原Gauss有序过程仍须匹配',
    '673双边界候选已明确且绝对可积；其与643原Gauss有序过程仍须匹配')
ledger=ledger.replace('群选择、共享边界与动态Gauss过程仍缺|',
    '673全双边界Haar对象及全部字段输送已接；群选择、物理正性与原动态过程仍缺|')
ledger=ledger.replace('动态Gauss正性、指定归一与全局测度仍缺|',
    '673变秩字典与两片时间的几乎处处有限测度已接；一般反射身份、动态正性、指定归一与连续局域性仍缺|')
ledger=ledger.replace('动态过程、机制选择及现实参数|',
    '673全配置质量／物理来源绝对可积；动态过程、机制选择及现实参数|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 673将670正则物理来源身份扩展为矩形K、L；积分变量可有不同尺寸，原观测行数保持。有限未归一泛函不再要求固定等秩片；这不是更改原可观测结构。
- 两时间片、m0=1时单位时间链路保证任意空间配置的H可逆，故两个群边界中的Wilson零集为Haar零测度。结合投影有界和669配置热矩，可定义完整有限积分而不加硬统一谱隙；该结论不免除局域性、背景导数和连续极限的谱控制。

## 本轮合并与下一项

673连接C01、C04、C14、C16、C17、C19、C22：原16通道、全部质量及物理来源采用同一固定尺寸局部Berezin表示；矩形字典不需要逆K。两个独立Gauss边界均保留，将一处设为单位时，另一半φ、E、U和来源必须同步输送，原H5距离及权重身份已实际核验。

原完整H_b标量热核与当前overlap物质权重的双边界候选有明确公式。原配置热矩给全Q、S9与Haar积分的绝对支配，可交换积分并取有限质量来源系数。这是共同有限泛函，不是已经归一的正物理态；尚未证明一般反射身份、非负／非零总权重，也未证明等于643的原有序CAR影响。

原文2017已有全拓扑部门测度构想；本次仅补本项目非零质量与物理来源的实际字典，不新领取拓扑理论。较强场样本中两份完整费米矩阵接近奇异，其极小Pfaffian相位／相对误差无可靠含义，失败诊断和审查前稿保留。一般Hermitian／实性初稿表述已撤回，未用数值修补冒充证明。

接[674](round674_drafts/STATUS.md)：检验完整平均后的物理正性与原过程身份，任何分步／抽样近似均明确误差范围，不用低维子群或仅电因子代替完整问题。旧空间成果直接继承，认知设计后置，统一目标不改。
'''
write('unified_physics_condition_ledger_673.md',ledger)
mapping={'joint_gauge_history_kernel':'joint_gauss_boundary_functional',
         '671':'672','672':'673','673':'674','3249':'3251','3251':'3253',
         '1325':'1328','1328':'1331','2389':'2403','2403':'2419'}
verify=replace((HERE/'verify_round672.py').read_text('utf8'),mapping)
verify=verify.replace('Verify independent signed half histories and original-group electric sewing.',
    'Verify rectangular physical sources and full double-boundary finite functional.')
verify=verify.replace('import joint_static_gauge_reflection as prior_model','import joint_gauge_history_kernel as prior_model')
verify=verify.replace('dynamic_commutator_probe','boundary_holonomy_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_673.md','round673_drafts/research_note_673_draft.md',
           'round673_drafts/final_review.txt','round674_drafts/STATUS.md',
           'round673_drafts/boundary_holonomy_entry.md','round673_drafts/boundary_holonomy_probe.py',
           'round673_drafts/boundary_holonomy_probe_results.json',
           'round673_drafts/first_boundary_attempt.py','round673_drafts/first_boundary_results.json',
           'round673_drafts/source_boundary_attempt.py','round673_drafts/source_boundary_results.json',
           'round673_drafts/complex_reflection_diagnostic.json','round673_drafts/research_note_673_before_review.md')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==16","assert text['display_formulas']==18")
write('verify_round673.py',verify)
publish=replace((HERE/'publish_round672.py').read_text('utf8'),mapping|{'## 318.':'## 319.','## 223.':'## 224.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第673轮完成：** [共同规范边界、变秩手征字典与有限动态泛函]({p}research_note_673.md)'
         '原物理来源接到矩形手征字典，全部边界／字段共同输送；原H_b热矩保证全配置与Haar候选绝对可积。'
         '未证明一般正性、归一或原过程身份，近零权重的相对数值判据已撤回。'
         '两组、十八式通过，最新673／3253，1331份编号科学文件、2419份保护证据。'
         '[核验]({p}research_round_673_checks.json)、[全条件账]({p}unified_physics_condition_ledger_673.md)。'
         '完整正物理过程、连续及量子GR仍开放。')
order=('**当前执行顺序（673后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[674完整边界平均与原有序过程]({p}round674_drafts/STATUS.md)，'
       '核同一物理半区泛函与全部Gauss条件，不以仅电因子或抽样替代；目标不改。')
'''+publish[b:]
publish=publish.replace('独立规范历史、原手征权重与电动能的共同拼接条件','共同规范边界、变秩手征字典与有限动态泛函')
publish=publish.replace('独立规范历史与原电动能的共同限制','原物理来源与完整边界积分的共同存在条件')
publish=publish.replace('固定边界核与完整Gauss过程的边界','有限动态泛函与正物理过程的验收区别')
write('publish_round673.py',publish)
write('postcheck_round673.py',replace((HERE/'postcheck_round672.py').read_text('utf8'),mapping))
write('round673_drafts/research_note_673_draft.md',(HERE/'research_note_673.md').read_text('utf8'))
review='''673 primary-agent proof/code/scope review; no independent agent.
Previous turn completed672 and executed673 entry: progress. Current navigation,
state, report, saved results and process list checked. No running Python found.
Goal intact, no new task/automation/install/image checks. Old evidence preserved.
382-386,425,522-523 inherited through667. No eliminated assumptions restored.
2017 primary source covers topological-sector measure definitions; no new claim
to invent that framework, prove all auxiliary normalizations, or remove locality
conditions. New bridge retains original mass and physical source dictionary.
Rectangular Kl is r x rv; u/v sizes ru+rv=n=2r, r multiple32. Triangular
transform still determinant1. Reordering c,d has sign(-1)^(r*rv)=1.
Physical map has fixed n rows, L has rv+r columns. Odd ru auxiliary integral
vanishes for every pure physical-source coefficient. Fixture ranks are algebra
checks, not fabricated physical topological backgrounds.
At two AP time slices with identity temporal links, Re X=Ws positive. Xpsi=0
forces all spatial shifts fix psi, hence Cspsi=0 and invertible gamma0T0psi=0.
This supplies one nonzero analytic determinant value for every spatial q.
Connected compact full boundary group => Wilson zero set Haar-null. Borel
extension on null set harmless for integrals; derivative/locality limits not
covered. Direct Phi/Pfaffian norm bounds independent of small Wilson gap.
Two seams a,b transform g_i*a*g_j^-1 and g_j*b*g_i^-1. Set g_j=a gives
q_j,E_j,U_j transformed and Omega=a*b; transport physical sources contragredient.
Full Hb scalar heat kernel used in analytic candidate, including original
potential and non-diagonal electric terms. It is not an identification with HF.
Heat-kernel CS plus gauge covariance and669 finite configuration moments prove
global absolute integrability of each finite source coefficient.
Full Haar integral not evaluated numerically. No positivity/nonzero Z claim.
Strong-sample relative-reflection audit failed at weights with4 numerical zero
modes. Diagnostic retained; reliable weight phases/ratios not inferred there.
Initial general Hermiticity/reality assertions were insufficiently proved and
removed, with prereview note and versions retained. Regular-sample reflection
identity is a test, not an all-sector theorem. No hidden symmetrization.
2 groups,18 equations. Next674 must address actual averaged physical functional
and original ordered process, not another thermal parameter scan.
'''
for name in ('research_note_673.md','joint_gauss_boundary_functional.py','joint_gauss_boundary_functional_results.json','unified_physics_condition_ledger_673.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round673_drafts/final_review.txt',review)
print('673 preparation completed; scope audited')
