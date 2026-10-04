"""Prepare682 publication, preserving all historical evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def replace(text,mapping):
    keys=sorted(mapping,key=len,reverse=True)
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in keys]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)
ledger=(HERE/'unified_physics_condition_ledger_681.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：681规范流、物理时间与正性接口',
    '# 联合条件总账：682量子体涨落与原边界全来源匹配')
ledger=ledger.replace('接[680全账](unified_physics_condition_ledger_680.md)，回填[681报告](research_note_681.md)。[结果](joint_gauge_flow_time_interface_results.json)、[核验](research_round_681_checks.json)。',
    '接[681全账](unified_physics_condition_ledger_681.md)，回填[682报告](research_note_682.md)。[结果](joint_bulk_fluctuation_matching_results.json)、[核验](research_round_682_checks.json)。')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '682原有限辅助体在OΣOᵀ=0下精确保全物理来源，允许非退化Berezin逆块；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('本原候选|',
    '本原候选；682非零辅助族可保原泛函，但没有动态RP结论|')
ledger=ledger.replace('体行列式与观测字典的有限背景响应须共同抵消；',
    '体行列式与观测字典的有限背景响应须共同抵消；682可见辅助变形还产生独立来源接触块，不能只配权重和作用；')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 682在678实际局部体中声明η二次配对，Σ=K⁻¹RK⁻ᵀ。保持所有关于Oz的条件多项式等价于OΣOᵀ=0；不把它说成完整Gauss积分偶然相等的必要条件。
- 原物理z来源TO中T可逆，奇异值1/2及1/4，故非零可见收缩产生非零物理来源接触。不能只凭体行列式匹配签收原观测。
- L=1的原局部规范协变不可见族包含可逆Σ，且OΣ非零、OΣOᵀ=0；所有有限质量与全部物理来源逐背景保持，原补偿及全部平均仍相同。无需坚持整个zz逆块为零。
- Σ为Berezin反对称逆块，不是正经典噪声；辅助表示自由不意味着正物理扩充。原泛函若有负方向，精确保它的任何表示也保该负方向。
- 完整正两层量子体对照给条件均值的RP负方向；补回条件涨落须同步边界Schur作用和来源接触。另行固定边界协方差的明确非局部补偿虽Euclidean正，仍可不RP。该对照不用于宣判原Grassmann/Gauss模型。
- 旧空间、热参考、仪器与四个原分支范围继承681及679表；没有新增维数定理或恢复已消去假设。

## 本轮合并与下一项

682联合连接C01、C02、C04、C14、C16、C20、C22：体扩充需共同匹配边界作用、所有来源、接触和权重。原辅助逆块可以非零／非退化，条件性原物理字典仍能精确保留；真正正时间实现没有因此自动获得。

原Q0正性、指定归一、H_F及记录身份、共同连续、量子GR和观测预测仍开放。不能把不可见表示变形当作修复原物理正性，更不能把另外的正自由模型拼成全部统一结果。

接[683](round683_drafts/STATUS.md)：检验原同源辅助扩充与补偿部门的真实反射、时间支撑及物理子代数，明确完整辅助RP作为充分路线的条件。保留可替代表示，不扫描δ或重复一般消元；目标不改。
'''
write('unified_physics_condition_ledger_682.md',ledger)
mapping={'joint_gauge_flow_time_interface':'joint_bulk_fluctuation_matching',
         '680':'681','681':'682','682':'683',
         '3267':'3269','3269':'3271','1352':'1355','1355':'1358',
         '2513':'2525','2525':'2539'}
verify=replace((HERE/'verify_round681.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_mass_reflection_congruence as prior_model',
                      'import joint_gauge_flow_time_interface as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_682.md','round682_drafts/research_note_682_draft.md',
       'round682_drafts/final_review.txt','round683_drafts/STATUS.md']
names+=['round682_drafts/'+p for p in ('conditional_bulk_probe.py','conditional_bulk_probe_results.json',
    'conditional_bulk_entry.md','check_entry.py','entry_checks.json',
    'body_fluctuation_probe.py','body_fluctuation_probe_results.json')]
assert len(names)==11
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==16")
verify=verify.replace('Verify physical gauge-flow time-interface counterexamples and their restricted scope.',
                     'Verify actual-body invisible deformations, complete sources and quantum-bulk controls.')
write('verify_round682.py',verify)
pub=replace((HERE/'publish_round681.py').read_text('utf8'),
            mapping|{'## 327.':'## 328.','## 232.':'## 233.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第682轮完成：** [量子体涨落、原边界权重与全部观测来源的共同匹配]({p}research_note_682.md)'
         '原678局部规范协变体允许非退化辅助Berezin逆块，并在可见收缩为零时精确保全部质量与来源；'
         '可见变形须另配作用、观测与接触项，体行列式相同仍不够。两组、十六式通过，最新682／3271，1358份编号科学文件、2539份保护证据。'
         '[核验]({p}research_round_682_checks.json)、[全条件账]({p}unified_physics_condition_ledger_682.md)。'
         '辅助表示不自动给RP；原正性、H_F身份、连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（682后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[683同源辅助扩充与真正的正时间实现]({p}round683_drafts/STATUS.md)，'
       '核原物理子代数、辅助反射和补偿，不以表示恒等式签收正性；目标不改。')
"""+pub[b:]
pub=pub.replace('规范流、物理时间与共同正性接口的边界','量子体涨落、原边界权重与全部观测来源的共同匹配')
pub=pub.replace('规范流不自动给共同物理时间','同源体扩充仍须接实际物理时间')
pub=pub.replace('旧空间合同复用，流后观测须核正性','旧空间合同复用，不把辅助涨落当物理噪声')
write('publish_round682.py',pub)
write('postcheck_round682.py',replace((HERE/'postcheck_round681.py').read_text('utf8'),mapping))
write('round682_drafts/research_note_682_draft.md',(HERE/'research_note_682.md').read_text('utf8'))
review='''682 primary-agent proof/code/scope review; no independent agent.
Prior turn published681 and executed682 entry: progress.
Navigation,681 report/results/postcheck,682 entry and processes inspected.
643/661/673/675/677-681 reused. No repeated space proofs or recovered assumptions.
Real Gaussian two-mode full quantum process is a declared control, not original Q0.
Conditional mean loses RP; true conditional covariance restores full original process.
Pinning a different boundary covariance requires another nonlocal counterterm.
Its explicit physical bulk RP witness matches the rational double-pole criterion;
read primary Arici et al1712.04308 Theorem3.7/Proposition3.6, no figure inspection.
Original678 added eta-eta R; K invertible, Sigma antisymmetric, not positive noise.
Berezin moment generating formula fixes conditional all-polynomial iff O Sigma O^T=0.
This is not necessity for accidental equality after physical/Gauss averaging.
Original physical T invertible, singular values1/2,1/4; visible contact cannot hide in T.
Full massive invisible identity proved polynomially, not by an inverse physical weight.
Two L1 families include full-rank mixed Sigma with O Sigma nonzero but two-leg zero.
Gauge covariance of K/O/Sigma/R checked under all original local gauge factors.
Finite physical range retained; no original Wilson hard gap or new physical particle.
Invisible family keeps old determinant compensation and all complete source averages.
Visible deformation changes N/Phi/contact, including nonzero two-source omission error.
Actual1536/1538-dimensional Pf compared to augmented Schur with original complex phase.
Small coefficients used only in numerical identity ratios, not physical normalization.
Source response must differentiate R/Sigma/background together.
No RP, originalHF, chiral continuum, quantumGR or full-goal completion claimed.
Next683 actual reflection and physical-algebra support; no parameter scanning.
Two groups16 equations; original history/entries retained.
'''
for name in ('research_note_682.md','joint_bulk_fluctuation_matching.py',
             'joint_bulk_fluctuation_matching_results.json','unified_physics_condition_ledger_682.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round682_drafts/final_review.txt',review)
print('682 preparation completed')

