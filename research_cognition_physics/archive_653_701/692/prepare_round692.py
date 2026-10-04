"""Prepare692 publication, preserving all691 and692-entry evidence."""
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
ledger=(HERE/'unified_physics_condition_ledger_691.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：691中性酉对称与完整正性的分次约化','# 联合条件总账：692原辅助轨道与全部物理来源的精确有限求积')
ledger=ledger.replace('接[690全账](unified_physics_condition_ledger_690.md)，回填[691报告](research_note_691.md)。[结果](joint_neutral_symmetry_reduction_results.json)、[核验](research_round_691_checks.json)。',
    '接[691全账](unified_physics_condition_ledger_691.md)，回填[692报告](research_note_692.md)。[结果](joint_auxiliary_orbit_quadrature_results.json)、[核验](research_round_692_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 692完整原平均中的辅助轨道与有限求积

- 原10维实向量在G下分6维颜色和4维弱部分。局部轨道由r=|P_c E|²标记，原球面边缘为Beta(3,2)，没有换测度或扩大群。
- 对全部Gauss不变、E独立物理插入F，在保留完整q、原H_b和双Haar平均的Psi_F中可同时转动每节点E与所有原字段／边界，故辅助角度消去。固定q端点的核不具备这个约化。
- 673／675同一辅助因子每节点E次数至多32；全平均局部群不变多项式仅含两个平方长度，限制球面后每个r次数至多16。无须对原128变量最高项重新逐项枚举。
- 同一九点Beta正权Gaussian规则精确保全部上述物理来源与归一；当前四节点为6561项完整规范／玻色泛函。该规则同时适用原有限质量，未删除辅助16通道，也未输入九态物理系统。
- 固定背景原128维Pfaffian的符号平均径向积分已独立复算；完整512维B来源沿局部规范轨道协变。未数值计算所有6561项或完整H_b／Haar积分。
- 正求积权不是RP证书；任意显含E的680扩充测试代数不受统一九点界覆盖，原合同仍保留。E不成为s记录；原四分支、612／646及旧空间结果不变。

## 本轮合并与下一项

C01物理来源、C14 Gauss、C19归一获得同一非扰动辅助消元规则，保留原模型参数和物理变量。尚缺完整剩余规范平均的符号、原H_F及共同尺度／量子几何。

接[693](round693_drafts/STATUS.md)：优先核剩余完整规范来源的符号控制或精确分块。不得把固定E核当已完成全平均；不同边界坐标选择须同步输送E，不无条件叠加约化。目标不改。
'''
write('unified_physics_condition_ledger_692.md',ledger)
mapping={'joint_neutral_symmetry_reduction':'joint_auxiliary_orbit_quadrature','690':'691','691':'692','692':'693',
         '3287':'3289','3289':'3291','1382':'1385','1385':'1388','2628':'2641','2641':'2655'}
verify=replace((HERE/'verify_round691.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_spatial_charge_dictionary as prior_model','import joint_neutral_symmetry_reduction as prior_model')
verify=verify.replace('Verify692 actual neutral symmetry, graded RP reduction and odd Gauss probe.','Verify692 original auxiliary orbit reduction and complete-source quadrature.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=('unified_physics_condition_ledger_692.md','round692_drafts/research_note_692_draft.md',
       'round692_drafts/final_review.txt','round693_drafts/STATUS.md',
       'round692_drafts/odd_gauss_source_probe.py','round692_drafts/odd_gauss_source_probe_results.json',
       'round692_drafts/odd_gauss_source_entry.md','round692_drafts/check_entry.py',
       'round692_drafts/entry_checks.json','round692_drafts/publish_entry.py','round692_drafts/navigation_checks.json')
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==14")
write('verify_round692.py',verify)
pub=replace((HERE/'publish_round691.py').read_text('utf8'),mapping|{'## 337.':'## 338.','## 242.':'## 243.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第692轮完成：** [完整物理泛函的辅助轨道与九点精确求积]({p}research_note_692.md)'
         '保留原H_b、全配置及双Gauss平均，对全部E独立Gauss物理来源，每节点S9精确化为九点正权求和。'
         '两组、十四式通过，最新692／3291，1388份编号科学文件、2655份保护证据。'
         '[核验]({p}research_round_692_checks.json)、[全条件账]({p}unified_physics_condition_ledger_692.md)。'
         '不替代固定背景球面积分，不凭求积正权领取RP；原H_F、共同连续及量子GR仍开放。')
order=('**当前执行顺序（692后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[693辅助消元后的原规范来源]({p}round693_drafts/STATUS.md)，'
       '优先核完整平均的符号控制或精确分块，不重复纯时间链与无控制抽样；目标不改。')
"""+pub[b:]
pub=pub.replace('中性相位、完整正性约化与奇Gauss测试','完整物理泛函的辅助轨道与九点精确求积')
pub=pub.replace('完整物理正性的中性分离与奇Gauss条件','完整物理泛函的辅助轨道与精确求积')
pub=pub.replace('旧空间合同复用，酉对称不冒充测量仪器','旧空间合同复用，辅助轨道代表不是物理参考')
write('publish_round692.py',pub)
write('postcheck_round692.py',replace((HERE/'postcheck_round691.py').read_text('utf8'),mapping))
with (HERE/'round692_drafts/research_note_692_draft.md').open('xb') as f:
    f.write((HERE/'research_note_692.md').read_bytes())
review='''692 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed691 and verified/published692 source entry: progress.
Read navigation,691 science/results and692 entry. No live Python process found.
Reuse614/653 originalG and real10 vector;657/675 sphere polynomials;673 full source separation.
OriginalG contains independently transitive SU3 and SU2 on unit C3 and C2.
Center(I3,I2,-1) produces rank6/rank4 real projectors in actualT representation.
Beta(3,2) follows6+4 isotropic Gaussian norm squares; normalized12r²(1-r).
All local E rotations absorbed only with complete q, original Hb, both Haar boundaries and Gauss-invariant F.
No fixed-q boundary kernel equivalence claimed; no independent deletion of another gauge seam after fixing E axes.
Local64 Grassmann row space bounds auxiliary E degree32; unchanged by arbitrary pure physical source order.
Off-sphere polynomial extension defined using the factored auxiliary expression, not unproved off-sphere N identity.
Complete-average local invariant polynomial generated by two squared norms, so radial degree<=16 each site.
Nine-node Gaussian quadrature exact through17; exact rational polynomial orthogonality and analytic positive weights.
Actual full512-dimensional baryon source covariance checked; actual128-mode sign-symmetrized radial Pfaffian integrated two ways.
All16 auxiliary channels/Hb/doubleHaar/normalization retained. Current four-site sum6561 is not evaluated in full.
Only E-independent Gauss physical sources covered; arbitrary E-dependent680 extended algebra kept separate.
Positive numerical cubature weights are not a new state positivity theorem. General RP/HF/continuum/GR still open.
Old382-386/425/522-523 conditions inherited; no restored384 Lipschitz;386/425 alternative;E not s.
Two new groups/14 displayed equations. Next693 full remaining gauge average, no precision-only repetition.
'''
for name in ('research_note_692.md','joint_auxiliary_orbit_quadrature.py','joint_auxiliary_orbit_quadrature_results.json','unified_physics_condition_ledger_692.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round692_drafts/final_review.txt',review)
print('692 publication preparation completed')
