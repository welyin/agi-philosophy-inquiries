"""Prepare681 publication without changing any frozen evidence."""
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

ledger=(HERE/'unified_physics_condition_ledger_680.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：680完整质量的正性约化与严格归一边界',
    '# 联合条件总账：681规范流、物理时间与正性接口')
ledger=ledger.replace('接[679全账](unified_physics_condition_ledger_679.md)，回填[680报告](research_note_680.md)。[结果](joint_mass_reflection_congruence_results.json)、[核验](research_round_680_checks.json)。',
    '接[680全账](unified_physics_condition_ledger_680.md)，回填[681报告](research_note_681.md)。[结果](joint_gauge_flow_time_interface_results.json)、[核验](research_round_681_checks.json)。')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '681排除仅反射协变就把正流时字段当作局部物理观测的接法；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('实际无限动力学存在|',
    '681辅助流的时间标签不自动保正时间代数；实际无限动力学存在|')
ledger=ledger.replace('该动态正性、指定非零归一与连续局域性仍缺|',
    '该动态正性、指定非零归一与连续局域性仍缺；681自由Abelian流后物理模式反例不判定本原候选|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 681对成熟规范流补上实际物理时间检验：完整Gaussian规范不变模式在任意正时空流后有显式负反射范数。是对指定观测映射的反例，不是原全Gauss/S9对象或边界手征理论的反例。
- 梯度流及衰减半空间EOM流满足同一有限谱矩判据；第二分支不等于文献完整有限slab边界条件。原论文背景场范围及动态行列式仍单列。
- 仅空间平滑保已有正时间代数和RP，不能制造原Q0正性；它不衰减空间均匀、时间变化的电历史，故不能直接替代原镜像背景衰减任务。
- 675的辅助E与实际物理观测继续区分。680允许E系数的大代数不能自动视为原s记录；E无关的q,Y子代数也对质量指数及逆封闭，故质量RP双向等价在真正物理子代数仍成立。
- 382—386、425、522—523按679表及本账直接复用；本轮给定空间作反例，并未新增维数或光滑坐标结论。

## 本轮合并与下一项

681联合检验C01、C04、C14、C16、C20：规范流、时间标签与正物理内积必须在实际观测映射上共同匹配。排除“反射协变即保RP”的直接搬用，保留仅空间处理和物理边界来源作为不同分支。未新增认知公理，未关闭原无质量动态RP。

四个原模型分支、指定归一、原H_F及记录身份、共同连续极限、量子GR和现实预测仍开放。条件账没有将自由规范诊断代替原全部物质积分。

接[682](round682_drafts/STATUS.md)：核完整量子体的条件涨落和全部带源Schur项，对比678确定性约束体。先给同一正过程内的精确检验，再判断原权重及来源能否匹配；不以一般Schur公式或另一个正模型签收原对象。
'''
write('unified_physics_condition_ledger_681.md',ledger)
mapping={'joint_mass_reflection_congruence':'joint_gauge_flow_time_interface',
         '679':'680','680':'681','681':'682',
         '3265':'3267','3267':'3269','1349':'1352','1352':'1355',
         '2501':'2513','2513':'2525'}
verify=replace((HERE/'verify_round680.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_physical_source_reflection as prior_model',
                      'import joint_mass_reflection_congruence as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_681.md','round681_drafts/research_note_681_draft.md',
       'round681_drafts/final_review.txt','round682_drafts/STATUS.md']
names+=['round681_drafts/'+s for s in ('flow_half_support_probe.py','flow_half_support_probe_results.json',
       'mature_flow_mapping_entry.md','check_entry.py','entry_checks.json')]
assert len(names)==9
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace('Verify full original mass congruence and exact normalization boundary.',
                     'Verify physical gauge-flow time-interface counterexamples and their restricted scope.')
write('verify_round681.py',verify)
pub=replace((HERE/'publish_round680.py').read_text('utf8'),
            mapping|{'## 326.':'## 327.','## 231.':'## 232.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第681轮完成：** [规范流、物理时间与共同正性接口的边界]({p}research_note_681.md)'
         '完整自由规范物理模式给任意正时空流的负反射范数；反射协变不足以把流后字段当作局部物理观测。'
         '仅空间流保已有正性，但不完成原电背景衰减。两组、十四式通过，最新681／3269，1355份编号科学文件、2525份保护证据。'
         '[核验]({p}research_round_681_checks.json)、[全条件账]({p}unified_physics_condition_ledger_681.md)。'
         '本轮不判定原全Gauss/S9候选；H_F身份、连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（681后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[682完整量子体、条件涨落与原边界来源]({p}round682_drafts/STATUS.md)，'
       '核真实体积分的涨落、权重及全部来源，和678约束体分别记账；目标不改。')
"""+pub[b:]
pub=pub.replace('完整物理代数中的质量约化与严格归一边界','规范流、物理时间与共同正性接口的边界')
pub=pub.replace('质量约化不替代原过程与空间接口','规范流不自动给共同物理时间')
pub=pub.replace('旧方向与热参考合同保持，归一门槛独立','旧空间合同复用，流后观测须核正性')
write('publish_round681.py',pub)
write('postcheck_round681.py',replace((HERE/'postcheck_round680.py').read_text('utf8'),mapping))
write('round681_drafts/research_note_681_draft.md',(HERE/'research_note_681.md').read_text('utf8'))
review='''681 primary-agent proof/code/scope review; no independent agent.
Previous reminder response was no new research evidence; not counted as a round.
Read four navigation heads,680 note/results,681 executed entry and historical scope.
Process command-line query denied; Get-Process fallback found no Python session.
Read primary sources1710.11618,1005.3751,2606.05306v2,2211.15176v2.
No claimed universal absence of a dynamic-gauge theorem from search results.
Countermodel: declared free Abelian transverse magnetic Fourier mode, theta-even.
Full Gaussian measure integrated. Not fixed E or original interacting Gauss/S9.
Analytic cosine Taylor remainder yields coefficient217/12, checked arithmetically.
m2 lower bound from[omega,2omega]; m4 upper bound via Gaussian/exponential tails.
epsilon squared6L2/(217U4) gives strict negative upper bound for any positive flow.
EOM branch is decaying half-space only, not whole finite-slab article.
Closed gradient formula includes spatial factor cancellation, checked by quadrature.
Finite periodic covariance positive and reflection commuting, but RP Gram negative.
Spatial smearing inherits existing RP for actual positive-time algebra and all moments.
Spatially uniform electric history remains; no automatic mirror-decoupling transfer.
No rejection of boundary-only slab theory, original Q0, HF, or cognitive principles.
Physical E-independent algebra distinction675 preserved;680 mass map closes there too.
No new space assumptions:382-386/425/522-523 reused; no Lipschitz reinstatement.
682 will compare full conditional fluctuations/sources, not assume new model identity.
Two groups14 display equations; original all evidence retained, goal unchanged.
'''
for name in ('research_note_681.md','joint_gauge_flow_time_interface.py',
             'joint_gauge_flow_time_interface_results.json','unified_physics_condition_ledger_681.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round681_drafts/final_review.txt',review)
print('681 preparation completed')

