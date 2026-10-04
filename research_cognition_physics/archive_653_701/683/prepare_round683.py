"""Prepare683 artifacts; preserve older evidence and publication history."""
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
ledger=(HERE/'unified_physics_condition_ledger_682.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：682量子体涨落与原边界全来源匹配',
    '# 联合条件总账：683同时间片读口、共同极限与辅助反射')
ledger=ledger.replace('接[681全账](unified_physics_condition_ledger_681.md)，回填[682报告](research_note_682.md)。[结果](joint_bulk_fluctuation_matching_results.json)、[核验](research_round_682_checks.json)。',
    '接[682全账](unified_physics_condition_ledger_682.md)，回填[683报告](research_note_683.md)。[结果](joint_time_local_source_interface_results.json)、[核验](research_round_683_checks.json)。')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '683同时间片辅助读口以层数一致O(a)全来源误差接回原极限，并保实际局部Gauss合同；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('实际无限动力学存在|',
    '683关闭辅助观测跨相邻时间片的支撑问题；实际无限动力学存在|')
ledger=ledger.replace('682非零辅助族可保原泛函，但没有动态RP结论|',
    '682非零辅助族可保原泛函，但没有动态RP结论；683原补偿G在指定组件逐点反射下有精确自由障碍，不能自动套成熟Wilson体证明|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 683只在观测中用2I替代端点B；原约束和动能不变，质量用同一新观测。有限调节不完全等于原调节，但全部固定阶未归一来源在原完整平均下相差O(a)，常数一致于L而非盒体积。
- 解析误差不要求统一Wilson谱隙；677原a→0、aL→∞极限和669全阶热矩直接复用。不是物理时间格距趋零，不提供归一化或背景导数极限。
- 原局部颜色、弱和超荷变换使新读口、质量及来源字典精确协变；正半区支撑仅在辅助变量中成立，积分后的裸Phi仍可非局部。
- 原物理反射由679继承至极限，有限调节的反射差为O(a)，没有签收有限辅助RP。
- 678实际L=1补偿G=K†K虽严格正定，在指定组件逐点时间反射下不协变。原自由2×2盒精确差范数8a(1+2a)>0。这只是该辅助反射接法失败；全辅助积分的行列式抵消及原物理候选保持，不能据此宣判原Gauss/S9不RP。
- 旧空间、热参考及实际仪器条件直接继承679逐项表；原四分支仍不混同。反射适配、对合和正性仍是下一项实际接口。

## 本轮合并与下一项

683合并C01、C04、C14、C16、C17、C20、C22中的时间支撑、规范来源、质量和完整平均条件；关闭辅助读口跨相邻时间片的问题，保留明确极限范围，并排除一条不满足作用反射合同的直接接法。

原Q0正性、指定归一、H_F及记录身份、共同连续、量子GR和观测预测仍开放。局部读口和收敛不等于正物理过程；不能把单独补偿部门的障碍升级为物理子代数反例。

接[684](round684_drafts/STATUS.md)：检验适配原K的实际辅助反射是否为反线性对合，并同步全部来源、补偿和正半区。具体选择失败后考察可替代实现或直接物理代数，目标不改。
'''
write('unified_physics_condition_ledger_683.md',ledger)
mapping={'joint_bulk_fluctuation_matching':'joint_time_local_source_interface',
         '681':'682','682':'683','683':'684','3269':'3271','3271':'3273',
         '1355':'1358','1358':'1361','2525':'2539','2539':'2555'}
verify=replace((HERE/'verify_round682.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_gauge_flow_time_interface as prior_model',
                      'import joint_bulk_fluctuation_matching as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_683.md','round683_drafts/research_note_683_draft.md',
       'round683_drafts/final_review.txt','round684_drafts/STATUS.md']
names += ['round683_drafts/'+p for p in ('time_local_source_probe.py',
    'time_local_source_probe_results.json','time_local_source_entry.md','check_entry.py','entry_checks.json',
    'gauge_reflection_probe.py','gauge_reflection_probe_results.json',
    'compensation_free_probe.py','compensation_free_probe_results.json')]
assert len(names)==13
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace('Verify actual-body invisible deformations, complete sources and quantum-bulk controls.',
    'Verify strict-time sources, actual gauge covariance and scoped compensation reflection failure.')
write('verify_round683.py',verify)
pub=replace((HERE/'publish_round682.py').read_text('utf8'),mapping|{'## 328.':'## 329.','## 233.':'## 234.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第683轮完成：** [同时间片物理读口、完整来源极限与辅助反射边界]({p}research_note_683.md)'
         '原辅助观测可改为严格同时间片读口，规范协变且以层数一致O(a)误差恢复全部原未归一来源；'
         '指定补偿反射有精确自由障碍。两组、十六式通过，最新683／3273，1361份编号科学文件、2555份保护证据。'
         '[核验]({p}research_round_683_checks.json)、[全条件账]({p}unified_physics_condition_ledger_683.md)。'
         '原Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（683后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[684辅助反射的对合与共同正时间条件]({p}round684_drafts/STATUS.md)，'
       '核实际反射、全部来源及补偿，不把单独辅助障碍当原物理反例；目标不改。')
"""+pub[b:]
pub=pub.replace('量子体涨落、原边界权重与全部观测来源的共同匹配','同时间片物理读口、完整来源极限与辅助反射边界')
pub=pub.replace('同源体扩充仍须接实际物理时间','同时间片读口接回原来源，辅助反射仍需核验')
pub=pub.replace('旧空间合同复用，不把辅助涨落当物理噪声','旧空间合同复用，不把辅助反射障碍扩为物理反例')
write('publish_round683.py',pub)
write('postcheck_round683.py',replace((HERE/'postcheck_round682.py').read_text('utf8'),mapping))
write('round683_drafts/research_note_683_draft.md',(HERE/'research_note_683.md').read_text('utf8'))
review='''683 primary-agent proof/code/scope review; no independent agent.
Prior goal turn published682 and executed683 entry: progress, not wait.
Authoritative navigation/report/results checked; no Python process found by Get-Process.
CIM enumeration denied; no old live handle restarted or terminated.
Strict observation modifies only endpoints B to2I, and mass together, at finite regulator.
Endpoint constraint solved exactly; norm delta uniform in L and Wilson gap.
Auxiliary parameter a is not physical lattice spacing; effective Phi remains nonlocal.
Finite fixed-order Pfaffian telescoping and669/673 moments give full-average O(a).
677 original limit reused after auxiliary integration; no volume-uniform or normalized claim.
Actual local gauge action checked for observation, effective dictionary and complete mass.
Original679 exact physical reflection yields an integrated O(a) defect for new regulator.
Small source diagnostics never divided by physical weight to infer positivity.
Original L1 compensation Gram formula computed from K, not assumed equal to Wilson5d.
Exact free original box proves norm8a(1+2a) obstruction for specified scalar reflection.
This is unintegrated action noninvariance, not full Gauss/S9 physical RP counterexample.
Compensation stays positive definite and its integrated determinant ratio remains1.
Primary KikukawaUsui1005.3751v3 sectionsII/VI read; no image inspection.
Alternative reflections, anti-linear square, joint sources and physical subalgebra remain open.
382-386/425/522-523 reused with no extra Lipschitz or stacked coordinate assumptions.
OriginalHF, fixed continuous matter, given classical geometry, chiral candidate remain distinct.
No whole-goal completion, new cognitive axiom, physical dimension or quantumGR claim.
Two coherent groups; sixteen equations; next684 entry tests actual alternative reflection.
'''
for name in ('research_note_683.md','joint_time_local_source_interface.py',
             'joint_time_local_source_interface_results.json','unified_physics_condition_ledger_683.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round683_drafts/final_review.txt',review)
print('683 preparation completed')
