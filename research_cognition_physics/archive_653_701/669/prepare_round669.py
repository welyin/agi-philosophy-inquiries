"""Prepare669 scientific ledger, review and publication; preserve all history."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_668.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：668原质量、物理观测与辅助正泛函',
    '# 联合条件总账：669原动态标量、手征物质与规范过程')
ledger=ledger.replace('接[667全账](unified_physics_condition_ledger_667.md)，回填[668报告](research_note_668.md)。[结果](joint_mass_auxiliary_reflection_results.json)、[核验](research_round_668_checks.json)。',
    '接[668全账](unified_physics_condition_ledger_668.md)，回填[669报告](research_note_669.md)。[结果](joint_dynamic_scalar_gauge_interface_results.json)、[核验](research_round_669_checks.json)。')
ledger=ledger.replace('668冻结标量不足以承载原s读口；真实共同instrument和动态引力反馈',
    '669原动态s函数可积，但真实instrument仍需同一H、条件后态和规范过程')
ledger=ledger.replace('668已接固定原质量与正归一；动态标量和一般规范相互作用测度仍缺，参考排序不是连续反常',
    '668固定质量已归一；669动态标量候选可积且未归一正，指定归一及一般规范测度仍缺')
ledger=ledger.replace('实际动态质量、机制选择、一般背景极点及现实参数',
    '669实际动态质量插入可积；同一规范传播、机制选择、一般背景极点及现实参数')
ledger=ledger.replace('668有限正泛函不等于原完整Gibbs；实际热态尺度极限及宇宙初态，不要求空Fock全局酉',
    '669已接原H_b配置热矩，动态归一多项式的指定零点仍待排除；原H_F态识别、尺度极限及宇宙初态仍缺')
ledger=ledger.replace('群选择及连续动力学|','669排除固定自由费米动能靠标量积分恢复局部Gauss；群选择及真实规范动力学仍缺|')
ledger=ledger.replace('667热迹、668质量权重／来源可排除指定错误改写',
    '667—668字典检查；669原质量与传播不同步会改变规范相关预测')
a=ledger.index('## 本轮合并与下一项')
ledger=ledger[:a]+'''## 本轮合并与下一项

669把C19原热参考的能量控制推进为实际配置矩，再接C17原动态质量及C03原读口的有限多时刻可积性。候选动态质量的全多项式未归一反射正性成立，归一是实非负多项式且除有限可能零点外严格正；不能假定指定λ=1已避开零点，亦不以改选弱耦合宣称统一完成。

更强的联合限制在C14／C15／C17：原玻色动态积分给规范不变相互作用乘子，不能修复固定自由费米动能的局部Gauss缺陷。该具体方案对所有λ均不满足所需局部规范合同。质量与传播必须实际使用同一链路；纯规范输送已复算，但非平坦动态规范测度、正性和原H_F映射尚未成立。

603、623、624、643原完整热／域／记录／有序历史全部继承，不重新构造路径积分。605—610已有热态全支撑、投影、电动能与局域组合条件保留，不能将常数自由证明用于一般规范背景。612时谱、653荷谱、553运行关系以及上表空间替代路线仍有效。

接[670](round670_drafts/STATUS.md)：检验真实非平坦规范链路中的手征核、原S⁹配对、物理观测及全部质量的共同协变字典，再核动态积分、反射和记录。同一Euclidean权重不代替221—222、230、506所需完整过程。认知设计后置；四分支尚未统一，目标不改。
'''
write('unified_physics_condition_ledger_669.md',ledger)
mapping={'joint_mass_auxiliary_reflection':'joint_dynamic_scalar_gauge_interface',
    '667':'668','668':'669','669':'670','3241':'3243','3243':'3245',
    '1313':'1316','1316':'1319','2342':'2354','2354':'2364'}
verify=replace((HERE/'verify_round668.py').read_text('utf8'),mapping)
verify=verify.replace('Verify actual original mass reflection, strict finite normalization and common sources.',
    'Verify actual scalar moment interface and joint mass/kinetic gauge requirement.')
verify=verify.replace('import joint_reference_spatial_process as prior_model','import joint_mass_auxiliary_reflection as prior_model')
verify=verify.replace('mass_observation_probe','dynamic_scalar_probe')
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==18")
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_669.md','round669_drafts/research_note_669_draft.md',
           'round669_drafts/final_review.txt','round670_drafts/STATUS.md',
           'round669_drafts/dynamic_scalar_entry.md','round669_drafts/dynamic_scalar_probe.py',
           'round669_drafts/dynamic_scalar_probe_results.json')
"""+verify[b:]
write('verify_round669.py',verify)
publish=replace((HERE/'publish_round668.py').read_text('utf8'),mapping|{'## 314.':'## 315.','## 219.':'## 220.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第669轮完成：** [原动态标量、手征物质与共同规范过程的条件]({p}research_note_669.md)'
         '原热过程控制全部有限时刻标量多项式；动态质量候选可积并有未归一反射正性。'
         '保自由费米动能只积分原玻色变量，不能恢复局部Gauss；质量和链路须共同匹配。'
         '两组、十八式通过，最新669／3245，1319份编号科学文件、2364份保护证据。'
         '[核验]({p}research_round_669_checks.json)、[全条件账]({p}unified_physics_condition_ledger_669.md)。'
         '指定归一、完整规范过程、连续及量子GR仍开放。')
order=('**当前执行顺序（669后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[670真实规范链路与手征物质]({p}round670_drafts/STATUS.md)，'
       '核非平坦背景的测度、质量和观测共同字典，再接动态正性与记录；目标不改。')
"""+publish[b:]
publish=publish.replace('原非零质量、物理观测与辅助测度的共同反射结构','原动态标量、手征物质与共同规范过程的条件')
publish=publish.replace('原质量与物理观测的共同正泛函','原动态标量与规范传播的联合条件')
publish=publish.replace('固定原质量的反射接口与动态标量缺口','原热标量可积不替代共同规范动力学')
write('publish_round669.py',publish)
write('postcheck_round669.py',replace((HERE/'postcheck_round668.py').read_text('utf8'),mapping))
write('round669_drafts/research_note_669_draft.md',(HERE/'research_note_669.md').read_text('utf8'))
review='''669 primary-agent proof/code audit; no independent agent review.
Previous goal turn completed668 and executed669 entry: progress. Read current
navigation/note/results/postcheck; no live Python. No goal edits, new tasks,
automation, images or subagents. Existing Python/NumPy only.
Original603/623/624/643 thermal/domain/history results reused, not rederived as
new research. New all-order configuration moment bridge is not inferred merely
from energy moments. Actual H5 gradient norm, quartic coercivity and fixed graph
positive coefficients give Gamma(W)<=C(1+W)^2. Weighted eigenfunction IMS with
compact cutoffs, inductively integrable tail derivative, gives all-order moment
polynomial. Scalar Hb has all positive heat traces; summation valid. Gauss W
moments restricted only to invariant W; no charged insertion treated physical.
Exact thermal loop measure retains original bosonic coordinates and links.
Stationary marginals plus Holder proves multi-time polynomial integrability,
not just a static potential integral. All masses linear in original x=phi/sqrtF.
Candidate adds free-Weyl/S9 to actual Hb loop and explicitly differs from H_F.
Finite Grassmann expansion integrable; even local reflection split gives full
unnormalized RP. Z real nonnegative polynomial, Z(0)=1; finitely many possible
real zeros remain. No all-lambda strictness or lambda=1 completion asserted.
Original scalar record bounded functions included but not actual instruments.
For time-independent spatial local gauge, original boson measure and local mass
give invariant interaction multiplier with body1. Its invertibility means it
cannot cancel noncovariance of fixed free physical kinetic kernel. Free AP
Weyl coordinate invertibility checked via symbol; total Jacobian is unity.
Specific candidate rejected for local Gauss at all lambda, not entire project.
Numerics actual full16 internal mass and local colour/weak/hypercharge. Pure
gauge simultaneous mass/kinetic transport recovers weight; not general gauge RP.
Potential sampling/recurrence diagnostic not claimed original spectrum/path
integral. Nonconstant mass entry replayed and not recounted as new science.
605-610 electric/projector/thermal limitations,612/653/553 and spatial old
interfaces preserved. Next670 actual nonflat gauge dictionary, no repeated
free spectrum, pure gauge-only round, design or precision optimization.
Two fresh groups,18 equations; unified full goal remains active.
'''
for name in ('research_note_669.md','joint_dynamic_scalar_gauge_interface.py',
             'joint_dynamic_scalar_gauge_interface_results.json','unified_physics_condition_ledger_669.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round669_drafts/final_review.txt',review)
print('669 report/review frozen; verification and publication prepared')
