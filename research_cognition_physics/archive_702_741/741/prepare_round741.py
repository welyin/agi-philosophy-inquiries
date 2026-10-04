"""Prepare the absolute-source first-order bridge and immutable evidence."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def remap(text,values):
    return re.sub('|'.join(re.escape(k) for k in sorted(values,key=len,reverse=True)),lambda m:values[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_740.md').read_text('utf8')
ledger='# 联合条件总账：741实际绝对来源与一阶共同反作用\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[739全账](unified_physics_condition_ledger_739.md)，回填[740报告](research_note_740.md)。[结果](joint_causal_eft_branch_results.json)、[核验](research_round_740_checks.json)。',
    '接[740全账](unified_physics_condition_ledger_740.md)，回填[741报告](research_note_741.md)。[结果](joint_absolute_source_development_results.json)、[核验](research_round_741_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**741当前增量：** 在原短时经典发展、固定参考／记录及735守恒处方下，完整绝对来源、全部线性场方程、初始约束与变背景Hadamard记录族共同构造，完整方程及约束有O(λ²)残差。首阶不依赖保持零阶准备的一阶辅助扩展，二阶仍依赖。精确解、λ=1误差、内部准备和物理连续未签收。\n\n'+marker)
updates={'C04':'741原完整部门的绝对首阶反作用及同一态族已有条件性构造',
         'C19':'741保原参考与记录；零阶准备仍影响首阶，辅助一阶扩展影响二阶',
         'C22':'741联合约束与完整绝对方程残差O(λ²)；精确半经典解及自主装置开放'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 741绝对首阶共同构造

- 原完整绝对源保真空、记录差及固定有限项，采用735联合Ward；不再仅靠同背景相对相消。
- 731原初始右逆与732全部约化主部直接复用，包含全部规范、五标量和几何。
- 所得光滑背景修正可接共同过去、正自对偶Hadamard参考及同一非选择记录；不暗中冻结后态。
- 完整变背景来源进入同一有限阶方程，方程和初始约束余项O(λ²)。只给已构造路径的有限阶界，不借未证的全源Banach适定。
- 零阶准备固定时首阶不依赖一阶辅助扩展；二阶历史差实际非零。参考选择和真实制备仍须交代。
- 740重求和增长极点和实际误差边界保留；不以本轮宣布精确非线性完成。

## 本轮合并与下一项

C04／C19／C22在一份条件性半经典近似族中连接。内部准备、量子约束、共同连续及唯一性仍开放。

接[742](round742_drafts/STATUS.md)：复用524、575、593、620、633—634及716，核现有记录的局域过程和共同来源条件，不新造仪器替代原对象。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_741.md',ledger)
write('round741_drafts/research_note_741_draft.md',(HERE/'research_note_741.md').read_text('utf8'))
write('round742_drafts/STATUS.md','''# 第742轮入口：现有记录的局域过程与共同来源

接[741](../research_note_741.md)、[条件账](../unified_physics_condition_ledger_741.md)。绝对一阶反作用和同一态族已条件性接通；数学参考与记录后态不等于实际内部准备过程。

1. 回查524系统—探针、575自主经典探针、593／620完整来源交换、633／634连续sterile记录及因果／界面限制、716 instrument和后续来源域。已有注能、类空传讯反例与记录噪声不重复计轮次。
2. 对接成熟局域测量方案的明确前提，核原读口、允许相互作用、末端记录和总Ward是否能共用；禁止只凭CP或可测就签收。
3. 先整合现有模型条件、指出尚缺同一作用或来源项，不优先发明新物种、新开关或全新装置。
4. 若必须扩展，明确区分原模型内实现、不同有效过程及额外建模输入；不把替代过程冒充原瞬时instrument的精确实现。
5. 准备、控制、物质与几何反作用一同记账；数学辅助过去不作为免费物理资源。
6. 不重复741一般有限阶Taylor展开或继续脉冲精度扫描；旧空间382—386、425、522—523，604及649／699与统一目标保持。
''')
write('round741_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/1311.7661',locations='Section2.4, proposition2.9 and smooth background families',
    use='Retarded comparison on compactly changed backgrounds; original matrix and reality structure mapped through730/734.'),dict(
    url='https://arxiv.org/html/1210.4031',locations='Proposition3.6 and parametrix construction',
    use='Smooth subtracted local sources; parameter regularity uses the jointly prepared family, not individually Hadamard states alone.')],
    new='Absolute source, full original linearized dynamics/initial constraints and compatible parameterized record-state family yield a common finite-order residual.',
    inherited='573 classical development;731 source right inverse;732 relative evolution;734 parameterized history and contacts;735 chosen conserving prescription.',
    excluded='Exact semiclassical solution, error relative to exact inverse, physical lambda=1 control, internal preparation, full quantum Gauss, original graph continuum.'),ensure_ascii=False,indent=2)+'\n')
write('round741_drafts/scope_and_dedup_review.md','''# 741范围与证明审查

上一目标轮完成740且执行741入口，属进展。本次按最新文件接续；详细进程接口拒绝访问，普通进程查询未见活跃Python。主代理审查，无子代理或图像检验。

- 731的任意光滑可容许源右逆直接使用；原颜色不变准备与singlet记录保证颜色源为零，未加入常颜色荷。
- 732只有相对源；735补绝对Ward；本轮新增真正变背景来源及同一记录族，不把旧波方程另称发现。
- 首阶Q在B0上计算；DQ[b]进入下一圈阶。740独立O(1)背景脉冲与本轮O(λ)反作用的两个计数轴分开。
- 仿射背景在紧条带、足够小参数下保时间片、正F和Lorentz型；共同过去输送保CAR正性、自对偶、Hadamard。
- 用同一半密度模式执行非选择记录；后态不假设Gaussian，不除以小概率。参考识别不是自主装置。
- 参数化来源估计沿已构造光滑路径，有限参数导数预算由高阶资料承担；未假定全未知来源无损失可微。
- 残差针对完整绝对方程和约束，仍非精确非线性解或λ=1精度证书。高圈、UV及全量子Gauss开放。
- 两种一阶辅助扩展保零阶准备，首阶背景不变但二阶实际历史项不同；零阶准备变化会改变首阶结果。
- 有限矩阵只校准完整原物种与态族，未计算连续重整化应力或联立PDE数值。
- 下一步先查现有局域操作及来源接口，复用524、575、633—634和716；不重复通用测量注能。
- 所有入口及冻结文件保留，目标与日程不变。
''')
main=('research_note_741.md','joint_absolute_source_development.py','joint_absolute_source_development_results.json','unified_physics_condition_ledger_741.md')
write('round741_drafts/final_review.txt','Primary-agent review only. Conditional absolute first-order backreaction and compatible state family; no exact nonlinear or physical lambda=1 conclusion.\n'+
    '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'741':'742','740':'741','739':'740','3419':'3421','3417':'3419',
         '1532':'1535','3359':'3373','3345':'3359','386':'387','291':'292',
         'joint_causal_eft_branch':'joint_absolute_source_development','total_pole_entry':'absolute_source_order_entry'}
publication=remap((HERE/'publish_round740.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第741轮完成：** [实际绝对来源与一阶共同反作用]({p}research_note_741.md)原全部线性场方程、初始约束及变背景记录态族与同一绝对来源接通，完整方程／约束残差O(λ²)。固定零阶准备时首阶不依赖辅助一阶扩展，二阶仍保历史差。两组、十八式通过，最新741／3421，1535份编号科学文件、3373份保护证据。[核验]({p}research_round_741_checks.json)、[条件账]({p}unified_physics_condition_ledger_741.md)。精确解、实际误差与内部准备仍开放。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（741后，优先于下方历史安排）：** 接[742现有记录的局域过程与共同来源]({p}round742_drafts/STATUS.md)，先回查524、575、593、620、633—634及716，核同一操作和总来源。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('原总反馈的增长极点与一致有限阶因果响应','实际绝对来源、联合初值与一阶共同反作用')
publication=publication.replace('旧空间合同保持，有限阶残差不替代完整解误差','旧空间合同保持，一阶绝对发展不替代内部准备')
write('publish_round741.py',publication)
write('postcheck_round741.py',remap((HERE/'postcheck_round740.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round740.py').read_text('utf8'),mapping)
verification=verification.replace('original_total_poles_and_causal_pulse_residual_checked=True','original_complete_source_hierarchy_and_reference_lift_checked=True')
write('verify_round741.py',verification)
print('Prepared741: full first-order source/state bridge and next local-process scope.')
