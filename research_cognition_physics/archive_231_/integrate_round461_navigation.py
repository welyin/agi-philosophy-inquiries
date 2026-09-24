"""Guarded navigation update after independent admission audit."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE,ROOT=HERE.parent,HERE.parent.parent
paths=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
       BASE/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md']
raw={p:p.read_bytes() for p in paths}
texts={p:v.decode('utf-8').replace('\r\n','\n') for p,v in raw.items()}
assert '**第461轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
for p in paths[:5]:
    prefix='research_cognition_physics/archive_231_/' if p==ROOT/'README.md' else ('' if p==HERE/'README.md' else 'archive_231_/')
    quote='> ' if p==ROOT/'README.md' else ''
    old=next(line for line in texts[p].splitlines() if '**第460轮完成：**' in line)
    new=quote+'**第461轮完成：** [独立未知主体加入的单一端口边界](%sresearch_note_461.md)对任意固定有限独立辅助，完整新主体总j＝1／2目标的最佳最坏权重恰为1／3，至少一个输入距任意目标支持态≥2／3；不是成功率1／3的无损接纳协议。保留两个j部门的84维相干接口可持续保存全部未知信息；删除部门相干会改变后继伙伴读数。6项检查，累计2199项；694份编号科学文件、731份保护证据（含一份461文字审阅前稿）。单一表示合同未升为认知公理，三维及GR仍开放。'%prefix
    texts[p]=texts[p].replace(old,old+'\n\n'+new,1).replace('231—460轮','231—461轮')
next_item='**下一项：丰富共同接口中的实际子接口与参与关系来源。** 460给精确保旧接口时的最小额外作用权限，461说明独立未知加入不能普遍强迫为单一G-qubit类型；后继应允许共同主体的多个表示、私有记忆和相干完整保留。优先核在旧处理持续运行时，哪些子接口参与交流能由内部实际记录决定，并给未知新主体及参考的完整后继保证；复用433—449的动态参与／内部载体及434—436的成熟模拟工具，只补尚缺的自然参与来源和真实功能接续。固定一张新接触表、再求更大不变子空间、重做一般高阶项编译均不算关闭来源缺口。不得把每层相同不可约表示、永久保持全部旧编码或免费恢复资源新增为认知必要性。用户的粒子猜想保留为同一动力学下稳定模式的后续解释检验，不把内部j或组织层数对应成物理粒子或空间维数；当前先落实空间的阶段目标保留。'
for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',HERE/'README.md'):
    old=next(line for line in texts[p].splitlines() if line.startswith('**下一项：递归主体的任务记忆怎样驱动对外响应。**'))
    texts[p]=texts[p].replace(old,next_item+'\n\n**459后的任务记忆响应问题（460—461已区分作用权限与接口类型）：** '+old.split('** ',1)[1],1)
p=BASE/'research_direction.md'
texts[p]=texts[p].replace('最新科学轮次与检查数为460／2193','最新科学轮次与检查数为461／2199')
p=BASE/'RESEARCH_STATE.md'
texts[p]=texts[p].replace('最新460轮及累计2193项见本文开头','最新461轮及累计2199项见本文开头')
p=BASE/'README.md'
texts[p]=texts[p].replace('当前复算与冻结入口：','当前复算与冻结入口：[461科学核验](archive_231_/verify_independent_subject_admission_round.py)、[461整合核验](archive_231_/verify_round461_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('第三阶段累计2193项；旧科学证据保持原字节','第三阶段累计2199项；旧科学证据保持原字节')
p=HERE/'README.md'
texts[p]=texts[p].replace('当前完成460轮','当前完成461轮').replace('当前460不作为预定终点','当前461不作为预定终点')
row=next(line for line in texts[p].splitlines() if line.startswith('| [460：'))
texts[p]=texts[p].replace(row,row+'\n| [461：未知主体的单一端口接纳边界](research_note_461.md) | 任意辅助最佳最坏权重1／3；参考扰动；更丰富相干接口 | [代码](independent_subject_admission_audit.py)、[结果](independent_subject_admission_audit_results.json)、[核验](research_round_461_checks.json)；6项检查；科学基线459，与460独立 |',1)
texts[p]=texts[p].replace('当前整合入口：','当前整合入口：[461整合核验](verify_round461_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('当前2193项科学检查；727份保护证据，其中本阶段编号科学文件691份','当前2199项科学检查；731份保护证据，其中本阶段编号科学文件694份')
texts[p]=texts[p].replace('另含[456公式排版修订]','另含[461文字审阅前稿](research_note_461_draft_before_text_review.md)一份（修公式及明确代码偏迹，未改科学结果）、[456公式排版修订]',1)
p=HERE/'spatial_premise_closure_audit.md'
assert '## 89.' not in texts[p]
texts[p]+='''
## 89. 第461轮：未知加入不能默认单一表示端口闭合

2026-09-24，独立基线459。两旧主体各有未知G⊗M，任意固定有限辅助τ与它们独立；辅助全部计入新共同主体，目标仅总j＝1／2。任意共同SU(2)不变酉守恒完整新主体的目标投影，因此时间、接触或旧码中间泄漏均不能提高目标权重。辅助初始部门权重只可能通过j＝1／2、3／2贡献；完整效果用前者一阶矩、后者二阶矩给出，任意跨j与重数相干已覆盖。

六份允许的独立平行G输入沿±X、±Y、±Z，其平均恰为triplet最大混态，目标权重平均q1/2／3＋q3/2／6≤1／3。故至少一个输入≤1／3。辅助I／2给F＝Ps＋Pt／3达到最佳最坏值；纯三体辅助将G与自身私有L纠缠亦达到，纯化载体全在辅助内部。该达到只属目标部门权重，未构造无损后选择接纳。

任意目标支持态与至少一个实际输出迹距离≥2／3，相应完整通道半diamond距离亦≥2／3。指定Bell参考输入后选择成功概率1／2，参考从I／2变diag(1／3,2／3)，扰动1／6；同一封闭不变规则不能让失败部门重新进入目标。单一表示、独立辅助及无余料外放都是额外合同，未推成所有未知主体不能组合。

两三体旧主体加单辅助可完整留在总j＝1／2、3／2的相干直和，重数各14，总维数84，包含原32维输入。全部固定两体交换保持它；原旧处理与新增接触全程同时运行，完整84维交织保留未知G、L及R。已有431私有关系读取作为该规则家族子例直接复用，不重算同一信号。

不能进一步只保留经典j标签：原允许独立A、B及I／2辅助输入，加入一个独立但显式同取向准备的新伙伴，实际singlet效果在是否删除旧总j相干时为0及1／18。整数分子2、分母36，且核原三体代码中的真实嵌入。此准备关系未宣称免费，内部j未解释为物理自旋。

6项检查、14公式通过，父级独立读证明与全代码并复算结果。累计2199项，694份编号科学文件、731份保护证据；额外保护461文字审阅前稿，主笔记仅修损坏公式及补代码偏迹解释，旧字节哈希完整保留。下一项回到丰富共同接口中的实际子接口和参与关系来源，不继续用固定单一G端口当普遍类型，也不把大不变子空间本身当完整认知。三维、粒子识别及GR仍未完成。
'''
for p in paths:
    assert p.read_bytes()==raw[p],f'Concurrent edit: {p}'
for p,s in texts.items():
    assert p.read_bytes()==raw[p],f'Concurrent edit: {p}'
    newline='\r\n' if b'\r\n' in raw[p] else '\n'
    p.write_bytes((s.rstrip()+'\n').replace('\n',newline).encode('utf-8'))
