"""Guarded navigation update for completed record-response and recursive-interface rounds."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE,ROOT=HERE.parent,HERE.parent.parent
PATHS=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
       BASE/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md']
original={p:p.read_bytes() for p in PATHS}
texts={p:b.decode('utf-8').replace('\r\n','\n') for p,b in original.items()}
assert '**第458轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summaries=[
'**第458轮完成：** [内部记录参与持续响应](%sresearch_note_458.md)由两体记录与目标的普通交换给参考一致的响应通道，以及适用全部共同态的激活上界；接收、记录和响应六边同时运行，严格共同窗口内两项信号均可辨。目标单体不必可辨，实例区别在目标—内部参考关系中；响应亦反作用于接收。6项检查，累计2181项；685份编号科学文件、721份保护证据（新增1份456纯文字修订）。未完成完整认知反馈或递归主体。',
'**第459轮完成：** [持续过程与递归接口的自相似](%sresearch_note_459.md)将任意有限交换网络及实际读者效果精确提升到复合块，保留未知私有信息与全部参考关联；11体复算继承431信号及严格时间窗。若始终保持各层编码且只用均匀层级接口，旧私有记忆不能驱动上层响应；未准备部门不能省略，有限接触预算下时间随规模增长。6项检查，累计2187项；688份编号科学文件、724份保护证据。接口及过程继承不等于完整功能闭合。'
]
for p in PATHS[:5]:
    prefix='research_cognition_physics/archive_231_/' if p==ROOT/'README.md' else ('' if p==HERE/'README.md' else 'archive_231_/')
    old=next(line for line in texts[p].splitlines() if '**第457轮完成：**' in line)
    added='\n\n'.join(('> ' if p==ROOT/'README.md' else '')+s%prefix for s in summaries)
    texts[p]=texts[p].replace(old,old+'\n\n'+added,1).replace('231—457轮','231—459轮')
next_item='**下一项：递归主体的任务记忆怎样驱动对外响应。** 458已给同一持续交换中的内部记录—响应机制，459已给实际过程与读取跨尺度继承，但459的均匀层级接口会把旧私有记忆隔离在外部交流之外。优先检验允许访问子接口或暂时改变祖先编码时，原有任务记忆能否真实影响上层响应，并保留未知信息、参考关联和后续可组合接口；不得把全部祖先编码永久不变加成认知公理。已有一般提升定理覆盖的换规模实例、三叉计数和431信号重算不另算新进展。需要给同一持续规则、真实参与者与资源账，不用外部按时切换、抽象模块命名或免费准备代替功能闭合。初始部门、共同辅助及接触来源仍是建模输入；功能自相似仍为新增候选，三维及完整GR未完成。'
for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',HERE/'README.md'):
    old=next(line for line in texts[p].splitlines() if line.startswith('**下一项：递归功能闭合中的最小实际功能链。**'))
    texts[p]=texts[p].replace(old,next_item+'\n\n**457后的功能链任务（458—459已有两项局部结果）：** '+old.split('** ',1)[1],1)
p=BASE/'research_direction.md'
texts[p]=texts[p].replace('最新科学轮次与检查数为457／2175','最新科学轮次与检查数为459／2187')
p=BASE/'RESEARCH_STATE.md'
texts[p]=texts[p].replace('最新457轮及累计2175项见本文开头','最新459轮及累计2187项见本文开头')
p=BASE/'README.md'
texts[p]=texts[p].replace('当前复算与冻结入口：','当前复算与冻结入口：[458科学核验](archive_231_/verify_continuous_record_response_round.py)、[459科学核验](archive_231_/verify_recursive_exchange_interface_round.py)、[459整合核验](archive_231_/verify_round459_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('第三阶段累计2175项；旧科学证据保持原字节','第三阶段累计2187项；旧科学证据保持原字节')
p=HERE/'README.md'
texts[p]=texts[p].replace('当前完成457轮','当前完成459轮').replace('当前457不作为预定终点','当前459不作为预定终点')
row=next(line for line in texts[p].splitlines() if line.startswith('| [457：'))
texts[p]=texts[p].replace(row,row+
'\n| [458：内部记录参与持续响应](research_note_458.md) | 全参考响应通道；实际记录—响应链；共同窗口及反作用 | [代码](continuous_record_response_audit.py)、[结果](continuous_record_response_audit_results.json)、[核验](research_round_458_checks.json)；6项检查；科学基线457 |'+
'\n| [459：持续过程与递归接口](research_note_459.md) | 任意有限网络与实际效果提升；私有记忆隔离边界；部门及成本条件 | [代码](recursive_exchange_interface_audit.py)、[结果](recursive_exchange_interface_audit_results.json)、[核验](research_round_459_checks.json)；6项检查；科学基线457，与458独立 |',1)
texts[p]=texts[p].replace('当前整合入口：','当前整合入口：[458整合核验](verify_round458_integration.py)、[459整合核验](verify_round459_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('当前2175项科学检查；717份保护证据，其中本阶段编号科学文件682份','当前2187项科学检查；724份保护证据，其中本阶段编号科学文件688份')
texts[p]=texts[p].replace('另含423的两份草稿','另含[456公式排版修订](research_note_456_text_v2.md)一份（仅补漏写的公式间隔反斜杠，旧原稿、旧结果均保留）、423的两份草稿',1)
p=HERE/'spatial_premise_closure_audit.md'
assert '## 86.' not in texts[p]
texts[p]+='''
## 86. 第458轮：记录参与后续响应的持续机制

2026-09-24，科学基线457。以已有431的物理记录34为M，5为未知目标，6为内部参考；持续响应K＝S35＋S45满足KPsinglet＝Psinglet。独立初始记录τp下，任意目标及参考的精确通道为退极化，系数1−32(1−p)sin²(3t／2)／27。该独立输入公式不能以共同动态中瞬时记录概率代入，替代存在关联的实际过程。

实际接收与响应由H＝S01＋S12＋S23＋S34＋S35＋S45全程同时运行。目标—参考singlet响应电流完全支撑在记录triplet部门，整数恒等式给范数sqrt(2)，对全部共同态有响应速度不超过sqrt(2)乘triplet概率；这给必要激活条件，非triplet必然响应。初始记录singlet，因此变化须由内部接收过程引入。未知源GL及全部R由整体等距保持，实际两项读取效果经Schur引理只依赖L及R。

100阶整数对易子与有理余项证明t＝3／2时记录、响应两项Y输入差都在1／4与3／10之间；|δ|≤1／64整窗内分别>7／32与>3／16。共同换轴不变实例中目标单体始终最大混态，可辨信息在目标与内部参考关系；不能宣称目标单体输出不同经典结果。响应K虽与记录效果对易，嵌套对易子平方范数96非零，实际接收差比原431严格减小>11／30。两对纯辅助、固定接触和读取仍为输入，没有完整学习、目标选择、无限存储或反馈闭合证明。

6项检查、13公式通过；累计2181项，685份编号科学文件、721份保护证据。增加的第4份保护证据仅为456文字v2，补公式间隔的漏反斜杠，科学原稿与旧结果原字节保留。父级与独立审查均核通道、激活界、整数证书及未知参考范围，无图像检验。三维及GR未完成。

## 87. 第459轮：过程继承与私有记忆参与是两项要求

2026-09-24，独立科学基线457。奇数原始单元的选定总自旋1／2部门分解为外露G与私有重数M。任意两块均匀全交叉交换等于标量加G间SWAP，私有M恒等；三体对三体情形复用453，不另宣称发现。由生成元等距交织得到任意有限网络、同时非对易接触、任意私有信息与参考关联的全时间过程继承。复合读者实际总自旋零投影也交织成两G的singlet效果；效应数学构造未充作探测装置或制备机制。

431实际五参与者被提升为3、3、3、1、1的11体、22条固定正交换，完整256列嵌入检验无泄漏，原始态演化复算两概率.169790866874与.866342012030。中心>2／3和整个旧窗口>1／2直接继承431严格证书；没有把同一信号重算计作新发现。全部五块都取三体时由定理覆盖15体36边，未假装全维数值求解。

递归三子块的G再选1／2部门，父M含新增关系L及全部子M；9体二层代码保留16维私有信息。但若固定层级只允许每个节点的均匀子块交换，且永久保持全部祖先代码，节点私有L动力学分解，各L不能驱动根G外部响应。该限制不适用于任意内部不变Hamiltonian，也未新增为认知公理。实际源关系读取可以访问其子G，因此过程提升并未同时证明所有祖先代码不变。

三体未准备部门的明确四体反例给相同接触在1／2部门实现交换，而3／2部门为整体相位；未知部门不能当成同一qubit端口。额外原始单元系数预算κ下，本均匀读者模板的时间至少3n／κ、接触数4n²；这是指定模型成本，不是热功或普适通信下界。三叉递归不推出几何分形维数或三维空间。

6项检查、15公式通过；两轮新增12项，累计2187项，688份编号科学文件、724份保护证据。父级独立审核证明和代码并复算保存结果；历史科学证据保持原字节。下一项不是更大块扫描，而是旧任务记忆如何在同一持续规则中影响上层实际响应，并保留后续可组合接口。完整递归认知、三维与GR仍未完成。
'''
for p in PATHS:
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
for p,s in texts.items():
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
    newline='\r\n' if b'\r\n' in original[p] else '\n'
    p.write_bytes((s.rstrip()+'\n').replace('\n',newline).encode('utf-8'))
    print(p.relative_to(ROOT))
