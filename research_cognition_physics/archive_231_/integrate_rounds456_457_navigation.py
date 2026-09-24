"""Guarded navigation update: two completed rounds and user functional recursion."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE,ROOT=HERE.parent,HERE.parent.parent
PATHS=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
       BASE/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md']
original={p:p.read_bytes() for p in PATHS}
texts={p:b.decode('utf-8').replace('\r\n','\n') for p,b in original.items()}
assert '**第456轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summaries=[
'**第456轮完成：** [持续转换与读者共同运行](%sresearch_note_456.md)完整分类与455转换对易的全部静态两体项，并证明这类项不能向外输出新逻辑信息；非对易七边持续正例给整窗读者概率差>1／4，同时严格证明旧目标码概率<3／5，不能继承原转换保证。6项检查，累计2169项；679份编号科学文件、714份保护证据。未把对易当公理，未排除其它共同实现；辅助、接触和自主交接仍有输入。',
'**第457轮完成：** [共同资源的内部浓缩与持续供给](%sresearch_note_457.md)独立从455基线给任意n对SU(2)不变资源的精确部门容量；四对首次出现指定p区间的额外守恒损失。三对质量0.7的资源经固定正交换，整窗目标singlet概率>0.8；整数证书和余料熵／关联账完整。6项检查，累计2175项；682份编号科学文件、717份保护证据。未实现资源自生、免费复原或供给到接入的自主交接。'
]
for p in PATHS[:5]:
    prefix='research_cognition_physics/archive_231_/' if p==ROOT/'README.md' else ('' if p==HERE/'README.md' else 'archive_231_/')
    old=next(line for line in texts[p].splitlines() if '**第455轮完成：**' in line)
    added='\n\n'.join(('> ' if p==ROOT/'README.md' else '')+s%prefix for s in summaries)
    recursion=('> ' if p==ROOT/'README.md' else '')+'**用户新增方向：递归功能闭合。** [原则候选与去重](%srecursive_functional_closure_proposal.md)：大主体须由子系统分工重新实现小主体的必要认知功能与反馈，并能继续加入更大整体。先检验功能和接口的自相似，不预设各层等维、同一数值Hamiltonian或几何分形维数；不默默改写F＋U＋C＋P。此为研究合同草案，不增加科学轮次。'%prefix
    texts[p]=texts[p].replace(old,old+'\n\n'+added+'\n\n'+recursion,1).replace('231—455轮','231—457轮')
next_item='**下一项：递归功能闭合中的最小实际功能链。** 按用户新补充，大主体的子系统分工应重现小主体必要功能，并使整体成为可再次组合的同类主体。复用451的完整闭环候选、221的过程闭合和456—457的接续／资源结果，优先把接收区别、保留任务记忆、内部更新和对外响应纳入同一持续过程；给各项实际承担者及后续接口，不仅命名模块。456已证明对易附加项不能输出新L、一个非对易正例会破坏旧转换保证；457已给混合资源的持续供给，但余料关联与供给—接入自动接续仍未闭合。后继只补真实缺失的接续或反馈，不以外部切换、一般门编译或新增矩阵分块充当主体形成。功能自相似作为新增候选记录，先不要求无限层、等维、同一数值生成元或特定几何维数。资源、控制和记录仍在整体内部，三维及完整GR均未完成。'
for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',HERE/'README.md'):
    old=next(line for line in texts[p].splitlines() if line.startswith('**下一项：在同一持续规则中接续未知主体转换与共同交互。**'))
    texts[p]=texts[p].replace(old,next_item+'\n\n**455后的持续接续任务（456—457已有两项局部结果）：** '+old.split('** ',1)[1],1)
p=BASE/'research_direction.md'
texts[p]=texts[p].replace('最新科学轮次与检查数为455／2163','最新科学轮次与检查数为457／2175')
p=BASE/'RESEARCH_STATE.md'
texts[p]=texts[p].replace('最新455轮及累计2163项见本文开头','最新457轮及累计2175项见本文开头')
p=BASE/'README.md'
texts[p]=texts[p].replace('当前复算与冻结入口：','当前复算与冻结入口：[456科学核验](archive_231_/verify_continuous_interface_join_round.py)、[457科学核验](archive_231_/verify_singlet_supply_capacity_round.py)、[457整合核验](archive_231_/verify_round457_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('第三阶段累计2163项；旧科学证据保持原字节','第三阶段累计2175项；旧科学证据保持原字节')
p=HERE/'README.md'
texts[p]=texts[p].replace('当前完成455轮','当前完成457轮').replace('当前455不作为预定终点','当前457不作为预定终点')
row=next(line for line in texts[p].splitlines() if line.startswith('| [455：'))
texts[p]=texts[p].replace(row,row+
'\n| [456：持续接入与读取的共同运行](research_note_456.md) | 对易族精确分类；非对易共同读取；旧转换保证失败的严格证书 | [代码](continuous_interface_join_audit.py)、[结果](continuous_interface_join_audit_results.json)、[核验](research_round_456_checks.json)；6项检查；科学基线455 |'+
'\n| [457：共同资源的持续内部供给](research_note_457.md) | SU(2)部门容量；固定交换浓缩；余料熵与整体交付误差 | [代码](singlet_supply_capacity_audit.py)、[结果](singlet_supply_capacity_audit_results.json)、[核验](research_round_457_checks.json)；6项检查；科学基线455，与456独立 |',1)
texts[p]=texts[p].replace('当前整合入口：','当前整合入口：[456整合核验](verify_round456_integration.py)、[457整合核验](verify_round457_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('当前2163项科学检查；711份保护证据，其中本阶段编号科学文件676份','当前2175项科学检查；717份保护证据，其中本阶段编号科学文件682份')
p=HERE/'spatial_premise_closure_audit.md'
assert '## 84.' not in texts[p]
texts[p]+='''
## 84. 第456轮：单独接口保证不能直接相加

2026-09-24，科学基线455。转换C为旧A三体分别与载体4的等强交换；任意完整原始静态两体附加K与C对易，当且仅当三条A—4权重相同，且每个外部r与A及4四条权重相同。三点Pauli手征项支撑正交给完整Frobenius平方和恒等式，N体系数3·2^(N−1)；七体21权重的整数Gram为192 BᵀB，秩11、核10。旧A代码上全部该类K把L与G及外部解耦，因此不能向外输出新的L信息。这是选定对易合同的边界，对易不是新认知公理，也不是保住455目标码的充分条件。

给明确非对易共同正例：C三边与431重标号读者四边同时全程开启，共七边等强正交换，初始辅助34和读者56各为纯singlet。共同SU(2)与Schur引理给全部t的压缩效果I_G⊗F_L(t)，读者结果及任意内部R只依赖ρ_LR；没有把整个读者态或主体旧L保持也一起宣称。完整七体输出等距，未知GLR信息全局保留。

100阶整数对易子与有理余项给t＝3／2时Y±输入读者差>1／3，整个|δ|≤1／24窗口差>1／4。对同一H及GL最大混态，四体目标码概率由另一严格证书给<3／5，而孤立455转换同刻>99／100；相对455理想输出的半diamond误差>2／5。这只否定具体直接并接保留旧保证，未排除其它非对易设计。初始两对纯资源、接触、读取及时间单位仍是输入，自主交接未完成。

6项检查、11公式通过；累计2169项，679份编号科学文件、714份保护证据。未重做432全部关系冻结、453固定旧码或一般模拟。三维与GR未完成。

## 85. 第457轮：持续供给的正例与余料义务

2026-09-24，独立科学基线455。n对τp初始资源、全部共同SU(2)不变酉、指定一对singlet目标下，按总自旋j分别取本部门最大本征值，容量为剩余2n−2体的j重数，给任意n精确最优。四对在p≥1／4时首次比无守恒全谱排序少7pq²(p−q)；p＝.7给.8722而非.9016。五对仅triplet最大混态，全谱支撑可进入singlet目标，不变酉却最多184／243；任意有限n、p<1还因最高自旋部门得成功率<1。256维显式不变酉达到性仅属于该较大控制类，不冒充固定两体实现。

同时给实际固定正交换供给：三对p＝.7，H＝2S01＋S23＋S45＋(S02＋S03＋S05＋S12＋S14＋S15)／5。整数定点14阶指数、12次平方及全程舍入上界证明f(93)>809／1000，整个|t−93|≤1／200上f>4／5，无纯辅助、测量、后选择或中间脉冲。具体图、强度、初始关系资源仍是输入，未达到一般不变酉最优的声明。

全六体保留未知信息，供给时旧主体GLR恒等不被触碰。目标改善时余料熵必须增加；指定实例余料增约.46081 bit，目标—余料互信息约.10971 bit。局部供给到455的通道误差1−f<1／5，保留全部内部余料时只保证与适当纯singlet乘积比较态距离≤sqrt(1−f)，二者不能混用。能量均值−1和方差6399／1250为所选模型守恒量，不是准备热功。供给H与455连续交接仍未构造。

父级已独立核部门容量证明、全维最优构造、整数误差传播、余料熵账和参考量词，并复算保存结果。6项检查、15公式通过；两轮新增12项，累计2175项；682份编号科学文件、717份保护证据。历史证据原字节保留，无图像检验。

用户同期补充大主体与小主体应功能自相似：小主体各认知功能由大主体内相应子系统分工承担，整体还须能再次组合。已另存recursive_functional_closure_proposal.md，回查SoCA L7、87—88、221、450—452及新结果；记为新增候选，未静默加强F＋U＋C＋P。功能重现须包含真实接续与反馈，不仅模块名称；暂不要求相同几何形状、完整态等维、相同数值Hamiltonian或无限层无成本扩张。下一项把最小接收—记忆—更新—响应链的实际分工放入同一持续规则。完整宏观主体、三维与GR仍未完成。
'''
for p in PATHS:
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
for p,s in texts.items():
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
    newline='\r\n' if b'\r\n' in original[p] else '\n'
    p.write_bytes((s.rstrip()+'\n').replace('\n',newline).encode('utf-8'))
    print(p.relative_to(ROOT))
