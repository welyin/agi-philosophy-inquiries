"""Guarded navigation update for rounds 454 and 455 and the Spark comparison."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE,ROOT=HERE.parent,HERE.parent.parent
PATHS=[ROOT/'README.md',BASE/'README.md',BASE/'research_direction.md',
       BASE/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md']
original={p:p.read_bytes() for p in PATHS}
texts={p:b.decode('utf-8').replace('\r\n','\n') for p,b in original.items()}
assert '**第454轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summaries=[
'**第454轮完成：** [扩大共同关系后的持续交换正例](%sresearch_note_454.md)在预备共同G-singlet的明确合同下，以六体七条等强正交换实现四维逻辑经一个共同模式的精确演化，端点给未知逻辑及任意参考的纠缠酉；实际局部读数在整个指定窗口的概率差大于1／2。6项检查，累计2157项；673份编号科学文件、708份保护证据。4→5扩大以预备G为前提，不覆盖453的任意未知G16维接口；共同资源、接触来源与读取权限仍为输入。',
'**第455轮完成：** [未知主体到四体接口的连续交换转换](%sresearch_note_455.md)独立从453基线出发：旧三体任意GL及参考，加一对内部singlet，只用三条持续等强交换精确把L送入434四体码，把未知G及全部关联迁移到第五体；旧逻辑全时不变。给完整通道的有限窗误差、指定闭合合同的最少两辅助与五体纯资源条件。6项检查，累计2163项；676份编号科学文件、711份保护证据。未实现辅助来源、永久交接或联合网络。'
]
for p in PATHS[:5]:
    prefix='research_cognition_physics/archive_231_/' if p==ROOT/'README.md' else ('' if p==HERE/'README.md' else 'archive_231_/')
    old=next(line for line in texts[p].splitlines() if '**第453轮完成：**' in line)
    added='\n\n'.join(('> ' if p==ROOT/'README.md' else '')+s%prefix for s in summaries)
    bridge=('> ' if p==ROOT/'README.md' else '')+'**用户Spark方向提示已记录：** [经济交流与认知接口对照]('+prefix+'spark_economy_cognition_bridge_review.md)将加入非吞并、双边互认、真实交付、局部摘要及内部资源账转为功能候选；区分零净账与零占用／风险、账本解锁与实际影响消失。此为来源审查，不计科学轮次、不推出复量子；Spark文件未改。'
    texts[p]=texts[p].replace(old,old+'\n\n'+added+'\n\n'+bridge,1).replace('231—453轮','231—455轮')
next_item='**下一项：在同一持续规则中接续未知主体转换与共同交互。** 455已经把未知GL及参考无损交付到434的同一四体接口，454已有预备共同G下的精确交互正例；它们尚不是一个连续运行的共同主体。优先核455三条转换作用不关闭时，结合既有434—436的实际任务，是否还能给原接口可读、对未知参考一致的联合功能，或由显式内部状态决定接触而完成交接；不按外部精确时序手动拼接两个Hamiltonian，也不把永久停机新增为认知必要条件。共同singlet及辅助不能默认为免费、用后自动复原或由低能谱自动准备；必须保留来源、暂时占用和后续任务的资源账。Spark的双边互认、加入非吞并、任务相关摘要作为设计参照，但积分、信用和账目不直接等同量子振幅或能量。避免返回时间扫描、一般门编译和仅矩阵分块。仍沿429自然交互主线；物理位置、三维及完整GR均未完成。'
for p in (BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',HERE/'README.md'):
    old=next(line for line in texts[p].splitlines() if line.startswith('**下一项：宏观主体的可用组合接口。**'))
    texts[p]=texts[p].replace(old,next_item+'\n\n**453后的接口任务（454—455已分别给共同交互与未知输入转换正例）：** '+old.split('** ',1)[1],1)
p=BASE/'research_direction.md'
texts[p]=texts[p].replace('最新科学轮次与检查数为453／2151','最新科学轮次与检查数为455／2163')
p=BASE/'RESEARCH_STATE.md'
texts[p]=texts[p].replace('最新453轮及累计2151项见本文开头','最新455轮及累计2163项见本文开头')
p=BASE/'README.md'
texts[p]=texts[p].replace('当前复算与冻结入口：','当前复算与冻结入口：[454科学核验](archive_231_/verify_shared_relation_revival_round.py)、[455科学核验](archive_231_/verify_relational_interface_conversion_round.py)、[455整合核验](archive_231_/verify_round455_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('第三阶段累计2151项；旧科学证据保持原字节','第三阶段累计2163项；旧科学证据保持原字节')
p=HERE/'README.md'
texts[p]=texts[p].replace('当前完成453轮','当前完成455轮').replace('当前453不作为预定终点','当前455不作为预定终点')
row=next(line for line in texts[p].splitlines() if line.startswith('| [453：'))
texts[p]=texts[p].replace(row,row+
'\n| [454：共同关系扩大后的持续交换](research_note_454.md) | 五维精确演化；端点逻辑纠缠；完整物理效果整窗；共同准备边界 | [代码](shared_relation_revival_audit.py)、[结果](shared_relation_revival_audit_results.json)、[核验](research_round_454_checks.json)；6项检查；科学基线453 |'+
'\n| [455：未知关系主体的接口转换](research_note_455.md) | 三边集体交换；未知GLR保留；四体码交付；完整通道误差与纯资源条件 | [代码](relational_interface_conversion_audit.py)、[结果](relational_interface_conversion_audit_results.json)、[核验](research_round_455_checks.json)；6项检查；科学基线453，与454独立 |',1)
texts[p]=texts[p].replace('当前整合入口：','当前整合入口：[454整合核验](verify_round454_integration.py)、[455整合核验](verify_round455_integration.py)；历史入口：',1)
texts[p]=texts[p].replace('当前2151项科学检查；705份保护证据，其中本阶段编号科学文件670份','当前2163项科学检查；711份保护证据，其中本阶段编号科学文件676份')
p=HERE/'spatial_premise_closure_audit.md'
assert '## 82.' not in texts[p]
texts[p]+='''
## 82. 第454轮：扩大共同关系的一份可用正例

2026-09-24，科学基线453。六qubit分成两个三体关系主体，先提供两G因子的共同singlet，未知两逻辑及内部参考任意。取两块各三条内部交换及跨块S36，七条非零权重均为正1且持续开启。旧四维逻辑与唯一两块各j＝3／2的共同singlet形成精确五维不变空间，最小性由H|11>的非零第五分量直接证明。整数未归一化码列的范数平方8、24、24、72，桥列范数平方576，全部生成元交织关系逐元素核对。

五维H的前三项为−1、1、1，末块为[[1／3,2sqrt(2)／3],[2sqrt(2)／3,17／3]]。T＝π／(2sqrt(2))时精确返回旧逻辑，门为diag(exp(iT),exp(−iT),exp(−iT),−exp(−3iT))，条件相位−1；|++>产生最大逻辑纠缠。算符恒等式覆盖全部未知逻辑与参考。过程中第五模式占用最高1／9，不能只算PHP。

完整64维上的B逻辑X效果给两输入中心读数约0.19715和0.80285。由于内部块和效果对易，单概率导数界1、差的导数界2；三角有理估计给|t−T|≤1／32上差严格大于10101／20000。全部未知参考通道误差≤4|δ|，旧码漏出≤8δ²／9；原始静态扰动范数η额外误差≤|t|η，η≤1／1024时整窗差仍大于1／2。实际读取权限与时序仍是输入，没有建造全层析器。

4→5以已准备的共同G-singlet为前提，不能视作无损扩展453的任意未知G16维码。G00属于总自旋1，六体封闭交换保其部门，不能凭此获得所需总自旋0准备；同一端点该输入最大旧码漏出约0.39149。正例同时改变准备合同和全程旧码合同，未证明任意独立主体已无条件融入。Heinz等及van Meter–Knill的已有交换编码工具明确归属，未声称新普适门定理。

6项检查、12公式，独立整数／解析／文献审核通过；累计2157项，673份编号科学文件、708份保护证据。三维、共同资源的自主来源与完整GR仍开放。

## 83. 第455轮：未知信息可以迁移到明确的内部载体

2026-09-24，与454独立，科学基线453。原三体主体G、L及内部R全部未知，可任意关联；提供独立两体singlet辅助3、4，采用H＝S04＋S14＋S24。旧代码上的精确生成元为I＋SWAP_G4，L全时不变；三边持续同强作用，无中途读取、后选择或普适编译。

t＝π／2时，前四物理体成为434同一四体singlet码，逐列无逻辑字典修正，第五物理体接收原G及全部GLR关联。完整输入／输出等距恒等式给未知信息保留。原三个逻辑效果全空间合法、与H对易、转换前后保持同一含义，不因此宣称已构造完整指针和触发器。

全过程目标码概率1−3cos²t／4；相对端点完整通道的半diamond误差精确为sqrt(3)|sinδ|／2，来自两等距相对重叠为同一标量而非样本拟合。|δ|≤1／8时码概率≥253／256，通道误差≤sqrt(3)／16。持续H在tπ又使码概率回到1／4，没有永久交接。

在闭合保留全部GLR且目标四体singlet的合同内，逻辑码2维之外至少还需2维载体，五原始qubit、即新增两体为最少。最小五体中独立SU(2)不变辅助τp的秩／谱给任意酉目标支撑概率至多max(p,(1−p)／3)；确定性须p1纯singlet。混辅助相对理想整个端点通道误差为1−p。这个纯度结论不排除更多内部熵汇或资源，不升级为普遍认知原则。G与参考Bell关联搬运时原辅助对可变成I／4，故没有免费归还一份singlet。H均值3／2、方差3／4的守恒只属于所选模型，不是准备资源最小热功证明。

用户本轮引入Spark，已读取其code/docs/product经济／共同体／解释权分册和母稿关键条款，另存不编号的spark_economy_cognition_bridge_review.md。加入非吞并、双边互认、理论路径与真实交付区别、有限局部摘要、内部资源账是有用功能候选。零净账不推出零占用或零信用风险，数字锁回滚不消除实际交付影响；经典协议结构本身也不推出复振幅和纠缠。Spark文件未修改，来源对照不增加科学检查数。

455经独立复算与完整解析审查，6项检查、13公式通过。两轮共新增12项，累计2163项；676份编号科学文件、711份保护证据，历史科学文件原字节保留。下一项检验同一持续规则下的转换与共同任务接续、辅助来源及复用；不能把两个指定Hamiltonian外部排程当成已得到自主主体，也不把永久停机新增为必要公理。三维与完整GR仍未完成。
'''
for p in PATHS:
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
for p,s in texts.items():
    assert p.read_bytes()==original[p],f'Concurrent edit: {p}'
    newline='\r\n' if b'\r\n' in original[p] else '\n'
    p.write_bytes((s.rstrip()+'\n').replace('\n',newline).encode('utf-8'))
    print(p.relative_to(ROOT))
