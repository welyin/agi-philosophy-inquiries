"""Preserve old science, create the 708 ledger, review and continuation entry."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def write(name,content):
    path=HERE/name
    path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:
        f.write(content)


ledger=(HERE/'unified_physics_condition_ledger_707.md').read_text('utf8')
head,rest=ledger.split('\n',1)
ledger='# 联合条件总账：708继续认知的共同涨落资料\n'+rest
ledger=ledger.replace('接[706全账](unified_physics_condition_ledger_706.md)，回填[707报告](research_note_707.md)。[结果](joint_geodesic_spatial_block_results.json)、[核验](research_round_707_checks.json)。',
    '接[707全账](unified_physics_condition_ledger_707.md)，回填[708报告](research_note_708.md)。[结果](joint_block_fluctuation_contract_results.json)、[核验](research_round_708_checks.json)。')
marker='## 当前共同对象及仍存在的分支'
assert marker in ledger
ledger=ledger.replace(marker,'**708当前增量：** 在原固定图／等正权／完整CAR模型内，继续认知的指定任务要求把集体读口预算、原动能主系数与Majorana差分耦合绑定到同一singlet二次量P₂。仅保嵌套中点与幅度不够；两矩也不自动闭合所有新读口。候选认知原则产生具体的资料保留约束，没有选定原群、曲率、Y或物理维数；完整粗量子态与动态共同实现仍开放。\n\n'+marker)
updates={
    'C02 区域组合':'708静态嵌套与保后续操作分开验收；同一中点／幅度纤维内读口预算仍可不同',
    'C03 事件记录':'708声明Z读口的概率仅依赖Z，原H注能却需P₂；实际自主装置未构造',
    'C04 内部演化':'708原全H的singlet双交换子确定粗动能资料，禁止仅验线性漂移即签收完整闭合',
    'C17 质量机制':'708全部原Majorana均匀—差分平方耦合为||Ds||²(P₂/n−Z²/n²)，与同一读口代价关联',
    'C18 参数':'708在指定任务内减少独立拟合预算与差分耦合的自由度；原Y、R及节点权仍为输入',
    'C20 尺度映射':'708保(n,A,m)不足以确定所列动力／资源资料；加入P₂可算本轮瞬时量，但新P₂读口再需P₄',
    'C21 主体存储':'708原势保证新的Z读口保有限能量形式域；代价全域无统一常数界，实际内部实现仍缺'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            c=row.split('|');c[2]+='；'+value;rows[i]='|'.join(c);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 708候选认知要求及其联合验收

用户确认少量认知假说可以作为共同选择条件。本轮采用的新增要求是：在原紧支撑光滑Gauss核及延伸的有限能量态上，指定集体singlet读口不仅保结果函数，也保其完整H能量代价及原质量差分资料。任务量词本身是明示输入，没有从FUCP推出。

- z=s/√F为原规范singlet嵌入坐标。原H⁵给|∇z|²=1+z²/R²；其余完整H项与z乘法对易。双交换子将原动能与仪器预算接通，无须去掉规范电项或费米子。
- 新读口L±=√(1/2±sinZ/4)的注能为hbar² λ(Z)(n+P₂/R²)/(2w)，非负；598原势控制其有限能量形式域，但不是全域有界注能仪器。
- 原全部节点质量的均匀—差分平方强度由目标坐标协方差与原Yukawa Gram矩阵收缩得到；Majorana部分只需P₂和Z。预算与该差分强度不能独立指定。
- 四点保持相同完整聚合向量、幅度及平均质量，读口预算与Majorana差分仍不同。说明当前配置资料不足，不否定保隐含记忆的量子通道或整个认知计划。
- P₂并非全部未来任务的终点：原双交换子给P_k与P_l的P_(k+l)项。八点见证保S与P₂而P₄不同。固定n可用Newton身份，不能宣称所有有限状态描述均失败；停止继续逐矩制造同类反例。
- 独立输入减少的准确范围是当前模型的指定读口预算／粗动能／差分质量共同系数；原群、图、曲率、Y、实际装置和原宇宙参考仍未选定。

## 本轮合并与下一项

C02／C03／C04／C17／C18／C20／C21在同一原模型中接受同一组任务量词的约束。原四分支及699对指定辅助候选的严格反例范围全部保留；原连续手征、量子几何及现实拟合仍开放。

接[709](round709_drafts/STATUS.md)：比较保条件性内部态／记忆与明确近似任务的可行方案，使原参考、嵌套资料、实际读口预算及质量响应共同输送。先给同一模型内可检验的条件，不把一般Petz恢复或经典条件期望直接当这里的实现。旧空间382—386、425、522—524照原条件复用，不恢复已删除假设。
'''
write('unified_physics_condition_ledger_708.md',ledger)
write('round708_drafts/research_note_708_draft.md',(HERE/'research_note_708.md').read_text('utf8'))
write('round709_drafts/STATUS.md','''# 第709轮入口：共同任务下的内部态与有效记忆

接[708](../research_note_708.md)及[全账](../unified_physics_condition_ledger_708.md)。708完成指定读口预算、原动能和Majorana差分耦合的共同P₂资料，反驳仅保中点／幅度即可支撑这些任务；继续逐矩追加不能代替完整粗过程。

用户已确认以少量认知假说协调缺口，具体工程设计仍后置。下一项优先检验可行连接，至少区分：

1. 完整微观资料保留的平凡基线；只保静态(n,A,m,P₂)的有限任务描述；带内部条件态／记忆的描述。任务族、态范围、精度与时间必须先写明，不能为了通过测试更换原H或参考。
2. 原全Gauss Gibbs态和真实有限记录既有结果可复用；热参考上的条件平均不自动给其它准备或受扰后的实际内部态。检查质量响应与读口预算是否共同接受同一个条件态，不能分别重拟合。
3. 对接成熟充分性／条件期望／投影动力学时，核量子态、子代数、定义域、参考与物理Gauss支持的实际映射。经典条件期望可作为指定任务的比较基线，不冒充完整量子通道。
4. 如只有旧式的直接推论，记入口不计整轮；不要重复一般偏迹、Gaussian记忆、固定图精度或继续逐矩反例。优先找可同时保若干部门且减少独立输入的构造。

382—386、425、522—524直接复用；384的额外Lipschitz不恢复，386／425替代，旧端点映射仍待与原共同模型接上。内部H⁵不当物理维数。四分支及699范围保持，统一目标未完成。
''')
literature=dict(date='2026-10-03',sources=[
    dict(url='https://arxiv.org/abs/math-ph/0412093',read_scope='primary abstract',
         supported='sufficiency relative to a specified state family; not automatically full future operations',
         not_claimed='no recovery theorem applied or recovered quantum channel constructed'),
    dict(url='https://arxiv.org/abs/0906.4865',read_scope='primary abstract',
         supported='equilibrium coarse information distinguished from conditional-expectation effective dynamics',
         not_claimed='no overdamped Langevin error theorem transferred to original quantum H')])
write('round708_drafts/literature_scope_audit.json',json.dumps(literature,ensure_ascii=False,indent=2)+'\n')
write('round708_drafts/scope_and_dedup_review.md','''# 708去重与主代理范围审查

上次目标轮次归类为进展：707及708入口已有核验科学证据。本次之前的用户确认改变已授权方法并已写入导航，本身不计科学轮次。开始时Get-Process未发现Python执行；未覆盖任何既有草稿或冻结结果。

回查574、598、636、640—642、707及708入口；搜索已有编号报告中的矩闭合、矩层级、carré、质量方差和主符号粗化词。此前已有一般记忆、量子充分性与有限矩不等于完整态，本轮不重报这些为新基础定理。实质增量是原全H的singlet预算、原全部Majorana差分系数及同一P₂的实际连接。

主代理检查：交换子符号；Kraus导数因子；原势对新读口形式域的控制；实正交CAR变换中Majorana的转置；四点纤维的规范区别和次浸没；八点见证与有限n Newton身份；D_s=0退化；等权输入与变权范围；配置系数不足不外推为任意量子粗化失败。数值不代替证明，没有独立代理审查、全Gibbs谱或图像检查。

已检查三组实际数值结果。常数R、Y、w及具体读口仍为输入。预算与质量的一致性是本模型内的共同限制，不宣称已由认知原则唯一选出物理。
''')
names=('research_note_708.md','joint_block_fluctuation_contract.py','joint_block_fluctuation_contract_results.json',
       'unified_physics_condition_ledger_708.md')
review='708主代理最终审查：上述范围检查已完成；没有独立代理审核。\n'
for name in names:review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round708_drafts/final_review.txt',review)
print('Prepared708 ledger, preserved draft, source/scope audit and709 continuation.')
