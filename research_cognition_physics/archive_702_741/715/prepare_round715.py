"""Prepare715 common-domain evidence and the smooth matter interface."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_714.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：715实际热读取与共同能源域\n'+rest
ledger=ledger.replace(
    '接[713全账](unified_physics_condition_ledger_713.md)，回填[714报告](research_note_714.md)。[结果](joint_center_thermal_reference_results.json)、[核验](research_round_714_checks.json)。',
    '接[714全账](unified_physics_condition_ledger_714.md)，回填[715报告](research_note_715.md)。[结果](joint_sharp_record_domain_results.json)、[核验](research_round_715_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**715当前增量：** 原完整Gauss热态的规范边缘在固定图上有严格正有限上下界；指定环路tanh读口锐化时，实际热平均注能为增益的线性阶。记录／后态有迹范数极限，但极限离开原有限能源域，共形来源同步发散。固定图结论不替代共同空间极限；下一项回到633光滑物质模的实际对象映射。\n\n'+marker)
updates={
 'C03 事件记录':'715原锐化读口虽有正常后态／记录极限，却失去有限能源；不能只验迹范数',
 'C19 参考态':'715完整Gibbs规范边缘与Haar可作严格正上下比较；不将有符号费米路径权当概率',
 'C20 尺度映射':'715固定图大增益的真实平均代价为线性阶，跨图常数未知；双重极限须另验',
 'C21 主体存储':'715有限结果字母表与后态收敛不能自动给有限记录资源；无限锐化路线限定失败',
 'C22 来源反作用':'715同一读取共形平均来源为负二倍注能，随锐化发散；与准备背景导数分开'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 714共同热参考与原带荷差量')
ledger=ledger[:at]+'''## 715实际参考上的读口与共同能源域

- 光滑椭圆全H、Gauss谱及643的矩阵核支配给原规范边缘rβ的正有限上下界；常数只在固定图成立。
- 原加权环字符密度g在零点严格正；真实Gibbs注能被Haar平均夹住，锐化增益z使均值为Θ(z)。
- 两支Kraus强收敛带来正常记录及后态迹范数收敛，但零水平上的跳跃使极限有无限形式能源。
- 同一共形来源同步失去有限性，后态收敛不能替代能源／来源收敛。
- 没有从固定图锐化推出空间连续失败，也没有选择不同参考或删去带荷物质。

## 本轮合并与下一项

C03／C19／C20／C22必须共用状态、操作与能源域。关闭指定原热参考上的无限锐化捷径，不继续优化锐化函数。

接[716](round716_drafts/STATUS.md)：复用633已成立的光滑连续物质读取，核原有限图sterile模、参考排序、质量和实际来源的映射；保留自由分支及因果实现限制。旧空间与四分支不变。
'''
write('unified_physics_condition_ledger_715.md',ledger)
write('round715_drafts/research_note_715_draft.md',(HERE/'research_note_715.md').read_text('utf8'))
write('round716_drafts/STATUS.md','''# 第716轮入口：光滑物质模与同一原过程

接[715](../research_note_715.md)、[全账](../unified_physics_condition_ledger_715.md)。指定原环路锐化在真实热参考上失去有限能源域，虽然记录后态收敛。停止继续扫描锐化形状和温度。

1. 回查633的原sterile光滑模、混合质量、Hadamard后态与来源，及其瞬时分布式读取的因果限制；不重新证明这些旧自由场结果。
2. 在原598完整有限图上寻找同一模和偶仪器的实际对应。全部动态标量、原32CAR、Gauss、质量和跳跃保持；不把冻结自由背景当全相互作用过程。
3. 666—667参考输送只保证已有过程的表示一致；当前需核原模、原热态、能源域及物理传播的共同匹配，保留604／613／690等已发现的物种或动力字典限制。
4. 若用额外跳跃／连续处方，逐项声明，不把canonical自由Dirac动力假装从原任意协变跳跃唯一得到。有限传播装置仍须内部资源，但先做对象与共同条件连接。

旧空间382—386、425、522—524及699的限定反例保持。目标是跨分支统一，不是单个仪器成本优化。
''')
write('round715_drafts/literature_scope_audit.json',json.dumps(dict(
    source='https://arxiv.org/pdf/1801.01335',primary_book='https://link.springer.com/book/10.1007/978-3-319-68903-6',
    read='Theorem II.1 and local elliptic regularity in its proof, pp35-37 of author PDF.',
    scope='Smooth elliptic semibounded self-adjoint operator on the full configuration bundle.',
    not_assumed='Positive fermionic path weights, a smooth gauge quotient, or spatial-continuum uniform bounds.',
    inherited='603 thermal trace; 623 common domain; 643 matrix kernel domination.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round715_drafts/scope_and_dedup_review.md','''# 715范围审查

旧热迹、共同域、数据处理和715入口最坏范数不重计。本轮补实际原热态的规范边缘界、锐化的平均线性阶与后态能源域不相容。

严格正热密度经物理本征函数及光滑Gauss包证明，不从有符号路径权读取。矩阵核标量支配控制非紧标量积分。相同H含全部物质。

有限z平均能源发散不足以推出极限无限能源；稿件另以规则零水平的Kraus跳跃证明局部H1失败。

数值为原正常Haar Gauss态及原系数校准，不是完整Gibbs数值。高阶图正则、全图热量和跨图常数解析／未知范围均注明。633自由连续读口及其因果限制直接继承。
''')
names=('research_note_715.md','joint_sharp_record_domain.py',
       'joint_sharp_record_domain_results.json','unified_physics_condition_ledger_715.md')
write('round715_drafts/final_review.txt',
      'Main-agent review; no independent review. Physical heat density positivity, domination, '
      'Haar coarea normalization, fixed-graph quantifiers and sharp-state domain reviewed.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared715 ledger, draft, reviews and716 entry.')
