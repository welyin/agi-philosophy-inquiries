"""Prepare 714 evidence and the next interface without changing frozen files."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_713.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：714中心信息、带荷关联与共同热参考\n'+rest
ledger=ledger.replace(
    '接[712全账](unified_physics_condition_ledger_712.md)，回填[713报告](research_note_713.md)。[结果](joint_holonomy_readout_scale_results.json)、[核验](research_round_713_checks.json)。',
    '接[713全账](unified_physics_condition_ledger_713.md)，回填[714报告](research_note_714.md)。[结果](joint_center_thermal_reference_results.json)、[核验](research_round_714_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**714当前增量：** 原整面中心变换的完整H差量仅含带荷梯度及双重态跨边跳跃；原Gibbs参考上的平均差量严格正。原环路读数、相对能源和几何响应共用同一热关联；来源差与代价的几何导数还相差参考响应项。没有计算原热环路均值或建立细化一致界，也未构造自主控制。\n\n'+marker)
updates={
 'C03 事件记录':'714原二元环路记录的对比给同一原Gibbs中心变换相对能源下界；尚未知原热均值是否非零',
 'C14 规范群':'714非平凡弱中心使原Q/L双重态切面跳跃变号；电与全部磁面项保持',
 'C19 参考态':'714原完整Gibbs上的中心差量严格正；全部原物质保持，同一带荷关联决定信号与相对代价',
 'C20 尺度映射':'714固定图严格正性不提供细化一致上下界；原热环路观测及共同尺度映射仍须匹配',
 'C22 来源反作用':'714中心前后平均来源差与代价的背景导数相差原热参考响应；不能只微分显式系数'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
at=ledger.index('## 713原环路、instrument和承载态的同一资源账')
ledger=ledger[:at]+'''## 714共同热参考与原带荷差量

- 中心变换保电动能、完整磁势和节点质量，改变原Log梯度及Q/L跨边跳跃；代价未另配一份H。
- 同一原Gibbs态与中心变换态的双向相对熵均等于βW_c；物理Gauss正常波包见证保证每个固定图上W_c>0。
- 原二元instrument要求βW_c≥2q artanh(q)，q=ηmβ。尚未求出mβ，不将713非热态的均值换成热均值。
- 用H与CHC的辅助插值，mβ与W_c分别为同一V的热交叉关联与自关联积分。中间H_s不是新物理公理。
- 几何平均来源差为Trρ∂ΔH；代价导数另有−βCov(G,ΔH)，参考响应不能省略。

## 本轮合并与下一项

C03／C14／C19／C20／C22共享原带荷奇部分V和同一热参考。没有自主控制协议、完整热均值数值或空间连续代价界；原四分支与旧空间范围保留。

接[715](round715_drafts/STATUS.md)：回查实际环路观测及原参考在共同物理尺度中的输送，区分固定图热差量、复合观测重整化与真正内部制备。停止重复一般被动性或温度扫描。
'''
write('unified_physics_condition_ledger_714.md',ledger)
write('round714_drafts/research_note_714_draft.md',(HERE/'research_note_714.md').read_text('utf8'))
write('round715_drafts/STATUS.md','''# 第715轮入口：原热环路信息与共同尺度映射

接[714](../research_note_714.md)、[全账](../unified_physics_condition_ledger_714.md)。原Gibbs中心变换有严格正平均差量，同一带荷关联约束原环路记录、相对能源与几何响应。固定图严格正不能替代细化一致界。

1. 先回查638—646、663—667、700—701及707—714，辨认原F_M与共同尺度观测／实际参考的对象映射，不再重证热迹、被动性或一般数据处理。
2. 当前原热均值mβ是否非零尚未证明；714的相对能源下界是条件性的，不能借用713非热包的数值。若原始字符在连续量子态中消失，需要说明重整化观测是否对应合法有限增益的实际读口。
3. 原完整带荷部门、正H、Gauss、参考来源必须共同保留；不得仅以无物质纯电模型、独立旋转器或任意小矩阵结果宣布原模型解决。
4. 优先给能排除错误连接、压缩共同条件的证明；若只是成熟公式推论，记入口，不计整轮。共同物理时间、真实操作和699的范围继续区分。

旧空间382—386、425、522—524及四分支保持。没有将当前工作缩为能源优化；最终仍需同一量子过程、时空、物质、引力与记录的共同连接。
''')
write('round714_drafts/literature_scope_audit.json',json.dumps(dict(
    primary_source='https://arxiv.org/pdf/1005.4495',read_equations='(2)-(5)',
    use='Established Gibbs unitary work/relative-entropy identity, not a new theorem.',
    not_imported='External driving does not supply an internal autonomous controller.',
    inherited='638 data processing; 603/623/703 common thermal domains and response.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round714_drafts/scope_and_dedup_review.md','''# 714范围与去重

前轮713及714入口有实质进展。本轮复用589一般梯度、598全部CAR、603热迹、623共同域、638相对熵与703热响应；不重计这些旧定理。

新增为原中心差量的具体带荷表达及物理Gauss非零见证，并将热环路读数、相对能源、原几何来源约束到同一V。非热入口等能状态不替换原Gibbs。

完整H热结论为解析；数值只计算原系数、来源、物理见证和条件读数函数。未算原热均值，没有小矩阵代替全H。无自主装置、细化统一能源界或动态GR结论。
''')
main=('research_note_714.md','joint_center_thermal_reference.py',
      'joint_center_thermal_reference_results.json','unified_physics_condition_ledger_714.md')
write('round714_drafts/final_review.txt',
      'Main-agent review; no independent reviewer. Domain, full cross gradients, Gauss witness, '
      'thermal-reference derivative and numerical scope reviewed. Sixteen equations; three groups.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
print('Prepared714 ledger, preserved draft, scope review and715 entry.')
