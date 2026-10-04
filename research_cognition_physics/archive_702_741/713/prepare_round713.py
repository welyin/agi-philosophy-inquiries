"""Prepare713 ledger/review; retain the initial diagnostic and entry."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as out:out.write(text)


ledger=(HERE/'unified_physics_condition_ledger_712.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：713共同尺度中的环路读口与状态电能\n'+rest
ledger=ledger.replace('接[711全账](unified_physics_condition_ledger_711.md)，回填[712报告](research_note_712.md)。[结果](joint_vertex_shared_evolution_results.json)、[核验](research_round_712_checks.json)。',
    '接[712全账](unified_physics_condition_ledger_712.md)，回填[713报告](research_note_713.md)。[结果](joint_holonomy_readout_scale_results.json)、[核验](research_round_713_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**713当前增量：** 原G上的平坦全局环可在单链趋单位时保整环对比；指定薄环instrument的完整H最坏注能随a^-2增长，固定横截面平均后有一致界。但任意正常态的非零平均原字符还要求裸电能至少随a^-4增长，明确正常Gauss态验证此限制。原几何来源同步，裸电能不当参考减除后的准备功，699反射型尚未沿该族输送。\n\n'+marker)
updates={
 'C02 区域组合':'713保全原图的平行闭环平均；不是617的表示切边，也未得到独立区域操作',
 'C03 事件记录':'713指定Gauss平方根仪器的原H注能精确；薄环最坏预算与横向平均预算分开',
 'C14 规范群':'713原L=(1,2)_{-3}字符下降到Z6商，弱中心holonomy在平坦细化中保留',
 'C19 参考态':'713任意正常有限电能态满足信号下界；显式正常Gauss态非Gibbs，共同参考减除仍开放',
 'C20 尺度映射':'713真实给定三维细化中，读取注能有界不自动保非零原字符信号的参考资源或量子连续极限',
 'C21 主体存储':'713承载态电能与读口预算须共同核算；集体instrument自主实现及记录总成本未构造',
 'C22 来源反作用':'713同一电主系数固定注能及共形-2、剪切来源；不任意减掉背景能源或其几何导数'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 713原环路、instrument和承载态的同一资源账

- 原L字符而非独立SU2模型；逐边Z6、局部Gauss、完整32模式相容。全局平坦holonomy不是瞬子数。
- 原完整H的所有乘法势及CAR保持；电双交换子决定指定集体仪器的注能。原非对角电项没有从H删除。
- 最坏预算精确为hbar² k_x N eta² C*/(4M)。固定物理横截面时它有限，但精确中心配置的实际注能为零，不能把最坏值当所有态的下界。
- 任意正常态信号m要求电能至少hbar² k_x NM C_R² m²/(4C*)。原Haar倾斜Gauss态给真实有限图信号，实际电能满足下界。
- 未减参考的电能不等于准备净功；波函数族每个有限图合法不等于已有共同连续态。参考、总能源和几何来源仍须共用。

## 本轮合并与下一项

C02／C03／C14／C19／C20／C22现由同一群、几何、原电项和状态信号限制。没有签收699负方向的连续存活或消失，也没有推出三维、面积熵或GR。

接[714](round714_drafts/STATUS.md)：回到原完整参考，核裸电能中哪些已属参考背景、哪些属实际状态差异，并保持几何变分和总H身份。旧固定图热态与熵不重证；不以随意真空减常数或局部顶点修补替代共同过程。
'''
write('unified_physics_condition_ledger_713.md',ledger)
write('round713_drafts/research_note_713_draft.md',(HERE/'research_note_713.md').read_text('utf8'))
write('round714_drafts/STATUS.md','''# 第714轮入口：共同参考、相对资源与几何来源

接[713](../research_note_713.md)、[全账](../unified_physics_condition_ledger_713.md)。原群平行环instrument在固定横截面下有一致注能界，但承载非零原字符均值的任何正常态有裸电能下界，显式原Gauss态也满足。这里不能把裸电能当准备功，亦不能无账减去它。

1. 回查603、623、637、666—667及702—706原完整Gibbs、参考运输、来源和状态收敛。停止重复固定图热迹／熵定理。
2. 明确同一个参考与测量任务，核E_E的截止增长哪些为参考已有背景，哪些为保持信号所需增量。原H中其它项、CAR和Gauss都保留；不能改为无物质纯电模型后宣称原系统已解。
3. 如使用参考能源或真空减除，必须同步处理其几何依赖与来源；这不是任意常数调零，也不能自动获得GR的宇宙学常数选择。
4. 原始Wilson字符的量子连续极限／重整化尚未构造。713的固定宽度仪器只控制一个实际读口预算；不要据此把699负方向归为UV伪影或宣布它连续保留。700—701的同一时间与来源量词继续适用。

先寻找能同时减少参考、资源与来源独立选择的真实合同；若只是旧公式推论，记入口，不计整轮。旧空间382—386、425、522—524和四分支范围保持，统一目标不变。
''')
write('round713_drafts/literature_scope_audit.json',json.dumps(dict(date='2026-10-03',sources=[
    dict(url='https://arxiv.org/abs/1006.4518',read_scope='abstract: Wilson flow and resolved gauge observables',
         reused='physical resolution context only',
         not_claimed='our transverse loop average is Wilson flow or a renormalized continuum composite field')],
    previous_results='568,574,589,598,617,637,699-701,708,712 and713 entry reused',
    no_new_universal_Casimir_or_uncertainty_theorem_claim=True),ensure_ascii=False,indent=2)+'\n')
write('round713_drafts/scope_and_dedup_review.md','''# 713主代理范围与去重审查

上一目标轮712和713入口属于进展。先读导航、报告、结果、核验与当前进程；无活动Python。回查700—701、699实际见证、617切边和637电能源，检索旧Wilson仪器／横向平均，未将一般Casimir或Luders乘积规则重新计为发现。

检查：原L字符下降到商群；弱中心holonomy不是瞬子；实际细化不等于形式切边；所有原质量和势保留；非对角电项只在双交换子中因方向支撑消去；collective instrument不是独立测量后平均；算符范数为最坏态注能而不是所有态的下界；精确中心是配置点，另有真正正常Gauss态；任意正常混合态电能界来自原全梯度分部积分；裸电能不等于参考减除后的制备功；固定宽度不保证Wilson重整化、参考态连续或699反射型结论。

三组初始结果保存在first_loop_results.json；补入第四组状态资源后保存完整结果。原Haar积分用两档规则及Casimir分部积分校对。二十式、四组检查，未计入口重复定理。未做独立代理、图像检验、安装、新任务、定时或目标变更。
''')
review='713主代理最终审查完成；没有独立代理审核。\n'
for name in ('research_note_713.md','joint_holonomy_readout_scale.py','joint_holonomy_readout_scale_results.json','unified_physics_condition_ledger_713.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round713_drafts/final_review.txt',review)
print('Prepared713 records and714 continuation.')
