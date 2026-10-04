"""Publish the full742 mathematical scope without inventing a new device."""
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

note=(HERE/'round742_drafts/research_note_742_working.md').read_text('utf8')
note=note.replace('# 742工作报告：','# 第742轮：',1)
intro=note.index('日期：');end=note.index('\n\n## 1.',intro)
note=note[:intro]+'''日期：2026-10-04。接[741](research_note_741.md)，保留[已执行入口](round742_drafts/research_note_742_working.md)。[代码](joint_record_gaussian_boundary.py)、[结果](joint_record_gaussian_boundary_results.json)、[核验](research_round_742_checks.json)、[条件账](unified_physics_condition_ledger_742.md)。一组原对象检查、十二式；主代理审查，无独立代理或图像检验。入口检查首次计入，正式累计742／3422。

**本轮完成的是一条有限定范围的实现反证：指定的纯二次费米扩张不能产生原实际非选择记录，且误差有正下界。未完成的是原完整量子相互作用中的内部局域实现。** 不把这二者或整个统一目标混为一谈。'''+note[end:]
note=note.replace('../','')
note=note.replace('目标不变；空间382—386、425、522—523及604、649／699原范围全部保持。',
    '接[743原量子相互作用与同一记录](round743_drafts/STATUS.md)。目标不变；空间382—386、425、522—523及604、649／699原范围全部保持。')
write('research_note_742.md',note)
write('round742_drafts/research_note_742_draft.md',note)

ledger=(HERE/'unified_physics_condition_ledger_741.md').read_text('utf8')
ledger='# 联合条件总账：742原记录与二次实现的严格边界\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[740全账](unified_physics_condition_ledger_740.md)，回填[741报告](research_note_741.md)。[结果](joint_absolute_source_development_results.json)、[核验](research_round_741_checks.json)。',
    '接[741全账](unified_physics_condition_ledger_741.md)，回填[742报告](research_note_742.md)。[结果](joint_record_gaussian_boundary_results.json)、[核验](research_round_742_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**742当前增量：** 原纯准自由参考上的非选择占据记录有四点Wick缺陷4p(1−p)；任何独立高斯探针、确定偶二次费米演化再丢弃探针的实现，输出迹距离至少2p(1−p)/7。此限制针对当前二次实现类，保741双线性来源及原量子玻色相互作用候选，不把成熟FLO的原始测量纳入排除。\n\n'+marker)
updates={'C03':'742原实际记录不属确定二次高斯扩张；限制类误差有正下界，完整相互作用实现仍待核',
         'C19':'742同两点高斯替代丢真实四点关联；参考、记录与高阶来源须共用同一后态',
         'C22':'742保741平均源接口；内部装置须超出背景二次扩张，未证明一定添加新物种'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 742原记录对内部实现的限制

- 原纯协方差提取已有伙伴模式，给非高斯四点缺陷和统一迹距离下界，不增加粒子。
- 633连续径向包p=1/2直接提供精确例；730变化准备的原128维对象独立校准，未冒充连续数值证明。
- 实现类包括任意二次耦合时间及任意独立高斯辅助费米模式；不包含外部随机控制、已给占据测量、非高斯资源或量子玻色相互作用。
- 741平均双线性来源保持，但不能把相同两点的高斯态当成同一完整记录后态。
- 不把一般Gaussian封闭作为原创定理，不重复633因果和634注能反例。

## 本轮合并与下一项

C03／C19／C22共同要求保留足以产生真实记录关联的内部量子过程。新增的是明确实现类的排除与定量间隔，未排除整个原模型。

接[743](round743_drafts/STATUS.md)：先查原量子Yukawa／Majorana及规范相互作用是否提供所需条件；实际读口、共同来源、因果过程及准备必须一起核，先不添加新装置。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_742.md',ledger)
write('round743_drafts/STATUS.md','''# 第743轮入口：原量子相互作用与同一记录

接[742](../research_note_742.md)、[条件账](../unified_physics_condition_ledger_742.md)。确定二次费米演化加独立高斯探针不能产生原实际记录；量子玻色—费米耦合尚未被排除，不能因此声称必须添加新物种。

1. 回查598原完整作用、602／604原量子来源、614／620相互作用来源、633—634原sterile记录、716实际仪器及735二次近似。
2. 区分原质量双线性及其量子标量系数、实际占据读口和必要的测量通道；超出高斯类仅是解除742限制，不等于目标操作已实现。
3. 先整合原耦合和保留关联，核其是否能产生原目标记录或一个明确不同的有效仪器；准备、传播、来源交换与末端记录共同入账。
4. 不用新四费米或指针耦合直接冒充原作用的推论；若需有效项，列出从原过程取得它的条件及误差。
5. 复用Fewster—Verch局域过程范围、575已有自主指针与633因果反例，不重新编号这些旧结果。不要只做高斯误差下界的常数优化。
6. 741首阶平均来源、空间382—386／425／522—523、604及649／699的范围保持。统一目标未改变。
''')
write('round742_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/quant-ph/0404180',locations='SectionsII,IV,V,VIII',
    use='Quadratic CAR transformations and Gaussian Wick closure; FLO also assumes occupation measurement and feedback as primitives, excluded from our restricted implementation class.'),dict(
    url='https://arxiv.org/html/1810.06512v3',locations='Section3 and conclusion',
    use='Local coupling and scattering define instruments under their assumptions; no arbitrary target-instrument realization in original matter or gravitational backreaction theorem.')],
    new='Actual original pure-reference partner, record fourth moment, and positive trace-distance obstruction for a restricted internal implementation.',
    inherited='633 nonGaussian warning and exact p=1/2 continuum example;634 causal/source restrictions;730 pure reference;735 quadratic scope;741 two-point source.',
    excluded='All FLO including primitive measurements, quantized original boson-fermion interactions, universal model impossibility, autonomous local implementation, new physical dimensional selection.'),ensure_ascii=False,indent=2)+'\n')
write('round742_drafts/scope_and_dedup_review.md','''# 742证明与范围审查

主代理审查，无子代理、无图像检查。工作入口保留；仅一组检查首次入正式累计，不将相同结果重复计轮次。

- Γ²=−I给伙伴正交、反对易及−A伙伴块；无需重新选择参考或添加物种。
- 非选择记录保四点、删两点交叉；真实后态四点1而Wick预测a²。p=0或1为退化例外，未纳入排除。
- 任意高斯候选四点差≤2δ；三个Pfaffian乘积差各≤4δ，得到14δ界。未宣称最优，数值代入不是区间认证。
- CAR乘积、确定二次演化、限制共同保高斯。即使不要求局域，指定类仍受限。
- Bravyi的FLO含原始测量和反馈；本轮没有排除它的定义。随机I／R分解也不自动给内部随机控制实现。
- 原物质量子标量乘费米双线性不属于这类纯费米二次Hamiltonian；下一轮应审查原耦合，不能把必要条件变成新装置需求。
- 633已说后态不必高斯，新增是实际伙伴、严格缺陷及正近似误差；连续p=1/2直接复用。有限128维校准不代替连续计算。
- 741只用后态双线性来源，其构造保持；噪声和多记录需要真实高阶关联。
- 保旧空间各桥接与目标；没有推出或否定全部量子引力。
''')
main=('research_note_742.md','joint_record_gaussian_boundary.py','joint_record_gaussian_boundary_results.json','unified_physics_condition_ledger_742.md')
write('round742_drafts/final_review.txt','Primary-agent review only. Restricted quadratic-dilation obstruction, not a no-go for the complete interacting model.\n'+
    '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'742':'743','741':'742','740':'741','3421':'3422','3419':'3421',
         '1535':'1538','3373':'3387','3359':'3373','387':'388','292':'293',
         'joint_absolute_source_development':'joint_record_gaussian_boundary',
         'absolute_source_order_entry':'gaussian_record_entry'}
publication=remap((HERE/'publish_round741.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第742轮完成：** [原记录与二次实现的严格边界]({p}research_note_742.md)同一纯参考的真实记录有非高斯四点关联，确定二次费米演化加独立高斯探针的误差至少2p(1−p)/7；633连续例为1/14。仅排除该实现类，保741平均来源及完整量子相互作用候选。一组、十二式通过，最新742／3422，1538份编号科学文件、3387份保护证据。[核验]({p}research_round_742_checks.json)、[条件账]({p}unified_physics_condition_ledger_742.md)。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（742后，优先于下方历史安排）：** 接[743原量子相互作用与同一记录]({p}round743_drafts/STATUS.md)，先查原Yukawa／Majorana耦合、实际读口和完整来源；不因高斯限制就添加新物种。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('实际绝对来源、联合初值与一阶共同反作用','原实际记录与二次实现的严格边界')
publication=publication.replace('旧空间合同保持，一阶绝对发展不替代内部准备','旧空间合同保持，平均来源不替代完整记录')
write('publish_round742.py',publication)
write('postcheck_round742.py',remap((HERE/'postcheck_round741.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round741.py').read_text('utf8'),mapping)
verification=verification.replace('==(2,0,0)','==(1,0,0)').replace('run=2,failures=0','run=1,failures=0')
verification=verification.replace("checks['display_formulas']==18","checks['display_formulas']==12")
verification=verification.replace('original_complete_source_hierarchy_and_reference_lift_checked=True','original_reference_record_partner_and_positive_Gaussian_gap_checked=True')
write('verify_round742.py',verification)
print('Prepared742 complete restricted-implementation result and743 original-coupling scope.')
