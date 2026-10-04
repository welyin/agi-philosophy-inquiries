"""Prepare706 publication; preserve frozen science and the parity addendum."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k)
          for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

names=('unified_physics_condition_ledger_706.md','round706_drafts/research_note_706_draft.md',
       'round706_drafts/final_review.txt','round707_drafts/STATUS.md',
       'round706_drafts/local_domain_entry.py','round706_drafts/local_domain_entry_results.json',
       'round706_drafts/local_domain_entry.md','round706_drafts/check_and_publish_entry.py',
       'round706_drafts/entry_checks.json','round706_drafts/prepare_entry_publication.py',
       'round706_drafts/literature_scope_audit.json','round706_drafts/parity_completion_705.md')
protected=2870+3+len(names);assert protected==2885
ledger=(HERE/'unified_physics_condition_ledger_705.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：705区域压缩、边界电荷与Gauss支持',
    '# 联合条件总账：706保边界局部近似与真实记录来源')
ledger=ledger.replace('接[704全账](unified_physics_condition_ledger_704.md)，回填[705报告](research_note_705.md)。[结果](joint_regional_charge_compression_results.json)、[核验](research_round_705_checks.json)。',
    '接[705全账](unified_physics_condition_ledger_705.md)，回填[706报告](research_note_706.md)。[结果](joint_local_history_source_limit_results.json)、[核验](research_round_706_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**706当前增量：** 保边界、按CAR宇称分块的局部重数通道，在原完整固定图的真实演化和记录后插入。来源无关的双向图范数接623—625，704夹权原热导数接准备／动力混合项，六类总阶数二标量响应及零阶cq态共同收敛。过程仍无限维，使用原H及原Gibbs准备，不是新的有限有效H，也不证明空间细化。705宇称补充另存，冻结历史不覆盖。\n\n## 当前共同对象及仍存在的分支')
updates={
'C01 量子对象':'706局部偶Kraus保原边界／Gauss，CAR宇称补充已核；整体无限维',
'C02 区域组合':'706同一严格局部族接原含跨边相互作用的真实历史',
'C03 事件记录':'706原有限历史的cq后态及总阶数二标量来源共同保留',
'C04 内部演化':'706各等待均用原完整H，不将分区误作独立演化',
'C19 参考态':'706共用原随准备参数变化的Gibbs态，记录后不重置为热态',
'C20 尺度映射':'706固定图保边界局部正近似已有实际来源极限；不是空间细化或整体有限化',
'C22 来源反作用':'706双向图范数及原夹权热导数共同给准备／动力混合项，来源Hessian保原身份'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|原固定图完整量子过程|'):
        c=row.split('|');c[2]+='；706保边界局部CP族与原真实历史、准备／动力来源共同极限';rows[i]='|'.join(c)
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            c=row.split('|');c[2]+='；'+value;rows[i]='|'.join(c);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';a=ledger.index('## 本轮合并与下一项')
ledger=ledger[:a]+'''## 706局部区域、真实历史与共同来源

- 705原按lambda回填补充为按lambda及局部CAR宇称回填，Kraus为偶；sigma不充当新边界荷，也不强制两侧宇称一致。见[补充](round706_drafts/parity_completion_705.md)，旧冻结科学文件保留。
- 706入口的固定辅助能量标签给局部等距及伴随的精确比较能源交织；正向和伴随图范数均已核。辅助能量不是物理H新增项或免费资源。
- 623原算符域与637区域比较R等价；原U通常不与R对易，但图范数等价足以给固定时间界。不能以单边形式下界冒充此步。
- 各次原记录后插入局部Phi，所有等待用原完整H，两CTP支共享辅助表示，初态为原Gibbs。跨边耦合未删，记录后不重热化，整个过程仍无限。
- 625双来源弱形式证明可接新局部等距：较晚来源移到bra，来源后不要求闭域。零、一、二阶标量形式共同收敛，不新增D(H²)要求。
- 704原A夹权准备导数经RA^-1输送到R夹权迹类；与上述弱算符极限配对，六类总阶数二来源共同匹配，包括准备／动力混合。
- 零阶cq态直接复用592／705；没有领取完整动态cq二阶迹范数可微、任意来源函数空间一致余项或仅凭初始尾的普适来源截止率。
- 数值用两节点256维条件径向／中性CAR诊断，保原测地边势、Majorana和实际跳跃；冻结链路不当正常全Gauss态。完整简并后只有两个不同局部近似；概率准不代表动力来源准。
- 旧空间382—386、425、522—524直接复用。384删去Lipschitz保持；386／425替代；523共同实现及524局域UV探针不重新列为缺口。实际端点／参考／仪器与本h、s、CAR／Gauss尺度身份仍待接。

## 本轮合并与下一项

C02／C03／C04／C19／C20／C22固定图的局部来源接口接通。保无限边界载体并使用原H是本族的范围，705有限容量权衡不撤销；不对空间极限、面积律、手征或量子GR销账。

接[707](round707_drafts/STATUS.md)：回到真实共同空间尺度，先回查637—646、663—667，区分区域合并、图细化、谱精度。701入口的对角精度选择不重复，四分支及699范围不变，先整合条件、认知设计后置。目标不变。
'''
write('unified_physics_condition_ledger_706.md',ledger)
mapping={'joint_regional_charge_compression':'joint_local_history_source_limit',
         '704':'705','705':'706','706':'707',
         '3315':'3318','3318':'3320','1424':'1427','1427':'1430','2855':'2870','2870':'2885'}
verify=replace((HERE/'verify_round705.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_preparation_history_limit as prior_model',
                      'import joint_regional_charge_compression as prior_model')
verify=verify.replace('regional charge support, finite-output tradeoff and local multiplicity channels.',
                      'local boundary-preserving histories and joint source jets.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==20","text['display_formulas']==18")
verify=verify.replace('==(3,0,0)','==(2,0,0)').replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=2,')
verify=verify.replace("'round706_drafts/attribution_correction_363.md'","'round705_drafts/attribution_correction_363.md'")
write('verify_round706.py',verify)
pub=replace((HERE/'publish_round705.py').read_text('utf8'),mapping|{'## 351.':'## 352.','## 256.':'## 257.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第706轮完成：** [保边界局部近似与真实记录来源]({p}research_note_706.md)'
         '同一保边界局部CP族接通原真实历史、准备／动力总阶数二来源及零阶cq态；'
         '仍使用原完整H及原Gibbs，整体无限维。'
         '两组、十八式通过，最新706／3320，1430份编号科学文件、2885份保护证据。'
         '[核验]({p}research_round_706_checks.json)、[全条件账]({p}unified_physics_condition_ledger_706.md)。'
         '705宇称补充已核，空间共同尺度仍开放，旧空间与目标不变。')
order=('**当前执行顺序（706后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[707实际共同空间尺度]({p}round707_drafts/STATUS.md)，'
       '回查637—646、663—667，区分区域合并、图细化与固定图谱精度；'
       '保留四分支及699范围，不重复固定图热迹、局部矩或容量反例。')
'''+pub[b:]
pub=pub.replace('区域压缩、边界电荷与Gauss支持','保边界局部近似与真实记录来源')
pub=pub.replace('区域、边界支持及有限容量共同验收','局部区域、真实记录与来源共同验收')
pub=pub.replace('旧空间接口继承，内部边界表示不当空间维数','旧空间接口继承，局部来源极限不当空间细化')
write('publish_round706.py',pub)
write('postcheck_round706.py',replace((HERE/'postcheck_round705.py').read_text('utf8'),mapping))
with (HERE/'round706_drafts/research_note_706_draft.md').open('xb') as f:
    f.write((HERE/'research_note_706.md').read_bytes())
review="""706 primary-agent mathematical/code/scope review; no independent agent.
No goal status or definition changed. Existing formal705 and706 entry make progress.
705 parity omission resolved in separate frozen addendum; no historical file changed.
Kraus reset within lambda and local CAR parity, identity on all boundary carriers.
Source-independent auxiliary energy construction and both adjoint graph limits reused.
Original H graph norm equivalent to regional comparison R via623/637; not inferred from one-sided form.
U need not commute R; uniform graph bounds and nonnegative commuting auxiliary energy suffice.
Actual local CP maps inserted after records; original full H and original Gibbs always retained.
Entire process infinite dimensional; not625 finite global process, not703 logarithmic transfer.
One insertion converges in Hilbert space; later source moved to bra for weak second forms.
No unproved source-after-insertion domain preservation, common D(H^2), or strong second derivatives.
Finite fixed source menu/time/history; uniform constants do not imply spatial uniformity.
704 weighted original preparation derivatives paired with bounded WOT form representatives.
Scalar joint jets total order2; no full cq dynamic second trace differentiability claimed.
Equal-history outcomes positive >=4^-q; different histories give complex CTP coefficients.
Zero-order cq convergence inherited592/705, not claimed as a new independent theorem.
Conditional2node256dim diagnostic freezes link/orientation; not full Gauss spectrum.
Original geodesic edge and neutral hopping retained, with JordanWigner CAR sign.
Local parity can change under hopping, while total parity and local Kraus parity are exact.
Original conformal powers T-6,W+6,V+2,B0,J-2 and Hessian inherited667.
Complete degeneracies mean four requested cutoffs yield only two distinct regulators.
Large intermediate source errors openly reported; zero full-cutoff difference only structural.
Independent exp finite differences validate mixed coefficient; equal-history normalization checked.
705 capacity no-go,699 candidate scope and all four branches remain.
Old382-386,425,522-524 retained; removedLipschitz not restored;386/425 alternatives.
No apparatus design, free reset, spatial continuum, Einstein action or quantumGR claimed.
Two scientific groups,18 formulas. Next707 returns to actual shared spatial mapping.
"""
for name in ('research_note_706.md','joint_local_history_source_limit.py','joint_local_history_source_limit_results.json',
             'unified_physics_condition_ledger_706.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round706_drafts/final_review.txt',review)
print('706 prepared; protected evidence',protected)
