"""Prepare734 causal reference/renormalized-response scope and publication."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_733.md').read_text('utf8')
ledger='# 联合条件总账：734同一过去准备的退迟来源与完整接触条件\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[732全账](unified_physics_condition_ledger_732.md)，回填[733报告](research_note_733.md)。[结果](joint_reference_polarization_boundary_results.json)、[核验](research_round_733_checks.json)。',
 '接[733全账](unified_physics_condition_ledger_733.md)，回填[734报告](research_note_734.md)。[结果](joint_retarded_reference_response_results.json)、[核验](research_round_734_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**734当前增量：** 原同一过去参考的退迟变分与完整来源接通，采用明确光滑局部减除合同。恢复背景后的历史源为光滑同背景态差，局部歧义相消；变化中的源须保态、算符、减除和有限项。联合Ward导数含背景接触及可能局部缺陷；无异常处方未由有限矩阵签收。\n\n'+marker)
updates={'C04':'734给条件性完整参考退迟源接口；局部方向可微不替代非线性反馈适定',
 'C19':'734共同过去与保记录重准备分开；恢复背景后仍可保非零光滑历史响应',
 'C22':'734共同源导数包含态、算符、减除和有限项；Ward导数须含背景接触及局部缺陷'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 734同一过去的完整来源响应

- 原光滑紧支背景族共用过去Hadamard态；完整原Dirac／Majorana传播给退迟Cauchy比较，不重选有限记录块。
- 背景恢复后的未来区域内，完整参考历史差是同背景光滑核，全部局部有限项相消；差可无限秩，与733进行中剪切非紧性相容。
- 变化中的来源按同一光滑局部减除处方共同变分态、源算符、减除核及有限项。测试来源方向可微不保证无导数损失的反馈映射界。
- 联合Ward响应含背景守恒算符导数；若所选处方还有局部缺陷，其导数也保留。628的反常类别消失不是任意parametrix自动守恒。
- 固定同一记录时，记录后态也随背景响应；固定读口来源的算符接触不能删除。
- 三组原完整矩阵核退迟重组、12生成元交换与真实记录后态来源；没有计算重整化连续真空应力或求解Einstein方程。
- 733式(10)裸qquad的格式勘误记在734尾部；冻结原件不覆盖。

## 本轮合并与下一项

C04／C19／C22由完整共同过去响应进一步连接。完整绝对源及反馈存在仍开放。

接[735](round735_drafts/STATUS.md)：核原完整二次场的局部Ward缺陷和共同允许调整，复用628反常类别、580／599有效项及734退迟核；不重算荷表或固定图接触。旧空间、604、649／699及统一目标保持。
'''
write('unified_physics_condition_ledger_734.md',ledger)
write('round734_drafts/research_note_734_draft.md',(HERE/'research_note_734.md').read_text('utf8'))
write('round735_drafts/STATUS.md','''# 第735轮入口：完整局部Ward处方与同一参考源

接[734](../research_note_734.md)、[全账](../unified_physics_condition_ledger_734.md)。同一过去的退迟参考和完整来源导数已给条件性接口，但绝对Ward归零还要核所选局部处方。

1. 先回查531、553、580／599—601、627—632及730—734。628已经关闭原普通spin＋G表示的反常类别，不能重算荷表当新结果；迹异常不同于规范／微分同胚缺陷。
2. 核原完整Dirac／Majorana二次场及变化质量的局部减除余项。源必须来自同一处方，不能各部门独立调有限系数。
3. 用成熟局部协变／共同有效作用的工具，明确哪些恒等式可直接继承，哪些需要原场的实际对象映射或局部调整。不得把矢量Dirac单独守恒结果直接套到全部手征质量。
4. 检验同一过去参考及真实记录差是否在调整后保持同一源字典；允许局部项改变绝对背景方程，但不能冒充改变远隔时间的退迟记忆。
5. 若文献和旧工作已足以关闭这一接口，立即转反馈的真实正则性／发展问题，不为常规Ward复述或参数扫描凑轮次。
6. 原经典背景是否满足绝对量子方程、全量子Gauss、装置来源、共同连续极限仍开放。旧空间和目标保持。
''')
write('round734_drafts/literature_scope_audit.json',json.dumps(dict(sources=[
 dict(url='https://arxiv.org/html/1311.7661',locations='Section2.4; Proposition2.9; Definition2.11; background smoothness near Definition2.13',
 use='Compact smooth background changes, spinor identification, retarded comparison and derivatives.',
 caution='Original finite matrix mass/self-dual extension retains730 assumptions. No full chiral Majorana interacting gravity construction imported.'),
 dict(url='https://arxiv.org/html/1210.4031',locations='Proposition3.6; Section3 smoothness; Section4 scope',
 use='Smooth Hadamard remainder, local Wick-source background regularity and renormalization choices.',
 caution='Separate conserved stress result has constant mass/background restrictions. Joint Ward-zero condition explicitly retained rather than assumed proved.'),
 dict(url='https://arxiv.org/html/1407.1994',locations='Section5',
 use='Local chiral divergence defects arise from parametrix remainders; local adjustments differ from genuine anomalies.',
 caution='628 anomaly-class cancellation is reused, not replaced by finite-matrix traces or treated as a complete local renormalization prescription.')],
 new='Original common-past reference, returned-background smooth source difference, complete point-split derivative and differentiated joint Ward mapped together; original full matrices exhibit necessary contacts.',
 excluded='Banach derivative bounds for nonlinear self-consistency; explicit absolute stress values; automatic Ward-zero scheme; continuum numerical Einstein solution.'
),ensure_ascii=False,indent=2)+'\n')
write('round734_drafts/scope_and_dedup_review.md','''# 734范围及推导审查

主代理审查，无独立代理，无图像检验。

- 上轮有正式733和734入口证据，属于进展；本轮开始未发现运行中的Python。
- 回查620接触警示、623—625固定图来源、628反常类别、632零导数接触和730—733，未把成熟Duhamel或一般记忆效应作为新定理。
- 新接口是同一原连续参考的退迟输送、恢复背景后的全参考光滑差、完整局部减除响应及联合Ward背景接触。
- 光滑联合参数依赖要求原共同过去和光滑局部减除处方。分别Hadamard不够；测试方向可微不等于无导数损失Banach可微。
- 返回背景后的差平滑可无限秩，不与733进行中跨背景非紧性混淆。
- 保留局部缺陷A(B)，只在声明相容零缺陷处方时签收齐次绝对Ward。628反常类别消失不自动确定任意局部parametrix。
- 12个生成元的接触通过实际标量／连接切向独立差分计算，对易子只是比较对象。实际态和Majorana均保留。
- 退迟核通过每步Frechet插入重组，与逐步协方差微分独立核对。真实记录后态同时保态和算符导数。
- 有限矩阵接触不是Hadamard反项；数值不作连续物理应力预测。测量时刻补偿和实际因果装置仍受634限制。
- 733格式错误登记勘误，未改冻结文件。734式19的制表符在冻结前修复。
- 反馈公式只定位下一阶，不证明绝对基准、半经典解存在或长期稳定。目标不改，旧空间接口保持。
''')
files=('research_note_734.md','joint_retarded_reference_response.py','joint_retarded_reference_response_results.json',
 'unified_physics_condition_ledger_734.md')
write('round734_drafts/final_review.txt','734 primary review; no independent agent review.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
 'Three checks pass. Declared smooth same-past reference response and joint background contacts. Local anomaly remainder retained; no absolute semiclassical closure claim. Goal remains open.\n')
v=(HERE/'verify_round733.py').read_text('utf8')
v=remap(v,{'joint_reference_polarization_boundary':'joint_retarded_reference_response',
 'joint_relative_source_development':'joint_reference_polarization_boundary',
 'round734':'round735','round733':'round734','round732':'round733',
 '_733':'_734','_732':'_733','range(584,733)':'range(584,734)',
 '3256':'3270','3270':'3284','round=733':'round=734','round=732':'round=733',
 '3404':'3407','1511':'1514','reference_embedding_entry':'causal_reference_entry','==18':'==20'})
write('verify_round734.py',v.replace('Reproduce733','Reproduce734'))
post=(HERE/'postcheck_round733.py').read_text('utf8')
post=remap(post,{'range(584,734)':'range(584,735)','round733':'round734','_733':'_734',
 '3270':'3284','第733':'第734','round=733':'round=734'}).replace('Check733','Check734')
write('postcheck_round734.py',post)
pub=(HERE/'publish_round733.py').read_text('utf8')
pub=remap(pub,{'734':'735','733':'734','732':'733','3404':'3407','3401':'3404',
 '1511':'1514','3270':'3284','## 379.':'## 380.','## 284.':'## 285.',
 'joint_reference_polarization_boundary':'joint_retarded_reference_response',
 '有限记录模式的非线性边界与完整参考响应':'同一过去准备的退迟来源与完整接触条件',
 '旧空间合同保持，完整参考响应不由有限记录替代':'旧空间合同保持，退迟来源与接触须共同变分'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第734轮完成：** [同一过去准备的退迟来源与完整接触条件]({p}research_note_734.md)原完整参考的退迟输送、恢复背景后的光滑源差与共同减除导数接通；联合Ward响应须含背景接触及可能局部缺陷。三组、二十式通过，最新734／3407，1514份编号科学文件、3284份保护证据。[核验]({p}research_round_734_checks.json)、[全条件账]({p}unified_physics_condition_ledger_734.md)。绝对Ward处方与非线性反馈仍开放。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（734后，优先于下方历史安排）：** 接[735完整局部Ward处方与参考来源]({p}round735_drafts/STATUS.md)，复用628反常类别、580／599有效项及734退迟响应，核同一允许局部调整；不重算荷表或固定图接触。旧空间、604、649／699及统一目标保持。'"
write('publish_round734.py','\n'.join(lines)+'\n')
