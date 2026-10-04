"""Prepare731 without overwriting frozen artifacts or changing the goal."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_730.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：731原物质上的联合初始约束补偿\n'+rest
ledger=ledger.replace(
 '接[729全账](unified_physics_condition_ledger_729.md)，回填[730报告](research_note_730.md)。[结果](joint_dynamic_continuum_reference_results.json)、[核验](research_round_730_checks.json)。',
 '接[730全账](unified_physics_condition_ledger_730.md)，回填[731报告](research_note_731.md)。[结果](joint_source_constraint_response_results.json)、[核验](research_round_731_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**731当前增量：** 原非平坦配置的电弱Gauss椭圆核为零，颜色保留总荷条件；原径向与曲率电动量菜单共同补足三个总动量。给定光滑源且C为正时，全部补偿能量进入新的唯一正共形解。一阶原量子来源有条件性初始接入口；数值是声明源，非同态自洽量子解。背景／参考响应、绝对基准及全时间发展仍开放。\n\n'+marker)
updates={
 'C04':'731给定源的原Gauss与Einstein初始约束共同可解，不等于完整量子反馈时间发展',
 'C19':'731原730全局颜色不变参考及sterile后态满足颜色源条件；实际绝对源仍须统一减除',
 'C22':'731全部电弱局部荷、三方向动量及能源共同补偿，保代价；完整自洽来源泛函及高阶有效项未关闭'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 731原初始约束的共同来源接入口

- 570核空间合同直接复用；原非平坦弱场消去电弱稳定子，得到全点Gauss右逆，颜色固定配置仍须零总荷。
- 三菜单保Gauss且张成三个总动量，补偿后用原572向量方程；有源玻色部分始终用协变动量，不漏A乘荷项。
- 原标量、规范配置及CMC保持，改变原共轭动量和引力自由资料；所有标量／电动能包括颜色二阶项留在原Hamiltonian。
- C减2ερ保持正时，572上下解与最大比值论证给新的唯一正ψ；新线性化正超解给来源的光滑初始响应。
- 原730参考与sterile记录的全局颜色不变性使其局部颜色平均源为零。任一确定的同态光滑来源可放入一阶初始系数。
- 15³配点只验证声明源、完整约束和响应，不冒充730连续应力计算、不证明绝对量子背景或有限强度自洽。
- 580／600高阶项、730准备响应、634守恒接口及全部旧空间限定结论保持。

## 本轮合并与下一项

C04／C19／C22在给定光滑初始源的范围内进一步合并。实际同态来源与修正背景的共同闭合未完成。

接[732](round732_drafts/STATUS.md)：优先核同一来源在变化背景所需的完整资料、参考处方、Ward及有效项，明确可受控展开的范围。停止补偿菜单、ε精度和椭圆求解器优化。旧空间、604、649／699及统一目标保持。
'''
write('unified_physics_condition_ledger_731.md',ledger)
write('round731_drafts/research_note_731_draft.md',(HERE/'research_note_731.md').read_text('utf8'))
write('round732_drafts/STATUS.md','''# 第732轮入口：实际同态来源与修正背景的共同闭合

接[731](../research_note_731.md)、[全账](../unified_physics_condition_ledger_731.md)。给定光滑源的原Gauss／三动量／完整Hamiltonian初始约束已共同可解；不能再用扫描ε或外给正能量增加轮次。

1. 先回查325、570—574、580—583、598—604、630—635、649—651、704及730—731。旧源跳接、准备响应、一般重整化常数不再重报。
2. 实际同态来源必须同时给荷、动量、能量、空间应力和标量质量力；共同Ward与初始资料依赖先明确。
3. 731初始右逆只对给定源成立。有限强度需在新背景重建同一参考／仪器及其源，不能把旧Hadamard态的数值常量原样搬入。
4. 区分原经典基线、相对记录源、绝对真空源及统一有限反项。背景变化可能引入更多时间导数；单一ψ初始不自动决定全部重整化源。
5. 可检验路线包括明确阶次的共同扰动展开、固定重整化处方的状态依赖连接，或有完整资料的局部半经典发展。不能把抽象可微假设当已经验证。
6. 先给真实共同接口的新结果，暂不设计自主认知装置；不把原输入的四维Einstein／SM解释成已由认知原则推出。旧空间和统一目标保持。
''')
write('round731_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(url='https://arxiv.org/html/2106.15027v2',
   use='Meta-Theorem1.1 and sections5-6: canonical matter data, matter constraints and CMC conformal framework.',
   caution='Classical structural framework, not a renormalized Dirac or self-consistent semiclassical existence theorem. Original equations checked explicitly here.')],
 inherited='570 stabilizer and Fredholm compatibility; 572 CMC Fourier solve/barriers/uniqueness; 580 and600 effective higher-derivative terms; 730 Hadamard/source response.',
 new='Original nonflat EW kernel removal, joint charge/momentum/energy repair, initial linear response and actual color-invariant source applicability.',
 excluded='Computed continuum quantum stress, absolute semiclassical solution, quantum Gauss projection, autonomous preparation or emergence of GR.'
),ensure_ascii=False,indent=2)+'\n')
write('round731_drafts/scope_and_dedup_review.md','''# 731范围、证明及代码审查

主代理审查；无独立代理、无图像检验。

- 电弱G的伴随号、Higgs五实切向和原D=partial−A cross逐项核对。原非平坦Ax使剩余电磁稳定子消失，颜色常数核仍保留。
- 补Gauss后单独玻色canonical与covariant不能混用；代码从p Dphi与E F直接算全动量。
- 入口第三方向与572径向两方向直接复用，列满秩后所有总动量消去。不称菜单唯一、最优或免费动作。
- 原标量动能、弱／圆／颜色电能、磁能均留在新约束。源码旧残差不含ερ，故新增独立物理残差并验证漏项会失败。
- 连续标量存在／唯一复用572；新的正超解身份给响应可逆，未把数值零阶势正性当额外一般公理。
- 配点选奇数15排除Nyquist一阶导数核；数值不证明连续误差。三幅度是非线性源、符号和解析导数校准，不开展无尽扫描。
- 代码来源明确外给。解析同态应用限既定光滑基背景来源的一阶约束系数；全局颜色不变性的证明与原730构造一致。
- 绝对基准／同背景相对源／变化背景自洽四者分开。Hadamard单态正则不足推出整个来源泛函的无导数损失Banach光滑性。
- 580／600的高阶项、Ward及空间应力／质量力保留到下一项。一般有效作用或约束工具不是本项目独创定理。
- 旧382—386、425、522—523直接按667／704继承，不将输入三维或Einstein当新生成结论。
''')
files=('research_note_731.md','joint_source_constraint_response.py','joint_source_constraint_response_results.json',
       'unified_physics_condition_ledger_731.md')
write('round731_drafts/final_review.txt','731 primary review; no independent agent review.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
 'Three checks pass. Original initial constraints jointly repaired with complete energy; numerical inputs are declared sources. Actual first-order interface is conditional, full state-dependent feedback and unified goal remain open.\n')

v=(HERE/'verify_round730.py').read_text('utf8')
v=remap(v,{'joint_dynamic_continuum_reference':'joint_source_constraint_response',
 'joint_curved_periodic_reference':'joint_dynamic_continuum_reference',
 'round731':'round732','round730':'round731','round729':'round730',
 '_730':'_731','_729':'_730','range(584,730)':'range(584,731)',
 '3214':'3228','3228':'3242','==18':'==20',
 'round=730':'round=731','round=729':'round=730',
 '3395':'3398','1502':'1505','varying_reference_entry':'gauge_counterflow_entry'})
v=v.replace('Reproduce730','Reproduce731')
write('verify_round731.py',v)
post=(HERE/'postcheck_round730.py').read_text('utf8')
post=remap(post,{'range(584,731)':'range(584,732)','round730':'round731','_730':'_731',
 '3228':'3242','第730':'第731','round=730':'round=731'}).replace('Check730','Check731')
write('postcheck_round731.py',post)
pub=(HERE/'publish_round730.py').read_text('utf8')
pub=remap(pub,{'731':'732','730':'731','729':'730','3395':'3398','3392':'3395',
 '1502':'1505','3228':'3242','## 376.':'## 377.','## 281.':'## 282.',
 'joint_dynamic_continuum_reference':'joint_source_constraint_response',
 '原动态背景上的连续参考、真实记录与同一来源':'原物质上的联合初始约束补偿与量子来源接口',
 '旧空间合同保持，背景连续参考不替代自洽引力':'旧空间合同保持，初始约束补偿不替代完整自洽发展'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第731轮完成：** [原物质上的联合初始约束补偿与量子来源接口]({p}research_note_731.md)原非平坦场允许全电弱Gauss及三方向动量共同补偿，全部代价放回新的正共形解；颜色条件和同态一阶范围明确。三组、二十式通过，最新731／3398，1505份编号科学文件、3242份保护证据。[核验]({p}research_round_731_checks.json)、[全条件账]({p}unified_physics_condition_ledger_731.md)。实际绝对量子源、自洽反馈与全时间发展仍开放。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（731后，优先于下方历史安排）：** 接[732实际同态来源与修正背景]({p}round732_drafts/STATUS.md)，复用原初始右逆，核完整源、参考处方、Ward及高阶资料的共同闭合。停止补偿菜单、ε和椭圆精度优化；旧空间、604、649／699及统一目标保持。'"
write('publish_round731.py','\n'.join(lines)+'\n')
