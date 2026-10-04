"""Prepare761: explicit relative scaling and a finite-time joint connection."""
import hashlib
import json
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_760.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：761相对阶次、有限时间量子源与实际反作用',1)
ledger=ledger.replace(ledger.split('\n')[2],
 '2026-10-04。接[760全账](unified_physics_condition_ledger_760.md)，回填[761报告](research_note_761.md)。[结果](joint_loop_scale_transport_results.json)、[核验](research_round_761_checks.json)。总目标保持，新增缩放另列为条件分支。',1)
updates={
 'C01':'761新增Hb,epsilon+epsilon B族在小epsilon保完整有限Gauss纤维及任意被动参考；未选单带',
 'C03':'761相同量子玻色基线比较下，原有限sin s历史的首阶概率由同源响应决定；内部终端仍缺',
 'C04':'761只在新增相对阶次内有固定T输运及O(sqrt epsilon)误差，原B主阶族的759限制保持',
 'C19':'761同一矩阵沿原标量背景流运输，源与Gauss余额相容；空间gamma仍外给',
 'C20':'761相对能源阶次不是空间网格细化，不提供epsilon=1误差或连续重整化映射',
 'C21':'761完整量子纤维保留，709中性荷压缩可兼容；不称所有内部信息已经典化',
 'C22':'761精确Duhamel给费米诱导响应和报告共同首阶，不能仅由O(sqrt epsilon)态差除epsilon领取系数'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'''

## 761新增条件分支K0-L：不能隐藏的相对阶次

本分支取H_epsilon=Hb,epsilon+epsilon B，并用i epsilon d/dt生成真实时间；原726/758—760族为Hb,epsilon+B。两者仅在epsilon=1同为原H，小参数路径不同。保全部B系数和边不等于相对相互作用强度未改变。新增阶次是工作输入，未从H1—H3推出。

固定原图、Gauss包、有限纤维和紧相时间管内，主符号为标量，故无需B能带隙；同一完整矩阵给有限时间输运及来源矩。由精确Duhamel得到相对于相同量子玻色基线的首阶背景响应，原sin s有限历史有同一个首阶概率系数。原合法中性相干对给非零后果。

这是物质/规范背景响应，不是动态时空度规。小参数下的条件性构造不证明epsilon=1近似好，不证明对应真实宇宙尺度，也不删除原K0及759/760各自边界。741平均源形式阶次可以对照，但固定图/连续场、绝对参考、重整化和全部玻色修正仍须共同接入。

## 下一项：检验阶次依据及同一尺度的来源

接[762入口](round762_drafts/STATUS.md)。不继续局部振子数值精度；先把相对阶次选择、原有限模式态及连续源映射放在同一条件表。能从已有资源/占据缩放推出多少，哪些仍是物理输入，逐项辨别。旧空间及604、649/699限定结论保持。
'''
write('unified_physics_condition_ledger_761.md',ledger)
write('round761_drafts/research_note_761_draft.md',(HERE/'research_note_761.md').read_text('utf8'))
write('round762_drafts/STATUS.md','''# 第762轮入口：相对阶次依据与同一尺度来源

接[761](../research_note_761.md)、[条件账](../unified_physics_condition_ledger_761.md)及[共同候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 761新增K0-L缩放Hb,epsilon+epsilon B；原主阶B族保留，759/760没有被推翻。只在epsilon=1相同，不得将不同路径偷换。
2. 新分支在固定图/紧相时间管有完整矩阵输运和同源首阶物质/规范响应；原有限sin s历史也有首阶系数。来源及报告的系数来自精确Duhamel，不是粗态误差除epsilon。
3. 参数的资源或物理尺度解释仍输入；H1—H3未推出它，epsilon=1误差未验。不是改变某条边得到容易的玩具结论，全部原B在解析构造中保留；振子只是局部校准。
4. 下一项先回查353—357、525、741、756、758、709和704/706；把已有资源/占据缩放与本分支的相对能源/时间计数对齐，避免另做普通平均场玩具。
5. 真正下一接口是同一原态、实际有限任务及连续来源/受约束几何的共同尺度；有限图输运不自动给空间细化或连续手征。
6. 如尝试将epsilon与741的lambda认同，必须区分平均源计数一致与全部态/算符/重整化映射已证明。保完整量子相干、同一参考与来源噪声，不另拟合噪声。
7. 原空间gamma仍外给；761线性化Gauss余额不能替代Einstein全部初始/演化约束。绝对玻色基准及同阶量子修正继续列账。
8. 不扫局部振子精度，不重做已完成压缩；不改目标、不建任务/调度、不做图像。旧382—386、425、522—523及604、649/699原边界保持。
''')
write('round761_drafts/scope_and_dedup_review.md','''# 761范围与去重审查

- 上一目标轮为进展：760完整编号交付、检查与出版；761工作回查也已保存。本轮读三导航、760全文与结果、入口和回查，无运行Python。
- 709已有真压缩；704/706已有固定图过程/来源收敛。353—357有其他平均场模型，525有有限自反馈模型；均不重跑或冒充当前全图。
- 新关键是参数合同不同：旧只缩玻色量子化而B主阶；741另以形式参数乘费米来源。本轮明确新增Hb,epsilon+epsilon B路径，只在epsilon=1等于旧H。
- 全部原B系数、规范、质量及边保留不等于相对强度未变。额外阶次列为工作输入，未由认知原则推出；不把旧759反例撤回。
- 原自伴过程、Gauss初值、紧相管及有限纤维给局部符号估计；主符号全标量，所有内部矩阵属于一个块，无需B能带隙。无界源需图权和微局部尾控制。
- 先用残差O(epsilon^(3/2))证明状态输运O(sqrt epsilon)。首阶响应另用精确Duhamel及二阶标量符号控制，不能拿状态误差除epsilon。
- 原有限记录历史用逐段伸缩差保时间次序。源系数和记录系数来自同一矩阵；非零例为原sterile相干对，不是757同二点tau族。
- Gauss分量在代表框架中的余额与完整群平均期望区分；非中心固定分量不是物理标量观测。
- 数值：精确算术身份、原32模式来源、明示的局部二次玻色/原中性CAR校准。未解原完整图轨道，更未解连续Einstein。
- 本轮关闭一个新增尺度分支的固定时间连接；相对尺度依据、epsilon=1准确性及Q到E仍开放，不改总目标。
''')
write('round761_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(title='A semiclassical Egorov theorem and quantum ergodicity for matrix valued operators',
 url='https://arxiv.org/html/math-ph/0204018',
 checked='Section3 equations3.3-3.8 and theorem3.2; scalar principal symbol case has one full matrix block, subprincipal unitary transport.',
 use='Local symbol-calculus input after phase-space localization. Global bounded derivative conditions are not asserted for the original noncompact curved model.'),
 dict(title='Semiclassical propagation of coherent states with spin-orbit interaction',
 url='https://arxiv.org/abs/math-ph/0403030',
 checked='Abstract distinguishes fixed spin and hbar*s fixed limits, with different effective dynamics.',
 use='Supports auditing the choice of scaling; not imported as an exact theorem for the entire Gauss graph.')],
 new='An explicitly different relative-scaling branch of the original finite graph, coupled to its Gauss packet and real scalar records. Exact Duhamel extracts a first source-response coefficient that the leading state error alone cannot resolve.',
 excluded='No universal scaling axiom, no old no-go reversal, no uniform space refinement, no epsilon=1 accuracy or full Einstein claim.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_761.md','joint_loop_scale_transport.py','joint_loop_scale_transport_results.json','unified_physics_condition_ledger_761.md')
write('round761_drafts/final_review.txt',
 'Primary-agent review only. A new asymptotic family is prominently declared; the old B-principal family and its no-go remain intact. Fixed compact phase-time tubes allow scalar-principal Egorov and the Gauss packet symbol lemma, with finite graph dependent constants. Residual epsilon^(3/2) yields finite-time norm error sqrt(epsilon), including specified source graph weights. First response is obtained independently via exact Duhamel and a commutator expansion through epsilon^2, so no invalid division of a sqrt(epsilon) error occurs. Finite scalar Kraus histories telescope to the same retarded source coefficient. Exact global Gauss expectations and gauge-frame components are distinguished. The nonzero original sterile-pair witness has unit-weight coefficient -0.0471463826 and requires the actual node weight for graph normalization. The oscillator calculation is explicitly only a frozen local coefficient calibration. No autonomous source measurement or apparatus, continuum/renormalization bridge, or full metric dynamics is claimed.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')

mapping={'761':'762','760':'761','759':'760','3467':'3470','3464':'3467',
         '1592':'1595','3627':'3637','3613':'3627','406':'407','311':'312',
         'joint_conditional_background':'joint_loop_scale_transport'}
pub=remap((HERE/'publish_round760.py').read_text('utf8'),mapping)
summary='**第761轮完成（新增缩放的条件结果）：** [有限时间输运与同源反作用]({p}research_note_761.md)在明示Hb,epsilon＋epsilon B分支中，完整矩阵、首阶背景响应及原字段报告有共同固定时间连接。旧主阶B族保留，未证epsilon=1或连续GR。三组、十五式通过，最新761／3470，1595份编号科学文件、3637份保护证据。[核验]({p}research_round_761_checks.json)、[条件账]({p}unified_physics_condition_ledger_761.md)。'
order='**当前执行顺序（761后，优先于下方历史安排）：** 接[762阶次依据与共同尺度]({p}round762_drafts/STATUS.md)，核相对能源阶次的资源解释及同一态/连续来源连接；不继续局部振子精度。[范围审计]({p}round761_drafts/scope_and_dedup_review.md)。旧空间、604、649／699、原缩放及总目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace("assert '**第761轮完成：**' not in text","assert '**第761轮完成（新增缩放的条件结果）：**' not in text")
pub=pub.replace('混合Gauss过程与条件来源流','有限时间输运与同源反作用').replace('条件相干与来源','新增阶次及同源响应')
write('publish_round761.py',pub)
post=remap((HERE/'postcheck_round760.py').read_text('utf8'),mapping)
post=post.replace("assert '**第761轮完成：**' in content","assert '**第761轮完成（新增缩放的条件结果）：**' in content")
write('postcheck_round761.py',post)
ver=remap((HERE/'verify_round760.py').read_text('utf8'),mapping)
start=ver.index('    entry=core.read(');end=ver.index('    result=model.run();',start)
ver=ver[:start]+ver[end:]
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_761.md','round761_drafts/research_note_761_draft.md',
 'round761_drafts/final_review.txt','round761_drafts/literature_scope_audit.json',
 'round761_drafts/scope_and_dedup_review.md','round762_drafts/STATUS.md',
 'round761_drafts/coarse_process_reuse_review.md')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace("checks['display_formulas']==16","checks['display_formulas']==15")
ver=ver.replace('exact_original_conditional_amplitude_representation_checked=True,local_density_matrix_current_summary_insufficient_on_stated_input_family=True,reduced_complexity_macroscopic_limit_proven=False,full_continuum_process_equivalence_proven=False',
 'new_relative_scaling_explicitly_registered=True,fixed_graph_transport_and_first_response_checked=True,old_scaling_not_refuted_or_replaced=True,epsilon_one_accuracy_proven=False,full_continuum_process_equivalence_proven=False')
write('verify_round761.py',ver)
print('Prepared761 explicitly conditional scaling branch and762 entry.')

