"""Prepare756 common correlated Hadamard/source connection; exclusive writes."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_755.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：756同一关联态、连续记录与来源涨落',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[755全账](unified_physics_condition_ledger_755.md)，回填[756报告](research_note_756.md)。[结果](joint_correlated_hadamard_noise_results.json)、[核验](research_round_756_checks.json)。联合目标保持。',1)
updates={
 'C01':'756原中性两模式正关联族具有同二点Hadamard背景扩展；不是全有限图Gauss态连续映射',
 'C03':'756同一p/q指定原中性占据联合概率及完整来源噪声差；未构造自主局域记录装置',
 'C11':'756光滑噪声差满足完整联合Ward，可继承732/754线性补偿与发展；不是全量子约束态',
 'C19':'756明确保共同Hadamard参考、模式及互补态输入；非热平衡态不由对易响应唯一决定噪声',
 'C20':'756给有限中性记录代数到连续背景态的正扩展，同一关联输送噪声；未等同全部过程/尺度',
 'C21':'756若保该联合任务和来源二阶统计，完整二点摘要仍须至少保q；未证明任意任务有限摘要',
 'C22':'756同全部均值/对易响应但不同联合噪声；同q与原Ys固定差核及条件联合响应，不能独立选择噪声'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n\n## 756同一关联态的连续共同实现\n\n'
ledger+='''- 沿H1—H3，原中性双模式τ(p,q)与同一互补态组成正Hadamard背景态族，全部二点与二次来源均值相同；这是一份指定代数的连接，不抹掉755全部夸克支撑限制。
- q−p²是同一联合记录的关联，也是完整来源噪声差的唯一可变态参数。原Majorana耦合使标量力差核严格非零，无需新耦合或独立噪声源。
- 差核是光滑有限模式乘积，保T/j/s完整Ward；不可把差核视为独立正环境。绝对噪声需要时空涂抹，初片积分只对光滑差做。
- 原线性联合K逐光滑系数输送差核，至少一项完整联合场协方差差非零。未证明纯度规分量非零，未数值求量子Einstein协方差。
- 记录相关与来源涨落不再是独立自由输入；维数、群、作用、参考/制备及有限图连续映射仍未推出。无全量子引力、自主装置或ε=1反馈结论。

## 本轮合并与下一项

一个正态族共同承担指定记录与同一背景上的完整来源统计，关闭该局部对象接口；并排除只存二点函数的任务摘要。接[757](round757_drafts/STATUS.md)，核实际instrument的条件后态及来源响应是否保持同一过程，不继续噪声谱/网格精度。旧空间382—386、425、522—523及604/649/699保持，应用目标不改。
'''
write('unified_physics_condition_ledger_756.md',ledger)
write('round756_drafts/research_note_756_draft.md',(HERE/'research_note_756.md').read_text('utf8'))
write('round757_drafts/STATUS.md','''# 第757轮入口：同一记录的条件后态与来源共同输送

接[756](../research_note_756.md)、[条件账](../unified_physics_condition_ledger_756.md)、[认知候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 756给同一有限中性态族、背景Hadamard扩展和完整来源噪声的明确连接；q固定记录与噪声，背景/制备仍输入。尚非完整原Gauss图过程的连续映射。
2. 下一项先核原中性占据instrument的条件后态、全部来源以及联合约束响应：能否用同一正过程同时输送，分支概率与未选择态必须一致。复用633/634/730和750，不能逐分支重选准备。
3. 若仅重述全方差公式或普通CP后态，存工作报告不编号。新增验收须给此原共同对象的实际过程映射，或严格排除某个明确接法。
4. 不把全空间分布模仪器当瞬时局域装置；内部制备与传播仍须交代。若必须回到743原量子标量相互作用，说明它与背景二次分支的连接尚缺，不假称已经关闭。
5. 不继续单独噪声谱、网格精度或中心标签，不增免费探针。H1—H3、物理输入与量子/有效边界分开记录。
6. 空间382—386、425、522—523及604/649/699范围保持；不改目标、不建应用任务/调度、不做图像。
''')
write('round756_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[
  dict(title='Stochastic Gravity: Theory and Applications',url='https://arxiv.org/html/0802.0658',
       use='Mature noise and induced response framework, inherited620/634; joint T/j/s Ward checked separately.'),
  dict(title='Fermionic Quasi-free States and Maps in Information Theory',url='https://arxiv.org/html/0709.1061',
       use='Finite CAR positivity and mode factorisation only; continuum regularity established by finite smooth perturbation.'),
  dict(title='The locally covariant Dirac field',url='https://arxiv.org/abs/0911.1304',
       use='Previously checked730 smooth local operations; current positive extension explicitly constructed.')],
 inherited='620/634 noise;730 Hadamard;732 joint Ward/response;735 subtraction;742 Gaussian boundary;755 state qualification.',
 new='Same finite neutral correlation fixes actual joint records and smooth complete-source noise on original753 background.',
 excluded='Full Gauss continuum lift, autonomous local apparatus, scalar quantum H equivalence, quantum Einstein covariance or unique physical theory.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_756.md','joint_correlated_hadamard_noise.py','joint_correlated_hadamard_noise_results.json','unified_physics_condition_ledger_756.md')
write('round756_drafts/final_review.txt',
 'Primary-agent review only. Positive even two-mode states, not necessarily even-sector pure support. Smooth finite-rank Hadamard construction preserves selfdual reality; entire two-point unchanged. One q fixes full joint noise difference, not independent positive noise. Full T/j/s Ward required. Spatial integration only of smooth difference; energy numerical component is mass pairing only. Conditional K maps smooth finite coefficients, not full UV noise or quantum gravity. At least one joint component changes, metric-only nonzero not proved. Original full Gauss lift and autonomous preparation remain open.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'756':'757','755':'756','754':'755','3452':'3455','3449':'3452',
 '1577':'1580','3562':'3574','3550':'3562','401':'402','306':'307',
 'joint_gaussian_physical_state_bridge':'joint_correlated_hadamard_noise'}
pub=remap((HERE/'publish_round755.py').read_text('utf8'),mapping)
summary='**第756轮完成：** [同一关联态、连续记录与来源涨落]({p}research_note_756.md)原中性正关联族接到753同背景Hadamard态：同二点及平均来源，联合记录q却决定完整光滑噪声差与条件联合响应。未完成全Gauss连续映射或量子引力。三组、十八式通过，最新756／3455，1580份编号科学文件、3574份保护证据。[核验]({p}research_round_756_checks.json)、[条件账]({p}unified_physics_condition_ledger_756.md)。'
order='**当前执行顺序（756后，优先于下方历史安排）：** 沿[认知共同候选]({p}round755_drafts/cognitive_joint_candidate_working_report.md)，接[757同一条件后态与来源]({p}round757_drafts/STATUS.md)，保原记录概率、完整来源与联合约束，不另选独立噪声；不继续谱或网格精度支线。[范围审计]({p}round756_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及应用目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('共同物理态、来源与联合记录','同一关联态、连续记录与来源涨落').replace('旧空间合同保持，Gauss物理态与来源任务','旧空间合同保持，同一记录与来源统计')
write('publish_round756.py',pub)
write('postcheck_round756.py',remap((HERE/'postcheck_round755.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round755.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_756.md','round756_drafts/research_note_756_draft.md',
 'round756_drafts/final_review.txt','round756_drafts/literature_scope_audit.json',
 'round756_drafts/scope_and_dedup_review.md','round757_drafts/STATUS.md',
 'round756_drafts/initial_source_before_tuple_fix.txt',
 'round756_drafts/source_before_component_label_review.txt',
 'round756_drafts/results_before_component_label_review.json')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace("checks['display_formulas']==16","checks['display_formulas']==18")
ver=ver.replace('finite_closed_Gauss_state_bridge_scope_checked=True,source_matching_not_claimed_as_process_matching=True,full_continuum_physical_state_constructed=False',
 'same_record_Hadamard_source_noise_connection_checked=True,conditional_joint_response_scope_checked=True,full_continuum_Gauss_state_constructed=False')
write('verify_round756.py',ver)
print('Prepared756: same correlated state, records and complete-source noise.')
