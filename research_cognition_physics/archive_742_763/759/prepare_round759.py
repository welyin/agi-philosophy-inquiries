"""Prepare759: a scoped dynamic closure test, with all prior evidence retained."""
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

ledger=(HERE/'unified_physics_condition_ledger_758.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：759确定物质背景的联合来源限制',1)
ledger=ledger.replace(ledger.split('\n')[2],
 '2026-10-04。接[758全账](unified_physics_condition_ledger_758.md)，回填[759报告](research_note_759.md)。[结果](joint_background_source_closure_results.json)、[核验](research_round_759_checks.json)。统一目标保持。',1)
updates={
 'C01':'759原最大混合稳定子纤维经758等距嵌入为正常严格Gauss态；不是外给不合法密度',
 'C03':'759有限菜单含动量/来源交叉矩；并未构造其读出装置或排除仅保有界字段报告的较弱近似',
 'C04':'759原完整演化的积分身份迫使单确定背景合同满足零源协方差；原合法输入严格违反',
 'C19':'759统一确定物质/规范参考轨道不能对758全部输入保持所列任务；外给空间度规未量子化',
 'C20':'759只排除固定T、统一联合矩及初始连续接入的确定轨道极限；保分布、记忆、初始层和量子慢变量',
 'C21':'759矩阵资料丰富而背景固定为一点仍可能丢联合任务，须一起保留必要关联',
 'C22':'759原质量与全部跳跃共用净系数kappa=Tr_K M²/d>0，质量源与背景动量不能按受检验合同分离'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'''

## 759原过程对确定背景合同的净限制

- H1/H2/H3没有新增公理；本轮检验其一个明确实现分支：确定物质/规范相轨道加连续内部量子矩阵，并在固定物理时间统一保有限动量/来源联合矩。
- 原H的两个积分矩身份推出必要源协方差为零。758合法最大混合稳定子纤维上，原质量/跳跃交叉迹为零且质量平方迹严格正，故此统一合同被排除。
- 结论不是由初始Taylor项交换极限所得，也不是把726/727的非零来源方差重报。它排除一个动态映射，保留原正量子过程和多背景/记忆/量子慢变量候选。
- 见证属于758全部允许态族；未证明757特定两态的完整动态极限。所列交叉矩并未有免费内部测量实现；仅保有界字段报告的更弱合同不被排除。
- 背景指物质与规范构形，外给空间度规仍固定。未反证GR、未推出引力量子性、未完成共同Q/E或连续手征极限。

## 下一项：保必要背景分布与关联

接[760](round760_drafts/STATUS.md)。优先检验现有矩阵输运/有序过程的原模型适用条件，明确必要分布、相干及来源如何共同保留。不把多分支当已完成正解；新谱隙、时间平均、缩放或噪声处方逐项登记。避免继续扫描反例常数，旧空间与604、649/699各自范围保持。
'''
write('unified_physics_condition_ledger_759.md',ledger)
write('round759_drafts/research_note_759_draft.md',(HERE/'research_note_759.md').read_text('utf8'))
write('round760_drafts/STATUS.md','''# 第760轮入口：共同背景分布、量子关联与实际过程

接[759限定反例](../research_note_759.md)、[759全账](../unified_physics_condition_ledger_759.md)、[共同候选](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 758初始Gauss矩阵桥保留；759排除对全部允许输入的单确定物质/规范相轨道合同：固定T、联合矩统一逼近、初始矩阵资料连续。没有排除原量子过程或全部宏观描述。
2. 下一项转向必要背景分布/相干的共同输送，不继续优化反例常数。分布或分支须从同一原过程得到，不能另外拟合噪声。
3. 优先回查525、555—558、643、704/706、726—730、741及756—758。成熟矩阵Egorov需实际谱重数/分离及观测类，不能以静态谱点或平坦分支替代原完整背景条件。
4. 保同一原未知输入、实际字段记录、完整后态与来源。数学谱分解不是免费可实施的能源测量；正单时分布也不等于正确多时记录。
5. 若采用时间粗化、初始层、有限误差或继续保量子慢变量，明确其范围和输入。不要为容易证明而更换原相互作用，或只给一个新两能级玩具却宣称完成原过程。
6. 原度规gamma仍是固定输入；物质目标的背景分支不是时空维数或量子引力证明。维数、群、作用、连续手征及几何反作用继续在共同主账内。
7. 不改目标、不建任务/调度、不做图像；空间382—386、425、522—523及604、649/699原边界保持。
''')
write('round759_drafts/scope_and_dedup_review.md','''# 759范围与去重审查

- 前一目标轮属于进展：758已正式发布，759合同回查已有保存。当前读取三导航、758笔记与结果及入口；无正在运行的Python，未重启旧计算。
- 350非仿射性、525冻结控制、726/727源涨落、728字符、729曲几何谱、730背景场都直接复用。本轮新增为原全图积分矩必要条件和同一合法输入对确定相轨道合同的违反。
- 不从758的t/sqrt(hbar)坏上界推出不可能；不把757初始四阶系数直接取hbar极限。固定时间反证依赖明示的有限菜单统一收敛和初始连续接入。
- 见证是758允许的最大混合稳定子纤维，不是757特定tau族。原源菜单含动量/来源交叉量，若只要原有界L概率，应另验较弱合同。
- 使用紧测试cutoff，未修改H。变量x是内部质量/场坐标，背景不是时空度规。反例不否定GR、全部半经典方法或统一目标。
- 全周期图适用性由解析CAR与稳定子证明；96模式是系数校准，不冒充全周期传播。所有原质量、Weyl方向、非零边和规范都保留。
- 后继转向原过程的分布/相干及来源输送，未宣布多背景候选已被证明。
''')
write('round759_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(title='A semiclassical Egorov theorem and quantum ergodicity for matrix valued operators',
 url='https://arxiv.org/html/math-ph/0204018',
 checked='H0 in section2: constant eigenvalue multiplicities and separation; theorem3.2: symbol regularity and block-diagonal observables.',
 use='A candidate transport route, not applied as a theorem for this full Gauss model.')],
 new='Exact integral joint-moment condition and original-center-sector trace witness against a stated deterministic material-background contract.',
 excluded='No unverified multiband limit, no universal no-go, no quantum spacetime claim, no endpoint derivative/semiclassical limit interchange.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_759.md','joint_background_source_closure.py','joint_background_source_closure_results.json','unified_physics_condition_ledger_759.md')
write('round759_drafts/final_review.txt',
 'Primary-agent review only. The localized C=chi B commutes exactly with B, avoiding an unbounded 1/hbar matrix commutator. The two actual integral moment identities pass to the explicitly required uniform finite-time menu limit. A deterministic phase point kills joint centered momentum/matrix terms, forcing zero symmetric covariance. The original maximal mixed stabilizer fiber has zero means and mass/hopping trace, with strictly positive mass square. Continuity at initial time gives a contradiction without exchanging Taylor limits. Witness is in the broad758 family, not the earlier special tau inputs. Material/gauge configuration is not spatial geometry. Multibranch, noisy, initial-layer, finite-tolerance and quantum-background alternatives stay open. No autonomous instruments or full evolution simulated.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'759':'760','758':'759','757':'758','3461':'3464','3458':'3461',
 '1586':'1589','3600':'3613','3589':'3600','404':'405','309':'310',
 'joint_gauss_matrix_background':'joint_background_source_closure'}
pub=remap((HERE/'publish_round758.py').read_text('utf8'),mapping)
summary='**第759轮完成（限定反例）：** [确定背景动态闭合的联合来源检验]({p}research_note_759.md)原合法混合态的质量源协方差严格正，排除对全部输入统一保指定联合矩的单确定物质相轨道合同；保多背景、记忆及量子慢变量。三组、十六式通过，最新759／3464，1589份编号科学文件、3613份保护证据。[核验]({p}research_round_759_checks.json)、[条件账]({p}unified_physics_condition_ledger_759.md)。不是对GR或全部宏观近似的反证。'
order='**当前执行顺序（759后，优先于下方历史安排）：** 接[760背景分布与共同过程]({p}round760_drafts/STATUS.md)，让原量子过程决定必要背景分布/相干及来源，先核成熟输运工具的真实前提；不继续反例常数或单项精度。[范围审计]({p}round759_drafts/scope_and_dedup_review.md)。认知共同候选、旧空间、604、649／699及目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('同一Gauss背景中的量子矩阵与初始记录','确定背景动态闭合的联合来源检验').replace('旧空间合同保持，量子矩阵初始桥','旧空间合同保持，确定物质背景限制')
# The summary now has a scope qualifier, and all guards must use that exact header.
pub=pub.replace("assert '**第759轮完成：**' not in text","assert '**第759轮完成（限定反例）：**' not in text")
write('publish_round759.py',pub)
post=remap((HERE/'postcheck_round758.py').read_text('utf8'),mapping)
post=post.replace("assert '**第759轮完成：**' in content","assert '**第759轮完成（限定反例）：**' in content")
write('postcheck_round759.py',post)
ver=remap((HERE/'verify_round758.py').read_text('utf8'),mapping)
start=ver.index('    import sys\n');end=ver.index('    result=model.run();',start)
ver=ver[:start]+'''    entry=core.read(HERE/'round759_drafts/dynamic_contract_entry_checks.json')
    assert core.digest(HERE/'round759_drafts/dynamic_contract_working_review.md')==entry['review_sha256']
'''+ver[end:]
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_759.md','round759_drafts/research_note_759_draft.md',
 'round759_drafts/final_review.txt','round759_drafts/literature_scope_audit.json',
 'round759_drafts/scope_and_dedup_review.md','round760_drafts/STATUS.md',
 'round759_drafts/dynamic_contract_working_review.md','round759_drafts/publish_dynamic_contract_entry.py',
 'round759_drafts/dynamic_contract_entry_checks.json','round759_drafts/build_check_log.md')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('original_Gauss_matrix_initial_bridge_checked=True,initial_source_record_and_spectral_scope_checked=True,finite_physical_time_equivalence_proven=False,full_continuum_process_equivalence_proven=False',
 'original_joint_integral_moment_identity_checked=True,stated_deterministic_material_background_contract_refuted=True,full_physical_theory_refuted=False,positive_dynamic_replacement_proven=False')
write('verify_round759.py',ver)
print('Prepared759 scoped dynamic closure result and760 entry.')

