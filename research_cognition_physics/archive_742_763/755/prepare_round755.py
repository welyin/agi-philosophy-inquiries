"""Prepare755 state qualification and joint task/source scope; exclusive writes."""
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


ledger=(HERE/'unified_physics_condition_ledger_754.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：755共同物理态、准自由来源与实际记录',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[754全账](unified_physics_condition_ledger_754.md)，回填[755报告](research_note_755.md)。[结果](joint_gaussian_physical_state_bridge_results.json)、[核验](research_round_755_checks.json)。联合目标保持。',1)
updates={
    'C01':'755闭整体全部夸克边缘的准自由忠实替换受中心支撑限制；只有中性纯Slater符号满足中心必要条件，不代表充分Gauss',
    'C03':'755原完整32模式singlet混合与准自由态同二点来源而两sterile实际联合记录不同；不能以来源匹配替代任务匹配',
    'C11':'755处理原内部规范Gauss物理态资格；平均约束及753/754经典补偿不等于量子支撑，未解引力量子约束',
    'C14':'755直接复用710全局颜色中心；有限中心保留不受753连续稳定子消失影响；开放区域须保边界载体',
    'C19':'755原有限图正温物理Gibbs的全部夸克边缘非准自由；相同协方差的Gaussian替代不保物理态，未计算实际热协方差',
    'C20':'755排除完整夸克态不变且只用玻色补偿的特定桥；背景Hadamard与有限图的其他关系/有效映射仍开放',
    'C22':'755原Gauss物理关联态可保所有同玻色态上的瞬时CAR二次来源，但Gauss平方及联合记录不同；无同753几何或全动态来源完成'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n\n## 755共同物理态与同源任务资格\n\n'
ledger+='''- H1—H3共同选择“物理态、来源与继续任务同时匹配”的检验；观察与替代解释沿[未编号工作报告](round755_drafts/cognitive_joint_candidate_working_report.md)，不把认知类比当Gauss定理。
- 710入口的全局颜色中心直接复用。原有限闭整体全部夸克的准自由边缘满足零中心支撑，当且仅当协方差是秩为3倍数的投影；这只是中心充分，非全Gauss充分。
- 原物理Gibbs忠实性使固定有限图正温夸克协方差严格介于0与I，因此真实夸克边缘不等于其准自由替代。未数值计算热谱或跨尺度一致界。
- 754的两模式空/满比较若被解释为原有限CAR中其余夸克完全不变的忠实桥，则至少一份态与任何物理完成的夸克边缘距离不小于1/2；玻色补偿不能改变全局中心。该有限桥限制不否定连续背景Hadamard构造。
- 原完整32模式的全空/全满均为完整局域规范singlet，混合态与Gamma(pI)同所有二点函数及瞬时二次来源。p=1/2时联合sterile记录距离1/2、颜色Gauss平方差8。
- 因而纯来源匹配有真实正面实例，却不是完整状态/过程匹配。替代候选保Gauss关联态及实际记录，准自由场仅在已证有效范围采用。
- 不采用710另加四费米作用，不改原Hamiltonian。未证明自治准备、同753量子几何、连续映射或全量子引力；空间、604、649/699保持。

## 本轮合并与下一项

本轮排除一个完整态忠实接法，保留并界定只匹配来源的竞争接法；没有减少群、维数、作用或参数输入。接[756](round756_drafts/STATUS.md)，围绕同一物理关联态与完整记录向连续来源的映射；不再展开中心标签枚举或局部Gaussian扫描。应用目标保持。
'''
write('unified_physics_condition_ledger_755.md',ledger)
write('round755_drafts/research_note_755_draft.md',(HERE/'research_note_755.md').read_text('utf8'))
write('round756_drafts/STATUS.md','''# 第756轮入口：同一物理关联态、实际记录与连续来源

接[755](../research_note_755.md)、[全条件账](../unified_physics_condition_ledger_755.md)、[共同候选工作报告](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 755排除的是原有限闭图全部夸克态忠实Gaussian替换，不否定局部边缘、连续Hadamard或有效来源描述。正面非Gaussian全Gauss态能保瞬时来源，却改变实际联合记录。
2. K0因此保同一物理关联态和全部实际记录；接回633/730非Gaussian后态、743原相互作用、724实际区域及704/706过程与来源保留，先核可共同输送的原对象。
3. 下一实质单元须给同一准备、相互作用、记录后态、全来源与几何接口的明确对应；若只能比较均值，必须声明任务和范围，不能再说同一量子过程已完成。
4. 不把中心投影即当完整局域Gauss，更不当引力量子态。有限图到连续手征场、物理内积及有效误差仍是核心缺口，604/649/699保持其限定范围。
5. 不继续中心标签、Gaussian谱或小矩阵精度枚举，不增加免费探针。仅成熟工具整理或旧结果推论存工作报告，不计科学轮次。
6. 空间382—386、425、522—523复用，不恢复Lipschitz。应用目标不改，不建任务、定时或图像检查。
''')
write('round755_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
        dict(title='Fermionic Quasi-free States and Maps in Information Theory',url='https://arxiv.org/html/0709.1061',
             checked='Section IV finite CAR density spectrum; paper explicitly treats finite dimensions.',
             use='Mature independent-mode formula; original Gauss centre support and pair mapping checked separately.'),
        dict(title='Reference frames, superselection rules, and quantum information',url='https://arxiv.org/html/quant-ph/0610030',
             checked='Group averaging, sectors, relational encoding.',
             use='Alternative modelling framework; statistical invariance is not strict Gauss support.')],
    inherited='709/710 closed centre;728 stabilizer characters;730 background scope;742 Wick limitation;754 classical source completion.',
    new='Finite original Gaussian-state qualification; fixed-complement pair distance; original complete local singlet source/record comparison.',
    excluded='Infinite continuum number constraint, replacing open regions with closed ones, quantum Einstein solution or autonomous measurement.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_755.md','joint_gaussian_physical_state_bridge.py','joint_gaussian_physical_state_bridge_results.json','unified_physics_condition_ledger_755.md')
write('round755_drafts/final_review.txt',
      'Primary-agent review only. Mature centre rule reused. All-quark finite closed-system marginal, not local/open region or unconstructed continuum number operator. Pure-neutral Gaussian classification is a centre test only. Gibbs covariance not numerically supplied. Full32 vacuum/filled states are local gauge singlets; same bilinears do not preserve joint records; original H does not preserve the two-vector code. No new action, actual753 quantum geometry, autonomy or full unification claim.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'755':'756','754':'755','753':'754','3449':'3452','3446':'3449',
         '1574':'1577','3550':'3562','3541':'3550','400':'401','305':'306',
         'joint_irreducible_source_completion':'joint_gaussian_physical_state_bridge'}
pub=remap((HERE/'publish_round754.py').read_text('utf8'),mapping)
summary='**第755轮完成：** [共同物理态、来源与联合记录]({p}research_note_755.md)原闭Gauss全部夸克态的忠实准自由替换受严格支撑限制；两模式填充而其余态不变的接法，至少一份边缘误差≥1/2。原完整singlet态可保瞬时二次来源，却改变联合记录；未完成连续或量子引力桥。三组、十六式通过，最新755／3452，1577份编号科学文件、3562份保护证据。[核验]({p}research_round_755_checks.json)、[条件账]({p}unified_physics_condition_ledger_755.md)。'
order='**当前执行顺序（755后，优先于下方历史安排）：** 按[认知观察与共同候选]({p}round755_drafts/cognitive_joint_candidate_working_report.md)，接[756同一关联态与连续来源]({p}round756_drafts/STATUS.md)，保原Gauss态、实际记录和全部来源，不以Gaussian摘要替换整个过程；不继续中心标签或谱精度支线。[范围审计]({p}round755_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及应用目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('全颜色来源与共同约束','共同物理态、来源与联合记录').replace('旧空间合同保持，量子来源与完整初始约束','旧空间合同保持，Gauss物理态与来源任务')
write('publish_round755.py',pub)
write('postcheck_round755.py',remap((HERE/'postcheck_round754.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round754.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_755.md','round755_drafts/research_note_755_draft.md',
       'round755_drafts/final_review.txt','round755_drafts/literature_scope_audit.json',
       'round755_drafts/scope_and_dedup_review.md','round756_drafts/STATUS.md',
       'round755_drafts/cognitive_joint_candidate_working_report.md',
       'round755_drafts/cognitive_joint_candidate_working_report_checks.json',
       'verify_cognitive_joint_candidate_review.py')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace("checks['display_formulas']==20","checks['display_formulas']==16")
ver=ver.replace('all_mode_color_inverse_and_joint_source_scope_checked=True,nonlinear_prescribed_initial_constraints_solved=True,self_consistent_quantum_feedback_proven=False',
                'finite_closed_Gauss_state_bridge_scope_checked=True,source_matching_not_claimed_as_process_matching=True,full_continuum_physical_state_constructed=False')
ver=ver.replace("    result=model.run();", "    import verify_cognitive_joint_candidate_review as working\n    assert working.verify()==core.read(working.TARGET)\n    result=model.run();")
write('verify_round755.py',ver)
print('Prepared755: state qualification and joint source/record comparison.')
