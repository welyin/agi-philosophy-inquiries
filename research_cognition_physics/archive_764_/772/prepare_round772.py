"""Freeze 772 evidence and prepare guarded publication; preserve old files."""
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def write(name, content):
    path=HERE/name; path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:
        f.write(content)


lines=(HERE/'unified_physics_condition_ledger_771.md').read_text('utf8').splitlines()
lines[0]='# 联合条件总账：772同一物理态的自由BRST配对与来源字典'
lines[2]='2026-10-04。接[771全账](unified_physics_condition_ledger_771.md)，回填[772报告](research_note_772.md)。[结果](joint_brst_quartet_matching_results.json)、[核验](research_round_772_checks.json)。同一物理态可有适配自由BRST表示；来源代表须与二阶均值共同变换。相互作用荷、局部插入与实际仪器仍开放。'
updates={
    'C01':'772保767物理态和原ghost极化，改选非物理扩展；没有重选物理真空或全部保留旧未约束核',
    'C11':'772平滑代表差使来源改变−P m；均值及其初值同步加m，原首阶约束响应保持',
    'C19':'772自由有限粒子BRST具有显式收缩及正物理商；相互作用荷和完整局部复合域未证',
    'C22':'772已定义规范不变观测的首阶系数在共同代表变换下相同；不等同于实际完整instrument或有限强度过程',
}
seen=set()
for i,line in enumerate(lines):
    for key,value in updates.items():
        if line.startswith('|'+key+' '):
            parts=line.split('|'); parts[2]+='；'+value; lines[i]='|'.join(parts); seen.add(key)
assert seen==set(updates)
write('unified_physics_condition_ledger_772.md','\n'.join(lines)+'''

## 772最新结论（优先于上述保留的771历史状态）

766—767原完整Cauchy规范像可配一份频率相容的各向同性对偶方向。原物理协方差保留，配对的非物理扩展与原核只差平滑项，保Hadamard、CCR、现实及BRST身份。在自由代数有限粒子域，配对数算符给显式收缩，closed向量的范数等于物理分量范数，零范数即exact。它未构造相互作用BRST荷、所有局部Wick乘积的共同闭域或有限耦合量子态。

扩展改变不是物理态重新准备。原ghost保持使差核两腿都为物理解；其物理商压缩为零，故可因子化为规范—物理平滑核。完整Noether身份给Delta J=−S''m；将均值、规范修正及初始资料一同输送，Pw+j不变，已有规范不变观测的首阶系数亦不变。不能只改来源而冻结二阶均值，也不能把所有未约束二点核都说成未改。

认知动机、原作用、四维、群/表示/参数及F>0保持。旧732/754响应直接应用不重计轮次，空间382—386/425/522—523、604和649/699保原范围。

## 下一项：773原有限阶相互作用插入与物理态/观测提升

用772适配自由表示和同一初始字典，核完整原三、四阶顶点所需的因果Wick插入及BRST身份。771一次来源处方要与新插入兼容，不能把自由正性当相互作用QME。保原记录、条件概率与资源账；不再优化自由配对矩阵，不先要求全阶或有限强度解。
''')
write('round772_drafts/research_note_772_draft.md',(HERE/'research_note_772.md').read_text('utf8'))
write('round773_drafts/STATUS.md','''# 第773轮入口：原相互作用的有限阶共同插入

接[772](../research_note_772.md)及[772条件账](../unified_physics_condition_ledger_772.md)。

1. 772已保原物理态构造自由配对表示，并证明来源与二阶均值/初值字典补偿；不把未约束旧核也说成未变。
2. 自由代数有限粒子域的正性不是全部局部复合域或相互作用幂零荷。接9807078及Fröb1803.10235的准确前提；不靠名称宣布全相互作用无反常。
3. 原有限阶需要三、四价顶点及局部反项，完整背景和低阶耦合留在传播子中。先核当前所需插入，不要求先证明全阶或有限耦合收敛。
4. 原771一次来源规范化保留；若共同插入迫使重新选有限项或代表，必须明示改动，并保同一态、观测及初值字典。只算平均来源不是原实际记录过程。
5. 直接回用732/754、734/735、623—625/704/756、768/769。空间382—386/425/522—523与604/649/699不重证、不扩大。四维、群、作用与参数仍输入。
6. 无新任务、调度、图像或目标改动。全相互作用QME、实际仪器、原Q连续映射、有限量子强度反馈仍开放。
''')
write('round772_drafts/scope_and_dedup_review.md','''# 772范围、对象与去重复核

- 上次目标轮完成771并保存772入口，为有效进展。本次核文件与进程，无运行Python；未改目标或使用代理。
- 原M逆、Hadamard物理态、Noether二次变分直接复用。新增是同态的适配扩展、自由代数BRST表示及来源/均值字典共同实现。
- L的各向同性由半个K项补成；K*QL=Q0而不是辅助正内积。全部伴随保原非正交资料基。
- 频率相容至平滑误差依赖旧全阶交织；不只用主符号。保原物理mu与ghost，不保证保旧未约束lambda。
- 物理准自由GNS允许混合态/表示倍增，空方向先商去。配对收缩保持频率及其radical，范数结论限代数有限粒子域。
- 不称相互作用荷已存在，不称所有局部Wick算符的共同闭域已完成，未直接签收9807078定理4的全部前提。
- 差核按算符指标与两个上指标分开转换；原ghost相同使物理字段差满足完整P两腿方程。因子化用原Cauchy左逆后再传播。
- 来源差为−P m，均值加m；必要时另解D0的规范修正。初值/齐次准备同步运输，不能分别重选成同一数值。
- 物理观测仅签收已定义的规范不变泛函之首阶表示相容，不签收完整多时instrument、概率正性或有限量子强度。
- 第一组是原完整主符号加有限平滑块类比；第二组为总次数不变的quartet子空间；第三组原势驻点不是753动态背景。数值均不承担无限维/全PDE证明。
- 无为重复旧反例或矩阵精度新增轮次；旧空间、604与649/699保原范围。主代理审查，无独立审稿。
''')
write('round772_drafts/literature_scope_audit.json',json.dumps(dict(sources=[
    dict(url='https://arxiv.org/html/1407.8079', location='Sections 2.7 and 3',
         use='Inherited linear BRST/physical phase space and two-point construction.',
         limitation='Does not by itself provide the interacting BRST charge on the original coupled model.'),
    dict(url='https://arxiv.org/html/hep-th/9807078', location='Section 4, equation (4.8) and Theorem 4',
         use='Target free positivity/null-exact condition and conditional formal physical-state deformation.',
         limitation='772 proves an adapted algebraic free representation; common local Wick domains and interacting nilpotent hermitian charge remain unproved.'),
], inherited_state_on_physical_algebra_preserved=True, unphysical_extension_changed=True,
    source_and_mean_dictionary_transported=True, algebraic_finite_particle_scope=True,
    interacting_QME_proven=False, full_Wick_domain_proven=False),ensure_ascii=False,indent=2))

mapping={'772':'773','771':'772','770':'771','3499':'3502','3496':'3499',
         '1625':'1628','1622':'1625','3773':'3789','3757':'3773',
         'joint_local_ward_repair':'joint_brst_quartet_matching','local_contact_entry':'first_response_entry'}
pattern=re.compile('|'.join(map(re.escape,sorted(mapping,key=len,reverse=True))))
convert=lambda s:pattern.sub(lambda m:mapping[m[0]],s)
verify=convert((HERE/'verify_round771.py').read_text('utf8'))
verify=verify.replace("checks['display_formulas'] == 19","checks['display_formulas'] == 16")
start,stop=verify.index('                original_same_free_state_preserved='),verify.index('                primary_code_and_note_review_completed=')
verify=verify[:start]+'''                original_physical_free_state_preserved=True,
                unphysical_extension_changed=True,
                algebraic_free_BRST_quartet_positivity_proven=True,
                source_and_mean_dictionary_matched=True,
                original_single_source_Ward_inherited=True,
                full_Wick_operator_domain_proven=False,
                interacting_QME_proven=False,
                finite_strength_self_consistency_proven=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
'''+verify[stop:]
write('verify_round772.py',verify)
publish=convert((HERE/'publish_round771.py').read_text('utf8'))
start,stop=publish.index('summary = '),publish.index('planned = ')
summary='**第772轮完成（同态BRST配对与来源字典）：** [自由物理表示与共同均值](%PRE%research_note_772.md)保原物理态，改选非物理配对扩展；自由有限粒子物理商为正，来源代表差由二阶均值及初值共同补偿。相互作用荷与实际仪器仍未证。三组、十六式通过，最新772／3502，1628份编号科学文件、3789份保护证据。[核验](%PRE%research_round_772_checks.json)、[条件账](%PRE%unified_physics_condition_ledger_772.md)。'.replace('%PRE%','{p}')
order='**当前执行顺序（772后，优先于下方历史安排）：** 接[773原有限阶共同插入](%PRE%round773_drafts/STATUS.md)，核因果Wick插入、BRST荷与物理态/记录字典；不继续自由配对矩阵或重做守恒来源响应。[范围审计](%PRE%round772_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'.replace('%PRE%','{p}')
publish=publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish=publish.replace('第772轮完成（完整一圈来源的局部Ward修复）','第772轮完成（同态BRST配对与来源字典）')
publish=publish.replace('## 417.','## 418.').replace('## 322.','## 323.')
publish=publish.replace('同态完整来源的局部Ward修复','同态自由BRST配对与来源字典')
publish=publish.replace('同态来源与联合守恒','自由物理表示与共同均值')
publish=publish.replace('旧空间合同保持，完整来源与共同响应','旧空间合同保持，同态表示与首阶均值')
write('publish_round772.py',publish)
post=convert((HERE/'postcheck_round771.py').read_text('utf8'))
post=post.replace('第772轮完成（完整一圈来源的局部Ward修复）','第772轮完成（同态BRST配对与来源字典）')
write('postcheck_round772.py',post)

names=['research_note_772.md','joint_brst_quartet_matching.py','joint_brst_quartet_matching_results.json','unified_physics_condition_ledger_772.md']
write('round772_drafts/final_review.txt','Primary-agent review; no independent reviewer.\n'
      +'Reviewed source equations, pairing, ghost signs, source/mean sign, scope and all three executed calibrations.\n'
      +'No claim of interacting QME, full Wick domain, actual records or finite-strength theory.\n'
      +'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared 772 evidence and guarded publication.')
