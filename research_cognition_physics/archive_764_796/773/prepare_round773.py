"""Freeze 773 and prepare guarded publication without altering past evidence."""
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
def write(name,content):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(content)

lines=(HERE/'unified_physics_condition_ledger_772.md').read_text('utf8').splitlines()
lines[0]='# 联合条件总账：773原物质参考的局部BRST修复工具'
lines[2]='2026-10-04。接[772全账](unified_physics_condition_ledger_772.md)，回填[773报告](research_note_773.md)。[结果](joint_material_local_brst_results.json)、[核验](research_round_773_checks.json)。正则原参考片给有限微分规范左逆；扩大背景jet系数类下有形式局部原始元。不是原全局UV类别或完整QME的完成。'
updates={
 'C01':'773的局部BV收缩与772的Fock表示分开；原物理准备保持，未把局部系数类扩张当全部相互作用态',
 'C11':'773可保已相容低阶反项再递推；771来源与真实反常首项的对象对应待核',
 'C19':'原Higgs、协变导数及颜色电磁曲率共同给有限微分左逆；仅753正则物质参考片，未新增参照物种',
 'C22':'允许逆参考Jacobian与点态Gram的形式局部类中，正ghost障碍可逐阶求原始元；原UV类别、跨片和实际插入仍开放',
}
seen=set()
for i,line in enumerate(lines):
    for key,value in updates.items():
        if line.startswith('|'+key+' '):
            parts=line.split('|');parts[2]+='；'+value;lines[i]='|'.join(parts);seen.add(key)
assert seen==set(updates)
write('unified_physics_condition_ledger_773.md','\n'.join(lines)+'''

## 773最新结论（优先于保留的772安排）

原753背景的满秩X提供微分同胚参数的局部左逆；H与D_xH组成非退化二列，原颜色全曲率又能辨认所有八个颜色方向。只用空间磁曲率有一维中心化子，已保留失败资料；补入原电曲率，不加新物种。时空/内部参数合起来有至多三阶微分左逆，vierbein冗余另有零阶左逆。

这是比“无连续背景稳定子”更强的局部结果。它使最小BV的字段、ghost和反场在原在壳背景处形成局部可收缩配对；系数容许逆参考Jacobian与逆点态Gram时，局部正ghost同调消失，并可逐涨落次数提升。每阶局部，不主张无穷和收敛或全式统一有限导数阶。

新增的是明确的补丁系数分支，不等同735严格UV多项式类别或771在更大背景类中的一次Wick处方。左逆在颜色幅度归零、参考退化时无统一控制；跨片、背景独立和完整相互作用QME均未签收。未往作用中实际加入这些反项，772同一物理态和原经典作用保持。

## 下一项：774真实有限阶插入与原来源首项

先核实际局部反常及771一次来源修复的低阶对应，再用773保首项递推检验共同处方。字段独立性、子插入、辅助/frame部门及同态初值字典须一起保留。若扩大类不符合原目标规范化条件，明确原分支与替代分支，不偷换成已经完成的全局理论。记录过程与原Q连续映射继续开放；旧空间、604、649/699及目标不变。
''')
write('round773_drafts/research_note_773_draft.md',(HERE/'research_note_773.md').read_text('utf8'))
write('round774_drafts/STATUS.md','''# 第774轮入口：局部修复分支与同一来源首项

接[773](../research_note_773.md)及[773全账](../unified_physics_condition_ledger_773.md)。

1. 773只在原753正则物质参考片构造有限微分规范左逆；原颜色磁场单独秩7，须保原电场。不要重复秩或配对矩阵实验。
2. 允许逆J与逆点态Gram的逐阶局部系数类是新分支；其正ghost原始元存在性不等于旧严格UV类别、全局拼接或实际QME已成立。
3. 优先核原有限阶作用反常的首项与771来源修复是否对应；如相容，可用773保首项递推。不得直接用772的Fock收缩作局部反项。
4. 原相互作用、反场/frame与ghost插入、字段独立性及子插入共同核。别把一份背景规范Ward等同全部量子BRST身份。
5. 保772实际物理准备、来源/均值/初值字典。局部方案变化可能改变物理有限项，必须明确比较，不能一律说无影响。
6. 旧空间382—386、425、522—523，604与649/699保持。无新任务、定时、图像或目标更改；实际记录、Q连续及有限强度反馈仍开放。
''')
write('round773_drafts/scope_and_dedup_review.md','''# 773范围与终审

- 先核导航、772与773工作报告及结果；上一轮为有效进展。详细进程查询被拒后使用普通进程列表，未见Python进程。无权限升级或新代理。
- 直接复用573/753参考秩、完整经典背景、765 PR=0、771来源和772物理态。新结果是有限微分左逆，而非再证无稳定子。
- 磁曲率单独秩7，初版失败保存；完整原电磁曲率恢复秩8。未增场、群、作用或物理耦合。
- L只反演点态矩阵，未用Green逆。计入X的度规变分、DH/F的连接变分及全部参数jet。
- 原规范参数和765协变参数通过点态字典转换；点态Gram的辅助伴随与密度正式转置不同。
- 含反场的局部拆分有约束，但互逆；约束理想在s0和收缩下保持。密度配对可检查其分部积分后的分解。
- 773局部jet收缩不是772代数Fock收缩。仅原在壳B0保证非线性经典微分提高涨落次数。
- 正ghost结论采用明确较宽的正则片系数类，允许逆参考Jacobians/Gram；不宣称严格多项式类别中消去原反常。
- 逐阶导数阶有限，无统一全阶导数界/级数收敛声明；全局参考退化和跨片仍待核。
- 真实量子反常未计算、共同QME/字段独立性未签收；保首项定理是条件结论，不能自动把771一次来源等同相互作用首项。
- 数值第一组为原初片和参数jets，第二组为通用离壳诊断，第三组为精确局部超多项式；没有把它们当原圈积分或完整动态过程。
- 主代理审查，无独立审稿或图像检查；目标和旧限定结论保持。
''')
write('round773_drafts/literature_scope_audit.json',json.dumps(dict(sources=[
 dict(url='https://arxiv.org/pdf/hep-th/0002245',location='Section 2.7, equations (2.48)-(2.51)',use='Contractible pairs including derivatives, antifields and local functionals modulo a total derivative.',limitation='The original finite differential left inverse and enlarged patch class are established separately here.'),
 dict(url='https://arxiv.org/pdf/1806.04695',location='Section 2.3, equations (2.21)-(2.23)',use='Comparison with ghost-exact compensator variables.',limitation='Its added Stueckelberg model does not prove original model existence, masses or global Abelianization.'),
 dict(url='https://arxiv.org/pdf/1803.10235',location='Theorems 3, 10-12',use='Local anomaly and conditional interacting-insertion framework.',limitation='No original loop anomaly or common time-ordered scheme is constructed in 773.'),
],extra_coefficient_class='Smooth finite background jets on the regular material patch, allowing inverse material Jacobian and pointwise Gram; formal fluctuation-degree completion.',
 global_polynomial_UV_theorem=False, original_QME_proven=False, free_physical_state_unchanged=True),ensure_ascii=False,indent=2))

mapping={'773':'774','772':'773','771':'772','3502':'3505','3499':'3502','1628':'1631','1625':'1628',
         '3789':'3806','3773':'3789','joint_brst_quartet_matching':'joint_material_local_brst',
         'first_response_entry':'local_insertion_entry'}
pattern=re.compile('|'.join(map(re.escape,sorted(mapping,key=len,reverse=True))))
convert=lambda s:pattern.sub(lambda m:mapping[m[0]],s)
verify=convert((HERE/'verify_round772.py').read_text('utf8'))
verify=verify.replace("checks['display_formulas'] == 16","checks['display_formulas'] == 14")
verify=verify.replace("             'round773_drafts/local_insertion_entry.py',\n",'')
verify=verify.replace("             'round773_drafts/local_insertion_entry_results.json',\n",'')
verify=verify.replace("'round773_drafts/entry_postpublication_checks.json')", "'round773_drafts/entry_postpublication_checks.json',\n             'round773_drafts/initial_magnetic_only_source.txt',\n             'round773_drafts/initial_magnetic_only_failure.json',\n             'round773_drafts/pre_review_results.json')")
start,stop=verify.index('                original_physical_free_state_preserved='),verify.index('                primary_code_and_note_review_completed=')
verify=verify[:start]+'''                original_physical_free_state_preserved=True,
                finite_differential_gauge_left_inverse_proven=True,
                local_positive_ghost_contraction_in_enlarged_patch_class=True,
                original_single_source_Ward_inherited=True,
                enlarged_coefficient_class_explicit=True,
                original_loop_anomaly_computed=False,
                interacting_QME_proven=False,
                strict_global_polynomial_counterterm_result=False,
                actual_instrument_or_Q_continuum_equivalence_proven=False,
'''+verify[stop:]
write('verify_round773.py',verify)
publish=convert((HERE/'publish_round772.py').read_text('utf8'))
start,stop=publish.index('summary = '),publish.index('planned = ')
summary='**第773轮完成（原物质参考与局部BRST）：** [局部规范左逆与形式修复]({p}research_note_773.md)原参考、Higgs及颜色电磁曲率给有限微分左逆；允许逆参考jet的正则片系数类中可逐阶求局部原始元。旧UV类别、共同QME和全局延伸未证。三组、十四式通过，最新773／3505，1631份编号科学文件、3806份保护证据。[核验]({p}research_round_773_checks.json)、[条件账]({p}unified_physics_condition_ledger_773.md)。'
order='**当前执行顺序（773后，优先于下方历史安排）：** 接[774局部反项与来源首项]({p}round774_drafts/STATUS.md)，核原有限阶插入和771处方是否能共同实现；不重复参考秩或自由配对。[范围审计]({p}round773_drafts/scope_and_dedup_review.md)。补丁与系数类新增条件明示，原物理态、旧空间及目标保持。'
publish=publish[:start]+'summary = '+repr(summary)+'\norder = '+repr(order)+'\n'+publish[stop:]
publish=publish.replace('第773轮完成（同态BRST配对与来源字典）','第773轮完成（原物质参考与局部BRST）')
publish=publish.replace('## 418.','## 419.').replace('## 323.','## 324.')
publish=publish.replace('同态自由BRST配对与来源字典','原物质参考与局部规范左逆')
publish=publish.replace('旧空间合同保持，同态表示与首阶均值','旧空间合同保持，局部参考与反项类别')
publish=publish.replace('自由物理表示与共同均值','局部规范左逆与形式修复')
write('publish_round773.py',publish)
post=convert((HERE/'postcheck_round772.py').read_text('utf8'))
post=post.replace('第773轮完成（同态BRST配对与来源字典）','第773轮完成（原物质参考与局部BRST）')
write('postcheck_round773.py',post)
names=['research_note_773.md','joint_material_local_brst.py','joint_material_local_brst_results.json','unified_physics_condition_ledger_773.md']
write('round773_drafts/final_review.txt','Primary-agent review; no independent reviewer.\n'
      +'Reviewed local extraction, metric/connection jets, density transpose, antifield splitting, formal-degree recursion and patch/coefficient boundaries.\n'
      +'No original loop anomaly, common QME, strict global UV class, full state or instrument claimed.\n'
      +'\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared 773 science and guarded publication.')
