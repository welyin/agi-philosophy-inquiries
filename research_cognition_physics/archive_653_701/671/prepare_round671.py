"""Prepare671 publication with inherited evidence and fresh scope audit."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_670.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：670非平坦规范、共同观测与无质量零模',
    '# 联合条件总账：671静态规范曲率与原联合反射泛函')
ledger=ledger.replace('接[669全账](unified_physics_condition_ledger_669.md)，回填[670报告](research_note_670.md)。[结果](joint_nonflat_mass_measure_results.json)、[核验](research_round_670_checks.json)。',
    '接[670全账](unified_physics_condition_ledger_670.md)，回填[671报告](research_note_671.md)。[结果](joint_static_gauge_reflection_results.json)、[核验](research_round_671_checks.json)。')
ledger=ledger.replace('670非平坦共同字典跨无质量零模|',
    '670正则字典；671静态空间曲率下的原联合未归一反射泛函|')
ledger=ledger.replace('670原观测／质量有无逆K的正则字典；原物理态与辅助泛函的同一完整表示仍缺',
    '670正则观测；671静态背景物理／辅助／质量共用反射泛函；原物理态与辅助过程同一性仍缺')
ledger=ledger.replace('670非交换非平坦背景共同协变已核；群选择及真实规范动力学仍缺',
    '670非平坦协变、671静态曲率正性已核；群选择与动态规范过程仍缺')
ledger=ledger.replace('动态规范正性、指定归一与全局测度仍缺|',
    '671静态背景全部多项式未归一反射正性已证；动态规范正性、指定归一与全局测度仍缺|')
ledger=ledger.replace('动态过程、机制选择及现实参数|',
    '671静态背景原质量反射延拓；动态过程、机制选择及现实参数|')
ledger=ledger.replace('原辅助泛函与物理态识别；669已接原H_b配置热矩，',
    '671静态辅助／物理共同正泛函不等于原Gibbs；669已接原H_b配置热矩，')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 671将自由作用核的反射输入扩展至原静态非交换空间规范背景，不要求空间Wilson项与差分对易。归一物理基泛函仍需D／Kₗ可逆；这与670未归一观测字典不需逆K的结论不同，不混写。

## 本轮合并与下一项

671合并C01、C14、C16、C17、C19：矩阵谱割线给原静态规范背景的正作用核，继承反射锥及原局部辅助配对，接到物理Weyl／E联合代数和全部原质量的未归一反射泛函。严格归一仍需对应整体权重非零，不能将668自由全λ结果扩大到任意曲率背景。

2026原作者报告仍将相应结论限定于弱耦合；605已核Creutz硬允许域与转移正性的准确范围，不重算为新研究。新静态证明保留原全部内部矩阵及真实空间曲率；有限Gram和球面配置核仅是实现检查，全多项式论证在671正文。

动态链路在X†X中增加时间交换项，静态时间谱推导不能直接使用。该断点不是动态理论反证。原643 Gauss热历史、624仪器、606—610投影和电动能身份保持；静态背景混合分布不能替代这些原动态对象。

接[672](round672_drafts/STATUS.md)：核反射兼容的时变链路、原共同来源与完整正时间拼接。继续寻找动态规范与手征测度的实际映射，不继续自由参数扫描或认知装置设计。旧空间路线与已消去条件完整继承，四分支未统一，目标不改。
'''
write('unified_physics_condition_ledger_671.md',ledger)
mapping={'joint_nonflat_mass_measure':'joint_static_gauge_reflection',
         '669':'670','670':'671','671':'672','3245':'3247','3247':'3249',
         '1319':'1322','1322':'1325','2364':'2376','2376':'2389'}
verify=replace((HERE/'verify_round670.py').read_text('utf8'),mapping)
verify=verify.replace('Verify nonflat shared dictionary and regular physical observations at massless zero.',
    'Verify static noncommuting spectral density and original joint reflection.')
verify=verify.replace('import joint_dynamic_scalar_gauge_interface as prior_model','import joint_nonflat_mass_measure as prior_model')
verify=verify.replace('nonflat_gauge_probe','magnetic_reflection_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_671.md','round671_drafts/research_note_671_draft.md',
           'round671_drafts/final_review.txt','round672_drafts/STATUS.md',
           'round671_drafts/magnetic_reflection_entry.md','round671_drafts/magnetic_reflection_probe.py',
           'round671_drafts/magnetic_reflection_probe_results.json','round671_drafts/mature_scope_audit.json',
           'round671_drafts/first_static_attempt.py','round671_drafts/first_static_diagnostic.json')
"""+verify[b:]
write('verify_round671.py',verify)
publish=replace((HERE/'publish_round670.py').read_text('utf8'),mapping|{'## 316.':'## 317.','## 221.':'## 222.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第671轮完成：** [静态规范曲率、原手征物质与共同反射泛函]({p}research_note_671.md)'
         '非交换矩阵谱将自由反射输入接到原静态空间规范曲率；物理Weyl、辅助及质量共用未归一正泛函。'
         '严格归一及动态规范过程仍单独验收；时间变化新增交换项已定位。'
         '两组、十八式通过，最新671／3249，1325份编号科学文件、2389份保护证据。'
         '[核验]({p}research_round_671_checks.json)、[全条件账]({p}unified_physics_condition_ledger_671.md)。'
         '原完整过程、动态规范正性、连续及量子GR仍开放。')
order=('**当前执行顺序（671后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[672动态规范链路与原共同反射过程]({p}round672_drafts/STATUS.md)，'
       '核时变链路的完整反射拼接及原过程映射；目标不改。')
'''+publish[b:]
publish=publish.replace('非平坦规范场中的共同观测，以及无质量零模处的正则延伸','静态规范曲率、原手征物质与共同反射泛函')
publish=publish.replace('共同规范字典消去无质量零模限制','静态曲率与原物理辅助质量的共同反射条件')
publish=publish.replace('正则观测字典与动态规范正性仍须区分','静态反射泛函与真实动态规范过程的边界')
write('publish_round671.py',publish)
write('postcheck_round671.py',replace((HERE/'postcheck_round670.py').read_text('utf8'),mapping))
write('round671_drafts/research_note_671_draft.md',(HERE/'research_note_671.md').read_text('utf8'))
review='''671 primary-agent proof/code audit; no independent agent review.
Previous goal turn completed670 and actually executed671 entry: progress.
Current state read; no Python job found. Goal intact; no new tasks or automations,
no image verification; existing Python/NumPy only. Historical evidence preserved.
2010 primary RP cone method and 2026 author slides scope rechecked.605 already
contains the Creutz admissibility restriction and is reused, not a new theorem.
New bridge: static spatial gauge links, identity temporal transport,0<m0<=1.
B>0 for m0<1, C anti-Hermitian, gamma0 commutes B and anticommutes C.
Hs Hermitian; X^dag X=1+Hs^2-2Bcos p including nonzero [B,C].
Eta proof regularization avoids real-axis branch singularities. Positive/negative
imaginary part of Q in the two strips locates cuts on imaginary axis; at pi+iE
Q positive. On the cut Q=LE^2+eta^2-sinh(E)^2 so density is a function of LE.
Both signed densities PSD on |LE|<sqrt(sinh^2-eta^2). Large arcs vanish for
integer separation>=1. AP image coefficients are positive in both directions.
Finite AP continuity first eta->0 then m0->1; requires final finite Wilson gap.
No all-real-temporal-momenta gap or commuting spatial matrices newly assumed.
Directed action cross sign checked against657; initial wrong orientation retained.
Near-degenerate residue diagnostics replaced by whole-cut matrix function rather
than relaxed tolerances. Zero density support handled explicitly; cut finite
difference near threshold refined once with unchanged error gate and convergence.
Cone argument gives full Dirac RP; original local pair and compact S9 retain it.
Original auxiliary marginal includes complete detKl/detS phase. Physical Weyl
normalized subalgebra uses D/Kl invertibility, not a restriction of670's regular
unnormalized map. Tensor/Schur argument and physical local even mass give joint
unnormalized RP. Strict normalization all backgrounds/couplings not asserted.
Actual noncommuting original16 representation, four time slices, mass and signed
auxiliary kernels checked; finite samples not full sphere integration or proof.
Dynamic X^dag X has explicit temporal commutators; static proof does not decide
their full RP.605-610,643,624, old space and alternative bridges remain inherited.
No equality with original H_F, instrument, continuum or quantum GR claimed.
Two groups,18 equations; next672 restores actual time dependence and joint process.
'''
for name in ('research_note_671.md','joint_static_gauge_reflection.py','joint_static_gauge_reflection_results.json','unified_physics_condition_ledger_671.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round671_drafts/final_review.txt',review)
print('671 report/review and publication prepared')
