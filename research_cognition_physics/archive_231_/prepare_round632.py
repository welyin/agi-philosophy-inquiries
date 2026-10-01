"""Freeze round 632: finite Weyl contacts and common background matching."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_tensor_stress_spectrum':'joint_background_contact_matching',
         '630':'631','631':'632','632':'633',
         '3130':'3133','3133':'3136','1202':'1205','1205':'1208',
         '2008':'2015','2015':'2022'}
verify=replace((HERE/'verify_round631.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original tensor stress cut and common curvature response.',
                      'Verify shared vacuum potential, local Weyl contacts and background conditions.')
verify=verify.replace("['display_formulas']==14","['display_formulas']==16")
write('verify_round632.py',verify)
publish=replace((HERE/'publish_round631.py').read_text('utf8'),
                mapping|{'## 277.':'## 278.','## 182.':'## 183.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第632轮完成：** [同一真空势、局部接触项与背景方程的共同匹配]({p}research_note_632.md)'
         '原完整质量的有限势给出零导数Weyl项与q—σ接触矩阵；'
         '声明的三个背景匹配条件同时约束物质驻定和几何零动量项。'
         '三组、十六式通过，最新632／3136，1208份编号科学文件、2022份保护证据。'
         '[核验]({p}research_round_632_checks.json)、[条件账]({p}unified_physics_condition_ledger_632.md)。'
         '背景条件为输入；全动量Ward、实际记录后态、图映射及GR仍开放。')
order=('**当前执行顺序（632后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[633连续记录更新、态与共同来源]({p}round633_drafts/STATUS.md)，'
       '复用旧记录结果，核同一连续物质的局域操作和来源；不继续真空参数扫描，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物质张量谱与曲率系数的共同连接','原真空势与几何接触项的共同匹配')
publish=publish.replace('应力自旋2不是引力子或空间生成','共同背景匹配不选择空间维数')
publish=publish.replace('原整代物质的张量应力谱与共同曲率系数','同一真空势、局部接触项与背景方程的共同匹配')
write('publish_round632.py',publish)
write('round632_drafts/research_note_632_draft.md',(HERE/'research_note_632.md').read_text('utf8'))
review="""632 primary-agent review. No independent agent review.
Previous goal turn made verified scientific progress: published 630 and 631.
README, direction, state, 631 note/receipt and frozen 632 entry were reread.
No post-October-1 Python process was observed at the initial process check.
Older Python processes were neither treated as terminated nor restarted.
Historical review found unnumbered round554 vacuum matching work; its generic
shift formula and inability of one constant to remove tadpoles are inherited,
not recast as new theorems. Its source files/results remain unchanged.
New content is original full-mass finite Weyl contact matrix and its relation
to the common background/source conditions following the actual 630/631 cuts.
Markkanen et al 1804.02020 section 3.2 and eq 88 were checked for finite
fermion MS potential and multiplicity. Only the flat constant-mass potential
is used as exact one-fermion-loop background expression; no generic curved
state potential or the full SM contribution is attributed to this paper.
Original all-generation wording is avoided: exactly the full one generation.
Weights (3,3,1,1/2,1/2) equal sixteen Weyl half weights, not spin/Nambu count.
The two neutrino singular-value squares include original Dirac/Majorana mixing.
All masses are independently checked against the original 64-dimensional BdG.
W=F^2 Vf is a rewriting of fixed-Einstein quantum variables, not a second
fixed-Jordan quantization or an assertion of automatic frame equivalence.
Finite Weyl identity has plus C sigma exp(4sigma) lambda^4 on the covariant
density side; V'=4V-C and V''=3V'-4C fix the Hessian differences.
q has lambda'=b and lambda''=1/6. Rank-one nonlocal kernel does not equal
complete zero-frequency Hessian. Local contacts have no new positive cut.
Full anomaly, derivative/curvature contacts, boundaries and other loops
are explicitly outside this zero-derivative calculation.
The .27 diagnostic point is never promoted to a background solution.
Background matching moves the evaluation to the original tree vacuum u.
MS scale mu^2=1 and finite on-shell matching are declared additional inputs.
Delta V_J=-W(u)-gradW(u).(x-u) uses gauge-invariant constant and two quadratic
terms. Uniqueness is only in this three-dimensional counterterm space.
Vacuum value plus two gradients fix all stated zero-derivative metric and
mixed contacts. They do not fix general-momentum Ward identities.
The shear volume Hessian is -2U for h_ij=2 zeta e_ij and tr e^2=1.
The scalar potential Hessian is computed with all original mass derivatives.
Eigenvalues normalized by TREE target metric are diagnostics, not physical
one-loop poles; wavefunction and other-loop corrections are still missing.
At epsilon=1 the inherited linear shift leaves the positive field domain,
so that old formula is not claimed to construct another vacuum.
Finite on-shell matching is a possible conditional branch, not an explanation
of Lambda, vacuum selection, cognitive principles or physical parameter values.
No global/pole/nonperturbative stability follows from local Hessian positivity.
Three test groups passed on first execution. Report has sixteen formulas.
Next 633 integrates actual record-state updates with the same continuous
matter and sources instead of further vacuum or source-precision optimization.
"""
for name in ('research_note_632.md','joint_background_contact_matching.py',
             'joint_background_contact_matching_results.json','unified_physics_condition_ledger_632.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round632_drafts/final_review.txt',review)
print('632 verification and publication prepared; science frozen')

