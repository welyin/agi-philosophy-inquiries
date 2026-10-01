"""Freeze the local area reduction and same-geometry process interface."""
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
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_relational_corner_gluing':'joint_area_reduction_process',
         '647':'648','648':'649','649':'650',
         '3179':'3182','3182':'3185','1253':'1256','1256':'1259',
         '2135':'2144','2144':'2155',
         'joint_interface_entry_audit.md':'physical_inner_product_entry.md',
         'corner_gluing_scope.md':'area_constraint_descent_probe.py'}
verify=replace((HERE/'verify_round648.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original corner source and the specified regular boost gluing obstruction.',
                     'Verify original area quotient, process obstruction and common geometry extension.')
verify=verify.replace("'round649_drafts/area_constraint_descent_probe.py')",
    "'round649_drafts/area_constraint_descent_probe.py','round649_drafts/area_constraint_descent_probe_results.json','round649_drafts/round648_formula14_erratum.md')")
needle="    assert all(c>=32 or c in (10,9) for c in (HERE/'research_note_649.md').read_bytes())"
assert needle in verify
verify=verify.replace(needle,needle+"\n    assert not re.search(r'\\\\[ \\t]+(?:stackrel|frac|sqrt|begin|end|tag|operatorname)\\b',(HERE/'research_note_649.md').read_text('utf8'))")
write('verify_round649.py',verify)
publish=replace((HERE/'publish_round648.py').read_text('utf8'),
                mapping|{'## 294.':'## 295.','## 199.':'## 200.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第649轮完成：** [面积匹配的物理内积、动力学下降与共同量子几何]({p}research_note_649.md)'
         '原面积匹配有局部正物理内积，但复制的原几何动能不保持零范数类；'
         '共享唯一原量子几何可扩充规范切口并保持完整原过程。'
         '三组、十四式通过，最新649／3185，1259份编号科学文件、2155份保护证据。'
         '[核验]({p}research_round_649_checks.json)、[条件账]({p}unified_physics_condition_ledger_649.md)、'
         '[648式14勘误]({p}round649_drafts/round648_formula14_erratum.md)。'
         '限声明的内部几何模型；完整引力约束、连续映射及统一目标仍开放。')
order=('**当前执行顺序（649后，优先于下方历史安排）：** 继续整合共同条件，认知设计后置。'
       '接[650真实引力约束、动态嵌入与区域生成元]({p}round650_drafts/STATUS.md)，'
       '返回原引力作用的实际边界条件，不以控制模型失败代替GR检验，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('物质参考、角点荷与区域匹配','面积约化与同一量子几何的共同过程')
publish=publish.replace('紧群区域拼接不能直接替代引力约化','正物理内积还须容纳原演化')
publish=publish.replace('原物质参考、引力角点与区域量子拼接的共同条件','面积匹配的物理内积、动力学下降与共同量子几何')
write('publish_round649.py',publish)
write('round649_drafts/research_note_649_draft.md',(HERE/'research_note_649.md').read_text('utf8'))
review="""649 primary-agent review; no independent agent review.
648 completed, published and postchecked; previous goal turn is progress.
Navigation confirms648 and2144 evidence files; no Python process found.
Read original590,593,606,617,648 and the saved649 entry/probe. No repeat of
ordinary compact Haar, uncertainty relations or generic encoding as new theorem.
The new Q is expressly area matching on TWO590 geometry-controller replicas.
Their positive Hamiltonian is not ADM and Q was not derived from gravity's
symplectic form. This exact scope is maintained throughout the proof.
The regular patch excludes dQ=0. Compact test support gives a smooth compact
pushforward in Q; its Fourier transform is integrable and coarea yields a
nonzero positive physical norm. No global critical-set uniqueness asserted.
Restriction is unbounded in old L2 norm. Generator descent is invoked only on
common smooth domains, where preserving the physical null class is necessary.
Q^2 chi has zero restriction. All original matter, matrix masses and potentials
retain Q^2, while positive geometry kinetic yields -hbar^2 n^T M n chi.
This is a null-class obstruction, distinct from593's normal-state energy test.
Changing positive surface density cannot cure the same restriction null class.
M positive semidefinite implies n^T M n=0 iff Mn=0; this is necessary only,
not sufficient without first-order drift and domain conditions.
Projected M is a changed theory. Its intrinsic weighted closed form is specified
on a regular bounded smooth subpatch with Dirichlet boundary and the old matter
form domain. No uniqueness, global completion or cognition derivation claimed.
The cross replica block is nonzero, norm1/(2I) at equal configurations.
Numerical coarea is a2D section of the actual original map, not a12D full-state
integral. Analytic coarea carries the full local statement. Gaussian delta is
only an integration regularization. Gradient/Hessian differences are independent.
During primary review a proposed matrix-exponential monotonicity argument was
rejected before freezing. Final proof uses Jacobi minors and a Duhamel bound:
||Q_x||<.24; off-diagonal correction <=.24^2 exp(.48)/6<.016. This avoids the
false claim that matrix exponential is operator monotone. Numerical values unchanged.
The common-geometry positive extension keeps exactly ONE original theta space.
617 J is independent of theta. Geometry gradients commute through J, and
uniform compact-theta bounds close the original full matter and geometry form.
The transported range domain gives exact selfadjoint intertwining; records,
passive references and original form sources follow. No extra controller or
independent gravity tensor factorization is supplied, and it remains non-ADM.
648 formula14 has a backslash-space typo. Frozen original retained; separate
erratum changes no science. New note has a stronger malformed-command check.
3 groups,14 formulas; all ordinary assertions-enabled code checks passed.
No images, no task creation, no goal change. Next returns to actual gravity
constraints and moving boundaries, not further controller optimization.
"""
for name in ('research_note_649.md','joint_area_reduction_process.py',
             'joint_area_reduction_process_results.json','unified_physics_condition_ledger_649.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round649_drafts/final_review.txt',review)
print('649 science frozen; verifier and publisher prepared')
