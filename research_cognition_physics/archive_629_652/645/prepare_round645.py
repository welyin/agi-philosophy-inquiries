"""Freeze the shared thermal-reference, mass-force and background result."""
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


mapping={'joint_region_replica_source':'joint_thermal_background_source',
         '643':'644','644':'645','645':'646',
         '3168':'3171','3171':'3174','1241':'1244','1244':'1247',
         '2104':'2111','2111':'2119'}
verify=replace((HERE/'verify_round644.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original full Gauss replica identities and64-CAR diagnostics.',
                      'Verify original continuum thermal mass and geometric source matching.')
verify=verify.replace("'round646_drafts/STATUS.md')\n    preserved=",
                      "'round646_drafts/STATUS.md','round645_drafts/reference_scale_entry.md')\n    preserved=")
verify=verify.replace("text['display_formulas']==14", "text['display_formulas']==16")
write('verify_round645.py',verify)
publish=replace((HERE/'publish_round644.py').read_text('utf8'),
                mapping|{'## 290.':'## 291.','## 195.':'## 196.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第645轮完成：** [共同热参考、原质量与几何背景的联合条件]({p}research_note_645.md)'
         '原全代质量的同一热函数共同固定标量力、应力和几何接触项；'
         '固定β的lapse必须输送局部温度，仅调真空势不能保持平直常场热背景。'
         '三组、十六式通过，最新645／3174，1247份编号科学文件、2119份保护证据。'
         '[核验]({p}research_round_645_checks.json)、[条件账]({p}unified_physics_condition_ledger_645.md)。'
         '限指定连续自由费米分支；图态映射、完整熵与GR仍未完成。')
order=('**当前执行顺序（645后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[646跨分支共同区域、尺度与几何资料]({p}round646_drafts/STATUS.md)，'
       '核实际共同映射；不继续热势或装置优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('完整Gauss区域态的二阶熵与共同来源','共同热参考的质量与几何背景条件')
publish=publish.replace('区域复制一致不等于面积熵或GR成立','有限热参考仍须匹配真实几何来源')
publish=publish.replace('同一Gauss区域态的费米复制、二阶熵与几何来源','共同热参考、原质量与几何背景的联合条件')
write('publish_round645.py',publish)
write('round645_drafts/research_note_645_draft.md',(HERE/'research_note_645.md').read_text('utf8'))
review="""645 primary-agent review; no independent agent review.
644 is complete, published and postverified. Entry audit is separately frozen.
Scope is the same original full-generation continuum mass branch as632, not
the full graph Gauss thermal integral. Its continuum 3+1 and state are inputs.
The positive masses at ORIGINAL tree background u were verified against all32
positive original64 BdG eigenvalues. g=4d totals32, no second Nambu multiplier.
Neutral eigenmasses retain Majorana/Dirac mixing and complex Yukawa moduli.
Thermal integrals are standard (Quiros201-203), not a newly derived theory.
Hydrostatic zero-derivative functional and redshift are imported with their
static, slowly varying equilibrium scope (Banerjee1.3-1.5 and2.5). No thermalization.
Thermal-minus-vacuum f is UV finite and independent of mu, but the632 vacuum
Weyl-log remains. No claim about finiteness of full regional entanglement.
Integrating pressure by parts, theta=rho-3p and mass derivative agree. Theta
is MINUS the Lorentz stress trace under (-+++), explicitly stated in note.
Original mass homogeneity gives phi.grad m=(M/F)m; this exact identity retains
the neutral mixing. It is the new original-model source connection, not the
old554/632 generic minimum shift formula. No new thermal minimum is claimed.
Fixed coordinate beta makes Tloc=T0/N. Lapse derivative is rho, spatial volume
derivative is -3p, uniform Weyl derivative is theta, NOT4f. Their difference
is rho+p. At constant scalar/geometry this is exact thermal density physics;
curved-space use is explicitly only zeroth derivative, not full determinant.
Finite thermal homogeneity yields full q-sigma Hessian including F1 contacts.
F2 alone is not a complete source Hessian nor a positive-noise claim.
Scalar gradient includes both original neutral eigenvalue and F derivatives.
q uses the original target-metric normalization at u, not630 diagnostic q=.27.
No-gauge/scalar-loop approximation is explicit; no high-temperature expansion.
Flat static compensation obstruction assumes constant scalar, specified Einstein
leading equation and ONLY vacuum potential additions. Such additions have
rho+p=0 and cannot cancel strictly positive thermal rho+p. Curved/dynamic
geometry, extra stresses and other routes remain open. No GR no-go asserted.
Infinite-domain quadrature192/256, massless count, thermodynamic/geometry/scalar
finite differences and one hydrostatic Ward direction are three check groups.
Numerical errors are diagnostic, identities are proved analytically.
No image checks, design optimization, goal change or stage-completion claim.
Next step returns to the actual cross-branch map, not more thermal scans.
"""
for name in ('research_note_645.md','joint_thermal_background_source.py',
             'joint_thermal_background_source_results.json','unified_physics_condition_ledger_645.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round645_drafts/final_review.txt',review)
print('645 science frozen; verifier and publisher prepared')
