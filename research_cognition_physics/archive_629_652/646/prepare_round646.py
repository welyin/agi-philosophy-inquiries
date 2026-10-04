"""Freeze same-spectrum continuum windows and the coincident source boundary."""
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


mapping={'joint_thermal_background_source':'joint_overlap_continuum_sources',
         '644':'645','645':'646','646':'647',
         '3171':'3174','3174':'3177','1244':'1247','1247':'1250',
         '2111':'2119','2119':'2127',
         'reference_scale_entry.md':'cross_branch_entry_audit.md'}
verify=replace((HERE/'verify_round645.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original continuum thermal mass and geometric source matching.',
                      'Verify original overlap continuum window and coincident energy source boundary.')
write('verify_round646.py',verify)
publish=replace((HERE/'publish_round645.py').read_text('utf8'),
                mapping|{'## 291.':'## 292.','## 196.':'## 197.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第646轮完成：** [原传播核的共同连续窗口与重合来源边界]({p}research_note_646.md)'
         '原自由overlap谱在非零物理时间窗口共同匹配关联及两阶能源插入；'
         '同一谱的重合二阶能源矩有发散下界，分离时间收敛不能代替接触匹配。'
         '三组、十六式通过，最新646／3177，1250份编号科学文件、2127份保护证据。'
         '[核验]({p}research_round_646_checks.json)、[条件账]({p}unified_physics_condition_ledger_646.md)。'
         '限自由核；完整Gauss区域态、手征连续和动态引力仍开放。')
order=('**当前执行顺序（646后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[647固定图过程与动态几何的共同条件]({p}round647_drafts/STATUS.md)，'
       '回填共同对象和来源，先查旧约束；不继续自由核精度优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('共同热参考的质量与几何背景条件','同一原谱的有效连续窗口与来源边界')
publish=publish.replace('有限热参考仍须匹配真实几何来源','分离时间连续不替代重合来源匹配')
publish=publish.replace('共同热参考、原质量与几何背景的联合条件','原传播核的共同连续窗口与重合来源边界')
write('publish_round646.py',publish)
write('round646_drafts/research_note_646_draft.md',(HERE/'research_note_646.md').read_text('utf8'))
review="""646 primary-agent review; no independent agent review.
645 was complete and published, followed by a2127-predecessor check of2119
protected files and7 navigation snapshots. Prior turn is progress, not a wait.
Entry audit records the613 precedent for645's general vacuum-pressure warning.
646 keeps the ORIGINAL612 free massless m0=1 inverse-overlap kernel. It does not
substitute the613 spatial Hamiltonian or an interacting Gauss state.
The torus has fixed physical side2pi and integer p, a=2pi/N. Tested time points
have integer t/a at every N, so Fourier comparisons are actual lattice data.
Off-grid time and0+ are explicitly the612 positive cyclic spectral extension,
NOT the raw lattice equal-time contact value. No equal-time definition changed.
Rational-y change of variable is derived from cosh^2(E)-B^2 and its Jacobian.
Integral gives delta=b/(2B), w+delta=1/2 exactly, not a numerical inference.
For aP<=1/2, b<=a^2P^2/2<=1/8 and B>=1. Threshold >=log(1/b).
Pole-energy error follows two independent sin/asinh remainder inequalities.
Uniform cut energy moments use sup x^j exp(-tau x) on x>=L, with tau L>=j.
Pole derivative bound is on0<=epsilon,r<=P and tau<=t<=T. Combined bound covers
j0,1,2 together; P,tau are matching windows, not new hard physical thresholds.
Constant lapse rescales the SAME cyclic Hamiltonian with the vector fixed;
its derivatives are energy insertions, not a full stress tensor or GR source.
At0+ positive threshold bound gives divergent second energy moment at fixed p.
The first moment tends0 for the cut by uniformly integrable rational-y density
and E<=log(1/b)+log(2+y^2). Zero moment equals1/2 exactly.
The noncommuting second-derivative limit is for this specified bare source,
not an assertion that renormalized physical stress or quantum gravity diverges.
Finite-a continuous spectrum and non-trace-class612 cyclic heat remain.
Window-weight convergence cannot delete Hilbert states or prove603 Gibbs/KMS.
Bandlimited observables are not compact spatial regions or617/644 entropy maps.
Quadrature smooths the threshold; E1+90 is computational only and a positive
analytic tail bound is retained. Tail rho<=sqrt(b)exp(-E/2) follows from E>=log8
and s^2<=1/4, with ample constant margin. j<=2 incomplete-gamma tail is exact.
Two quadrature orders and an independent Fourier integral validate the code;
numerical samples do not prove sup bounds or divergence, the note proves them.
Known reflection positivity is cited with free/non-gauge scope and not claimed new.
3 groups and16 formulas; no images or design optimization. Goal remains active.
Next returns to the full joint model and geometry, not higher-order source scans.
"""
for name in ('research_note_646.md','joint_overlap_continuum_sources.py',
             'joint_overlap_continuum_sources_results.json','unified_physics_condition_ledger_646.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round646_drafts/final_review.txt',review)
print('646 science frozen; verifier and publisher prepared')
