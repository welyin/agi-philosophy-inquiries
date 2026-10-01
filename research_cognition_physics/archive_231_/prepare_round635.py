"""Freeze 635: original mass and target-gradient geometric entropy."""
import hashlib,re
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
mapping={'joint_record_ward_boundary':'joint_entropy_geometric_matching',
         '633':'634','634':'635','635':'636',
         '3140':'3143','3143':'3146','1211':'1214','1214':'1217',
         '2029':'2036','2036':'2043'}
verify=replace((HERE/'verify_round634.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original record charge noise and distributional source boundary.',
                      'Verify original mass spectra, replica UV entropy and target-gradient term.')
verify=verify.replace("['display_formulas']==14","['display_formulas']==12")
write('verify_round635.py',verify)
publish=replace((HERE/'publish_round634.py').read_text('utf8'),
                mapping|{'## 280.':'## 281.','## 185.':'## 186.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第635轮完成：** [原物种、曲目标梯度与区域熵的共同几何系数]({p}research_note_635.md)'
         '原五标量和整代费米质量共同匹配面积与R项；'
         '原RS项强制给梯度几何熵，同场值同面积不再由一个常数系数覆盖。'
         '三组、十二式通过，最新635／3146，1217份编号科学文件、2043份保护证据。'
         '[核验]({p}research_round_635_checks.json)、[条件账]({p}unified_physics_condition_ledger_635.md)。'
         '限部分一圈及局部UV；规范边界、有限熵、Newton值和GR仍开放。')
order=('**当前执行顺序（635后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[636原Gauss区域、边界通量与规范熵]({p}round636_drafts/STATUS.md)，'
       '回查617原商群切分及旧熵结果，对接共同态与边界模式；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原记录四动量与几何守恒的共同边界','原物种和曲目标的共同几何熵')
publish=publish.replace('记录来源匹配不完成空间或GR推导','几何熵匹配不选择维数或Newton常数')
publish=publish.replace('记录条件化、守恒来源与几何约束的共同边界','原物种、曲目标梯度与区域熵的共同几何系数')
write('publish_round635.py',publish)
write('round635_drafts/research_note_635_draft.md',(HERE/'research_note_635.md').read_text('utf8'))
review="""635 primary-agent review. No independent agent review.
Previous goal turn made scientific progress: 633-634 published and postverified.
Latest README, direction, state, 634 note/receipt and 635 entry were reread.
No October-1 Python worker was observed; older processes were left untouched.
An initial guessed heat-kernel filename was absent; the 599 note supplied the
actual joint_fermion_scalar_loop_matching.py path, which was then read.
313-315 already prove generic scalar entropy/G matching and splitting.
Those theorems are inherited, not counted anew. 580/599 supply original
scalar Jacobi and fermion curvature coefficients; 631/632 scope retained.
Solodukhin 1104.3712 sections 3.12, 3.14-3.17 were read for first replica
variation and spin weights. The historical review's old gauge-puzzle wording
is not treated as current: Donnelly-Wall 1506.05792 is explicitly acknowledged.
The new mapping uses original five scalar Hessian and complete one-generation
masses, with seven Dirac plus two Majorana species and correct half weights.
Tree Hessian, not the corrected one-loop Hessian, enters the one-loop determinant.
Gauge fields and ghosts are NOT integrated; 21 is a partial determinant weight,
not total physical SM degrees of freedom or a prediction about the universe.
Independent S2 Dirac heat spectrum checks -R/3 and original R*m^2 coefficient.
Proper-time integrals are checked against E1 series; large-s integration is
only used on the constant free planar branch, not arbitrary curved backgrounds.
The mixed gradient entropy uses the original H5 Jacobi trace trA+2S/3,
which matches inherited RS through the first conical variation.
Tangential background gradients only, normal derivatives zero, planar
surface, zero extrinsic curvature, no outer boundary are explicit conditions.
Nonconstant backgrounds are only used for local UV coefficients; no complete
off-shell replica state, full finite entropy or global stability claimed.
Pure-curvature terms are not discarded: their linear response vanishes for
this flat planar sample, and must be restored for general surfaces.
The two backgrounds have the SAME local field and masses at x=0, so one
field-only or constant area coefficient cannot match the gradient correction.
This comparison is of local jets, not assumed global background solutions.
Derivative field redefinitions require changing area, sources and state;
the original surface is not reused after erasing RS in the action.
Jordan/Einstien equality is a rewrite of SAME quantum variables and local term,
not independent quantizations or a complete frame-anomaly equivalence claim.
Finite G remains input; no universal physical cutoff or entropy balance derived.
Hadamard records inherit local counterterms, but finite entropy differences
are not claimed to follow solely from the Hadamard two-point property.
Three groups passed first execution; 12 formulas checked without image rendering.
Next 636 integrates the actual quotient gauge region and edge entropy after
history audit; it does not open cognition hardware design. Goal unchanged.
"""
for name in ('research_note_635.md','joint_entropy_geometric_matching.py',
             'joint_entropy_geometric_matching_results.json','unified_physics_condition_ledger_635.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round635_drafts/final_review.txt',review)
print('635 science frozen; verifier and publisher prepared')

