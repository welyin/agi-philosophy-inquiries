"""Freeze the original corner/source and specified boost-gluing interface."""
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


mapping={'joint_relational_readout_source':'joint_relational_corner_gluing',
         '646':'647','647':'648','648':'649',
         '3177':'3179','3179':'3182','1250':'1253','1253':'1256',
         '2127':'2135','2135':'2144',
         'relational_quantum_entry.md':'joint_interface_entry_audit.md',
         '2026-10-01':'2026-10-02','20261001':'20261002'}
verify=replace((HERE/'verify_round647.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original material reference, readout and joint location sources.',
                     'Verify original corner source and the specified regular boost gluing obstruction.')
verify=verify.replace('==(2,0,0)','==(3,0,0)').replace('fresh_tests=dict(run=2,','fresh_tests=dict(run=3,')
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==14")
verify=verify.replace("'round648_drafts/joint_interface_entry_audit.md')",
                      "'round648_drafts/joint_interface_entry_audit.md','round648_drafts/corner_gluing_scope.md')")
write('verify_round648.py',verify)
publish=replace((HERE/'publish_round647.py').read_text('utf8'),
                mapping|{'## 293.':'## 294.','## 198.':'## 199.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第648轮完成：** [原物质参考、引力角点与区域量子拼接的共同条件]({p}research_note_648.md)'
         '原F同时固定角点面积荷与标量来源，原物质参考给同源局部角点；'
         '指定正则boost拼接无非零正常不变态，有限规范区间近似没有强极限。'
         '三组、十四式通过，最新648／3182，1256份编号科学文件、2144份保护证据。'
         '[核验]({p}research_round_648_checks.json)、[条件账]({p}unified_physics_condition_ledger_648.md)。'
         '结论限经典角点及声明的量子接法，未推出完整引力量子态或统计面积律。')
order=('**当前执行顺序（648后，优先于下方历史安排）：** 先整合共同条件，认知系统设计后置。'
       '接[649物理区域内积与同一过程的下降]({p}round649_drafts/STATUS.md)，'
       '核非紧约化的替代方案与原过程是否相容；不扫描规范区间精度，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物质参考与关系定位的共同来源','物质参考、角点荷与区域匹配')
publish=publish.replace('读取、体积与定位必须共同输送','紧群区域拼接不能直接替代引力约化')
publish=publish.replace('原物质参考、关系读取与几何定位来源','原物质参考、引力角点与区域量子拼接的共同条件')
write('publish_round648.py',publish)
write('round648_drafts/research_note_648_draft.md',(HERE/'research_note_648.md').read_text('utf8'))
review="""648 primary-agent review; no independent agent review.
647 was published and postchecked; previous goal turn is verified progress.
Current navigation and evidence confirm647, with no Python process observed.
Historical search covered references, field-dependent regions and corner terms.
554 joint-reference uncertainty and606 general encoding derivatives were not
repeated.617 uses compact internal Haar,618 explicitly excluded corners.
Original four-dimensional two-derivative FR/2, F>0 and five scalars are inputs.
Jubb et al corner normalization is R/2; orientation and angle origin are fixed.
The Einstein frame transforms a2D area by F and leaves the normal angle fixed.
Curvature contraction independently gives -F with epsilon_ab epsilon^ab=-2.
Noether unit-boost density is F sqrt(q_J). Terms xi grad F vanish at the corner.
Original gauge background remains; pure diffeo YM corner term has xi=0. Internal
Gauss charge still matches separately. Fermion/torsion corners are not proved.
Corner variations retain angle, induced metric and all scalar derivatives.
Fixed-embedding formulas do not replace647 moving-location derivatives.
No independent corner counterterm or new cognition axiom has been introduced.
573 actual h,s intersection has tangent x,z at the stated original point. The
two time gradients remain same oriented and timelike at lambda=.005. Gram
and orthonormal-frame rapidities agree. Lambda selects boundaries, not dynamics.
Classical source is on shell; diagnostic variations need not solve constraints.
This ordinary local corner is not declared a horizon. Wald charge normalization
does not derive statistical entropy or the637-644 regional entropy equality.
Donnelly-Freidel is used for the interface, not claimed to supply gravity Hilbert
space. Only a fixed-normal R boost subgroup is analyzed, not full SL(2,R).
The L2(R) regular frame representation is an explicit candidate, not a deduction
from original quantum gravity. No normalized Haar exists for this noncompact R.
Invariant L2(R2) functions are zero by u=a+b and infinite a volume. This rejects
the specific direct617 extension, not all noncompact constraint quantization.
J_L normalized windows have exact fixed-t defect and J_L/J_4L distance1, proving
absence of a strong limit despite pointwise-in-t invariance. Cutoffs and Gaussian
test functions are diagnostics, not finite-energy physical gravity states.
Positive area restrictions within that kinematical space cannot create normal
invariants; alternative representations and rigged inner products remain open.
Unnormalized group averaging is only the next candidate, not a proven original
model solution. Mere generic Fourier reduction must not be called another round.
All3 numerical groups passed on first ordinary assertions-enabled execution.
14 equations,3 groups, no images, unchanged goal; cognition design stays later.
"""
for name in ('research_note_648.md','joint_relational_corner_gluing.py',
             'joint_relational_corner_gluing_results.json','unified_physics_condition_ledger_648.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round648_drafts/final_review.txt',review)
print('648 science frozen; verifier and publisher prepared')
