"""Freeze 638: same original record, reference and regional relative entropy."""
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
mapping={'joint_full_thermal_region_entropy':'joint_record_relative_entropy',
         '636':'637','637':'638','638':'639',
         '3149':'3151','3151':'3154','1220':'1223','1223':'1226',
         '2050':'2059','2059':'2067'}
verify=replace((HERE/'verify_round637.py').read_text('utf8'),mapping)
verify=verify.replace('Verify complete-model regional entropy and original quotient heat bounds.',
                      'Verify original record reference support and common relative entropy limits.')
verify=verify.replace("(2,0,0)","(3,0,0)").replace("fresh_tests=dict(run=2,","fresh_tests=dict(run=3,")
verify=verify.replace("['display_formulas']==14","['display_formulas']==16")
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_638.md','round638_drafts/research_note_638_draft.md',
           'round638_drafts/final_review.txt','round639_drafts/STATUS.md',
           'round638_drafts/relative_entropy_entry.md')
"""+verify[end:]
write('verify_round638.py',verify)
publish=replace((HERE/'publish_round637.py').read_text('utf8'),
                mapping|{'## 283.':'## 284.','## 188.':'## 189.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第638轮完成：** [原记录、共同参考与区域相对熵的联合近似]({p}research_note_638.md)'
         '原完整Gibbs与记录共用信息—注能预算；'
         '有限参考配原读取产生无限相对熵，共同保Gauss粗化恢复正确区域极限。'
         '三组、十六式通过，最新638／3154，1226份编号科学文件、2067份保护证据。'
         '[核验]({p}research_round_638_checks.json)、[条件账]({p}unified_physics_condition_ledger_638.md)。'
         '限固定图；统一空间尺度、连续几何熵及GR仍开放。')
order=('**当前执行顺序（638后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[639固定物理支撑、共同几何与信息预算]({p}round639_drafts/STATUS.md)，'
       '核原区域读取能否保尺度统一预算及几何来源；因果实现单列，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原完整能源与区域熵的共同成立','原记录、参考和区域相对熵的共同近似')
publish=publish.replace('固定图有限区域熵不等于连续面积律','谱相对熵极限不替代空间连续匹配')
publish=publish.replace('原完整Gauss热态、区域熵与同尺度近似的共同能源条件','原记录、共同参考与区域相对熵的联合近似')
write('publish_round638.py',publish)
write('round638_drafts/research_note_638_draft.md',(HERE/'research_note_638.md').read_text('utf8'))
review="""638 primary-agent review. No independent agent review.
Previous goal turn made progress: 636 and 637 published and postverified.
README, direction, state, 637 note/results, 638 entry and relevant history read.
Get-Process found no Python worker at entry; no existing process was changed.
305/317 generic entropy identities and 592 finite invariant-space obstruction
are inherited. 578/579 and 613 scale/subtraction boundaries remain in force.
Witten 1803.04993 and Shirokov 2205.10341 Eq1/Theorem2 were checked;
the latter HTML path failed, then its actual PDF text was read without images.
Original full Gibbs is faithful on the complete Gauss Hilbert space.
Original L+/- equals a two-unitary average when outcomes are discarded.
Full-energy control and finite entropy justify entropy concavity and the
same-H Gibbs relative entropy identity; no infinite-minus-infinite step.
Record retention adds the actual Holevo information and cannot be dropped.
Original point-volume budget is not claimed uniform under spatial refinement.
Any nonzero finite reference band plus original uncut L leaks: the proof
uses positivity and the absence of eigenvectors of original L+ multiplication.
This obstruction is compatible with simultaneous ordinary entropy convergence.
The correct comparison applies the SAME CPTP Gauss reset to the actual two
states; it does not re-run the original process on an artificially truncated
reference or claim a free physical reset/cooling operation.
Nested reset composition plus joint lower semicontinuity proves global
relative entropy convergence; fixed regional convergence uses Shirokov Thm2.
No monotonicity is claimed for the regional sequence when maps do not commute.
Numeric matrices are inner-product Gram matrices of original compact Gauss
packets; phase span and common Gram span are NOT called actual H eigenbands.
The explicitly mixed reference zeta is non-Gibbs. Full Gibbs results analytic.
Original five-coordinate inverse metric and random-unitary kernel checked.
Positive physical leakage 8.88323616e-5 yields infinite relative entropy by
support, not a logarithm floor or floating point infinity calculation.
Both 96/160 quadratures agree; tiny positive reference eigenvalue is retained.
Three groups passed first execution, sixteen equations and links checked.
No spatial continuum, boost modular flow, Newton value, GR or autonomous
instrument derived. Next639 tests common physical support and geometry,
with distributed causality explicitly left as an independent requirement.
"""
for name in ('research_note_638.md','joint_record_relative_entropy.py',
             'joint_record_relative_entropy_results.json','unified_physics_condition_ledger_638.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round638_drafts/final_review.txt',review)
print('638 science frozen; verifier and publisher prepared')
