"""Freeze the full original Gauss heat representation and conditional checks."""
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

mapping={'joint_tree_matter_source':'joint_gauss_fermion_influence',
         '641':'642','642':'643','643':'644',
         '3162':'3165','3165':'3168','1235':'1238','1238':'1241',
         '2089':'2097','2097':'2104'}
verify=replace((HERE/'verify_round642.py').read_text('utf8'),mapping)
verify=verify.replace('Verify full quotient tree transport and same-state source obstruction.',
                      'Verify full Gauss heat influence and original32-CAR conditional weights.')
verify=verify.replace("['display_formulas']==16","['display_formulas']==14")
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_643.md','round643_drafts/research_note_643_draft.md',
           'round643_drafts/final_review.txt','round644_drafts/STATUS.md')
"""+verify[end:]
write('verify_round643.py',verify)
publish=replace((HERE/'publish_round642.py').read_text('utf8'),
                mapping|{'## 288.':'## 289.','## 193.':'## 194.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第643轮完成：** [完整Gauss物质、共同热参考与来源的有序历史表示]({p}research_note_643.md)'
         '原全模型满足矩阵热桥与Gauss闭合的共同可积条件；原32模式验证有序费米权重及来源插入。'
         '三组、十四式通过，最新643／3168，1241份编号科学文件、2104份保护证据。'
         '[核验]({p}research_round_643_checks.json)、[全条件账]({p}unified_physics_condition_ledger_643.md)。'
         '数值为条件因子，未积分全图；连续手征、动态几何及GR仍开放。')
order=('**当前执行顺序（643后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[644共同区域态、物质与几何熵]({p}round644_drafts/STATUS.md)，'
       '核同一热历史的区域拼接、规范边界及来源，不继续条件因子精度扫描；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('完整物质的空间表示与共同来源','完整Gauss物质的共同热历史与来源')
publish=publish.replace('完整Gauss表示不等于物质来源可删除','完整热桥表示不替代连续几何推导')
publish=publish.replace('完整规范物质的共同表示、空间消元与几何来源','完整Gauss物质、共同热参考与来源的有序历史表示')
write('publish_round643.py',publish)
write('round643_drafts/research_note_643_draft.md',(HERE/'research_note_643.md').read_text('utf8'))
review="""643 primary-agent review; no independent agent review.
642 completion and current navigation were checked before643; no Python job was live.
Historical search distinguished529 scalar FK,603 trace comparison,615 auxiliary
Pfaffian and622-625 response from the full original matrix influence now derived.
Guneysu theorem2.11 and2.13 conditions checked: the COMPLETE matrix potential
V I+B has bounded negative part using original W/2 domination; positive part is
smooth locally integrable. B alone is not assumed bounded or Kato-small globally.
Original non-diagonal electric kinetic term defines a complete compact-group
metric and target Laplacians. Heat kernels are relative to original fixed measure,
with the original kinetic normalization. No lost gamma-dependent volume constant.
The bridge runs from q to g^-1 q, with Gamma(g) multiplying the same conditional
Fock propagator. Kernel orientation and late-left ordering follow operator trace.
Scalar kernel Cauchy-Schwarz plus gauge invariance and603 integrability dominate
the absolute Gauss path integral. No pathwise positivity or sign-free Monte Carlo
claim. Full nonlinear fields, gauge links and hopping remain in the analytic formula.
Nambu generator includes antisymmetric Delta and correct conjugate lower block.
Original h is traceless, full R has determinant1; no missing scalar normal-order
factor. The square relation does not select a principal root; original Fock lift
and analytic continuation fix its sign. No continuum chiral Pfaffian imported.
Numerics independently combine24 number-conserving modes with256-dimensional
8-mode paired Fock trace; minors implement the exterior representation. All32
original onsite modes retained. These are CONDITIONAL path factors, not a full
spatial hopping calculation, bosonic bridge integral or physical Gauss average.
Gauge covariance, Z6 representative independence, original Majorana contribution
and noncommuting order are checked. Complex contour factors do not define a
classical real-time Wiener measure or actual record probabilities by themselves.
Source insertion and actual finite difference are only the conditional fermion
lapse derivative. Full geometric derivative retains kinetic and scalar terms by
operator Duhamel. No Brownian measure score assumed across different diffusion
coefficients. Zero weights handled by derivative of the trace, not its logarithm.
Actual bounded record histories use the SAME full H and physical Gibbs state.
Positive-real-time regulators and trace-class continuity establish the real-time
boundary, without claiming numerical full-model record probabilities.
Three groups pass. All27 conditions backfilled; bookkeeping not a separate round.
No new cognition axiom, design task, goal change, continuum, chirality or GR claim.
"""
for name in ('research_note_643.md','joint_gauss_fermion_influence.py',
             'joint_gauss_fermion_influence_results.json','unified_physics_condition_ledger_643.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round643_drafts/final_review.txt',review)
print('643 science frozen; verifier and publisher prepared')
