"""Prepare audited664 and preserve historical navigation without changing goals."""
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


mapping={'joint_spinorial_geometry_sources':'joint_connection_matter_matching',
    '662':'663','663':'664','664':'665','3229':'3231','3231':'3234',
    '1298':'1301','1301':'1304','2289':'2299','2299':'2306'}
verify=replace((HERE/'verify_round663.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original spinorial local lapse, moving geometry and full source dictionary.',
    'Verify original scalar/CAR connection branches, matching and common sources.')
a=verify.index('    import importlib.util');b=verify.index('    result=core.read(model.TARGET)',a)
verify=verify[:a]+"""    import joint_spinorial_geometry_sources as prior_model
    assert prior_model.run()==core.read(prior_model.TARGET)
"""+verify[b:]
verify=verify.replace("==(2,0,0)","==(3,0,0)").replace('fresh_tests=dict(run=2,','fresh_tests=dict(run=3,')
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==16")
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_664.md','round664_drafts/research_note_664_draft.md',
           'round664_drafts/final_review.txt','round665_drafts/STATUS.md')
"""+verify[b:]
write('verify_round664.py',verify)
publish=replace((HERE/'publish_round663.py').read_text('utf8'),mapping|{'## 309.':'## 310.','## 214.':'## 215.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第664轮完成：** [原物质、联络选择与共同几何来源的匹配条件]({p}research_note_664.md)'
         '原Jordan最小独立联络同时改变标量目标和总轴流作用；指定匹配项恢复原经典公共作用。'
         '原CAR核共同接触、排序及几何来源。三组、十六式通过，最新664／3234，1304份编号科学文件、2306份保护证据。'
         '[核验]({p}research_round_664_checks.json)、[条件账]({p}unified_physics_condition_ledger_664.md)。'
         '保留无挠基线；量子测度、连续映射及量子GR仍开放。')
order=('**当前执行顺序（664后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[665同一CAR过程与联络辅助权重]({p}round665_drafts/STATUS.md)，'
       '核接触作用、有序权重、辅助测度及共同来源；目标不改。')
"""+publish[b:]
publish=publish.replace('完整旋量的局部时间组合、空间标架与共同几何来源','原物质、联络选择与共同几何来源的匹配条件')
publish=publish.replace('完整旋量的局部来源与跨背景几何组合','原物质与联络分支的共同匹配')
publish=publish.replace('空间组合须同时输送度规、标架和来源','联络选择同时约束标量几何、物质与来源')
write('publish_round664.py',publish)
write('postcheck_round664.py',replace((HERE/'postcheck_round663.py').read_text('utf8'),mapping))
write('round664_drafts/research_note_664_draft.md',(HERE/'research_note_664.md').read_text('utf8'))
review="""664 primary-agent proof/code audit. No independent agent review.
Previous answer only restated execution order: no progress. This unit executes
the actual original-model connection matching rather than another plan.
Read root/research navigation,663 report/results/postcheck,664 saved STATUS,
358,375,580,598-600,619,643 and actual original F,K,representation,mass,CAR code.
No Python process handle existed at the check. No restart was inferred.
No goal edits, app tasks, automation, image review or agents created.
Source: Karananas et al 2106.13811, equations21,28,34,46-50. Mostly-plus
Lorentz and their axial-current convention stated; global axial sign cancels.
Minimal independent connection in JORDAN action is explicit extra input.
Hermitian spin kinetic, no Holst/nonminimal torsion/gauge-curvature replacement.
Boundary divergence removed only in bulk/compact support or matched boundaries.
Elimination gives +3(dF)^2/4F-3J5^2/16F; conformal factors cancel radial metric
piece and leave -3J5E^2/16. Majorana mass does not source the connection.
Five-field Cartan target scalar curvature independently differentiated from g.
Constant original curvature versus varying Cartan curvature excludes ordinary
scalar-coordinate equivalence with fixed Einstein geometry, not all field changes.
Finite radial boundary distance does not by itself prove non-self-adjointness.
Einstein-frame minimalization is a different declared branch; spin forces neither.
Counterterms cancel in common classical metric/matter variables, not auxiliary
connection observables, quantum determinants or all boundary structures.
All32 physical CAR modes used for currents/internal gauge commutation; no doubling.
Normal-order contraction is 2N. Finite-cell vacuum order is declared, not a
continuum renormalization theorem. Total current must retain cross species.
Original eight-lepton conditional sector retains Dirac and Majorana mass.
Quark vacuum is invariant for this onsite sector. No full Gauss thermal claim.
Charged singlet witness is not mislabelled a globally closed physical state.
Connected neutral diagonal identity rules out a quadratic-CAR replacement.
Hamiltonian sign is minus Lorentz nondifferentiated contact Lagrangian.
w=epsilon^3 exp(6theta); c held fixed; source from same finite thermal trace.
Three groups,16 numbered equations; derivation covers general stated bulk
identities, numerics cover actual project matrices and specified diagnostic values.
Retain original metric branch and all frozen history; no GR-generation claim.
Next665 tests actual auxiliary quantum weight and common sources; no design work.
"""
for name in ('research_note_664.md','joint_connection_matter_matching.py',
             'joint_connection_matter_matching_results.json','unified_physics_condition_ledger_664.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round664_drafts/final_review.txt',review)
print('664 note/review frozen; verification and publication prepared')
