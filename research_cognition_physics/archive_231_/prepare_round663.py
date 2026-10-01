"""Freeze663 actual spin/frame/source matching and prepare publication."""
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


mapping={'joint_local_lapse_fermion_current':'joint_spinorial_geometry_sources',
    '661':'662','662':'663','663':'664','3227':'3229','3229':'3231',
    '1295':'1298','1298':'1301','2279':'2289','2289':'2299'}
verify=replace((HERE/'verify_round662.py').read_text('utf8'),mapping)
verify=verify.replace('Verify full original local-lapse CAR current, mixed electric term and record interface.',
    'Verify original spinorial local lapse, moving geometry and full source dictionary.')
verify=verify.replace('local_fermion_lapse_probe','spinorial_local_source_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_663.md','round663_drafts/research_note_663_draft.md',
           'round663_drafts/final_review.txt','round664_drafts/STATUS.md',
           'round663_drafts/local_spin_geometry_entry.md','round663_drafts/spinorial_local_source_probe.py',
           'round663_drafts/spinorial_local_source_probe_results.json')
"""+verify[b:]
write('verify_round663.py',verify)
publish=replace((HERE/'publish_round662.py').read_text('utf8'),mapping|{'## 308.':'## 309.','## 213.':'## 214.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第663轮完成：** [完整旋量的局部时间组合、空间标架与共同几何来源]({p}research_note_663.md)'
         '原连续完整物质的局部时间组合接到spin切向输送；同步几何字典恢复有限组合及来源导数。'
         '两组、十四式通过，最新663／3231，1301份编号科学文件、2299份保护证据。'
         '[核验]({p}research_round_663_checks.json)、[条件账]({p}unified_physics_condition_ledger_663.md)。'
         '给定几何与spin结构；真实图态映射及量子GR仍开放。')
order=('**当前执行顺序（663后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[664原物质的联络来源与几何变分]({p}round664_drafts/STATUS.md)，'
       '核无挠／独立联络分支、完整来源和有效相互作用；目标不改。')
"""+publish[b:]
publish=publish.replace('完整费米物质、局部时间来源与原记录的共同接口','完整旋量的局部时间组合、空间标架与共同几何来源')
publish=publish.replace('完整局部时间来源与原记录的共同域','完整旋量的局部来源与跨背景几何组合')
publish=publish.replace('完整费米来源不能由纯玻色流代替','空间组合须同时输送度规、标架和来源')
write('publish_round663.py',publish)
write('postcheck_round663.py',replace((HERE/'postcheck_round662.py').read_text('utf8'),mapping))
write('round663_drafts/research_note_663_draft.md',(HERE/'research_note_663.md').read_text('utf8'))
review="""663 primary-agent proof/code audit, no independent agent review.
Previous goal turn completed662 and saved actual663 probe: progress.
Read root/research navigation,662 report/results/postcheck,663 saved entry,
375,600/601,604,611/612/613,619 and actual original model definitions.
Ordinary process query found no Python; no inferred live handle or restart.
No goals, app tasks, automation, visual inspection or subagents created.
Full32CAR/64Nambu Gamma and mass inherited604, not another Clifford discovery.
Naive lattice doubling and exact-time Gibbs mismatch remain excluded contracts.
Positive Hermitian spatial Gamma convention fixed. Kosmann spatial Clifford
translation c=iGamma; independent rigid rotation fixes sign and lift direction.
Local commutator /i yields i K_v, not silently the opposite Hermitian D_v.
Compatible spin/gauge connection and given positive metric are explicit inputs.
Mass cancels even for coordinate-dependent B with the same anticommutation.
Core differential identity is not regulated QFT commutator or anomaly proof.
Half-density volume factor is distinct from four-dimensional Weyl weight.
Affine vector bracket is (BA-AB)x. Fixed-g defect -S([sym A,sym B]) verified.
Counterexample restricted to fixed-background endomorphism composition.
Actual repair uses g'=F^T g F and R=sqrt(g) F inv(sqrt(g')), not fitted spin.
Second composition uses updated metric. Local continuous Spin branch only.
At flat base, derivative of finite pullback recovers original K(Ax), proving
repair concerns the same tangent dictionary. General frame choices not unique.
Full original mass, nonorthogonal lapse symbols and bracket transformed together.
Metric-dependent U requires derivative commutator; finite difference independently
checks its nonzero contribution. Bare-source comparisons alone are not physical.
Formal quantum states, actual original scalar instruments and graph-continuum
maps remain open. No complete Einstein/ADM dynamics follows from a cocycle.
2 groups,14 equations; existing entry probe reproduced without extra count.
Next664 audits original spin source and metric/independent connection branches;
does not force Einstein-Cartan from spin or design cognitive devices.
"""
for name in ('research_note_663.md','joint_spinorial_geometry_sources.py',
             'joint_spinorial_geometry_sources_results.json','unified_physics_condition_ledger_663.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round663_drafts/final_review.txt',review)
print('663 note/review frozen; checks/publication prepared')
