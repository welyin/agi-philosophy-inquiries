"""Prepare662 verification and append-only publication; no frozen edits."""
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


mapping={'joint_physical_auxiliary_state':'joint_local_lapse_fermion_current',
    '660':'661','661':'662','662':'663','3225':'3227','3227':'3229',
    '1292':'1295','1295':'1298','2269':'2279','2279':'2289'}
verify=replace((HERE/'verify_round661.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original physical/auxiliary reflection state and common source.',
    'Verify full original local-lapse CAR current, mixed electric term and record interface.')
verify=verify.replace('observable_interface_probe','local_fermion_lapse_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_662.md','round662_drafts/research_note_662_draft.md',
           'round662_drafts/final_review.txt','round663_drafts/STATUS.md',
           'round662_drafts/local_source_entry.md','round662_drafts/local_fermion_lapse_probe.py',
           'round662_drafts/local_fermion_lapse_probe_results.json')
"""+verify[b:]
write('verify_round662.py',verify)
publish=replace((HERE/'publish_round661.py').read_text('utf8'),mapping|{'## 307.':'## 308.','## 212.':'## 213.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第662轮完成：** [完整费米物质、局部时间来源与原记录的共同接口]({p}research_note_662.md)'
         '原完整固定图的正局部lapse族共用算符域、Gauss与原记录；完整能源流保留电微分及CAR矩阵补项。'
         '两组、十四式通过，最新662／3229，1298份编号科学文件、2289份保护证据。'
         '[核验]({p}research_round_662_checks.json)、[条件账]({p}unified_physics_condition_ledger_662.md)。'
         '局部分配、连续手征映射及量子GR仍开放。')
order=('**当前执行顺序（662后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[663局部来源与连续手征、几何的共同字典]({p}round663_drafts/STATUS.md)，'
       '核同一传播、态及几何来源的实际映射；目标不改。')
"""+publish[b:]
publish=publish.replace('原物理Weyl观测与辅助测度的共同反射泛函','完整费米物质、局部时间来源与原记录的共同接口')
publish=publish.replace('原物理观测与辅助测度的共同内积','完整局部时间来源与原记录的共同域')
publish=publish.replace('原物理局部代数与共同背景必须同时保留','完整费米来源不能由纯玻色流代替')
write('publish_round662.py',publish)
write('postcheck_round662.py',replace((HERE/'postcheck_round661.py').read_text('utf8'),mapping))
write('round662_drafts/research_note_662_draft.md',(HERE/'research_note_662.md').read_text('utf8'))
review="""662 primary-agent analytic/code review; no independent agent review.
Previous goal response only reaffirmed integration priority: no science progress.
Current turn resumed actual662 entry, read all navigation,661 note/results,
585/587/588/589/598/601/623/624/650 and authoritative original655 code.
Detailed process query denied by OS; ordinary Get-Process returned no Python.
No unverified process was restarted. No goal/task/automation/image changes.
Selected nonnegative local allocations are explicit inputs, not uniqueness.
Full non-diagonal electric terms remain positive outgoing-link clusters;
different clusters use disjoint factors and strongly commute. Do not assume
each non-Abelian electric generator or shape-dependent block commutes.
Positive scalar potential Shubin theorem applied before matrix perturbation.
Lapse derivative of onsite potential has N_i^2, not N_i. Domains proved with
graph norms and complete factors; C1 history uses actual Kato hypotheses.
H[N] semibounded, not falsely positive before scalar shift. Fixed finite graph.
CAR mass is not globally bounded over configuration; original W^(1/4) controls.
Full bracket sign checked directly: [H_N,H_M]/i, three terms, unordered pairs.
Weak commutator on common D avoids assuming unbounded strong products exist.
Original local mass shares scalar-kinetic lapse, so its mixed terms cancel.
Record simplification restricted to inherited no-scalar-dependence hopping.
New CAR part commutes with original s instrument, but full dynamics, thermal
state, conditional outcomes and current variance may still change.
Instant variance-injection equality only on core/sufficient moment domain;
finite A^2 moment not asserted to imply finite squared current.
Numerics: actual64 modes, independent endpoint gauge transformation; neutral
four-mode16-Fock identity includes normal-order scalar. No full Gauss integral.
Electric test uses original physical hypercharges and exact finite Fourier
polynomials. Nambu coefficient vectors are not physical many-body states.
Both mixed and matrix omission produce nonzero errors; endpoint-half electric
allocation supplies zero-mixed-term positive control. No new current fitting.
Mature algebra/domain theorems are tools, not claimed general discoveries.
2 substantive groups,14 numbered equations. Entry repeat not extra test group.
Original585 ambiguity,588 closure restriction and612 spectrum remain intact.
Next663 targets real common propagation/source dictionary; no device design.
"""
for name in ('research_note_662.md','joint_local_lapse_fermion_current.py',
             'joint_local_lapse_fermion_current_results.json','unified_physics_condition_ledger_662.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round662_drafts/final_review.txt',review)
print('662 note/review frozen; verification and publication prepared')
