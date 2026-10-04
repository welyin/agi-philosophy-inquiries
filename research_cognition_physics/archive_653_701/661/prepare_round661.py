"""Freeze original physical/auxiliary joint interface, retaining all history."""
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


mapping={'joint_local_mirror_process':'joint_physical_auxiliary_state',
    '659':'660','660':'661','661':'662','3221':'3225','3225':'3227',
    '1289':'1292','1292':'1295','2256':'2269','2269':'2279'}
verify=replace((HERE/'verify_round660.py').read_text('utf8'),mapping)
verify=verify.replace('Verify full onsite mirror completion, original observations and source/RP interface.',
    'Verify original physical/auxiliary reflection state and common source.')
verify=verify.replace('mirror_schur_probe','observable_interface_probe')
verify=verify.replace('(4,0,0)','(2,0,0)').replace('run=4','run=2')
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==14")
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_661.md','round661_drafts/research_note_661_draft.md',
           'round661_drafts/final_review.txt','round662_drafts/STATUS.md',
           'round661_drafts/observable_interface_entry.md','round661_drafts/observable_interface_probe.py',
           'round661_drafts/observable_interface_probe_results.json')
"""+verify[b:]
write('verify_round661.py',verify)
publish=replace((HERE/'publish_round660.py').read_text('utf8'),mapping|{'## 306.':'## 307.','## 211.':'## 212.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第661轮完成：** [原物理Weyl观测与辅助测度的共同反射泛函]({p}research_note_661.md)'
         '原局部物理观测与E测度共用有限自由归一反射内积；不必先局部实现全部镜像变量。'
         '共同背景下的联合观测响应必须同时保留两部门。'
         '两组、十四式通过，最新661／3227，1295份编号科学文件、2279份保护证据。'
         '[核验]({p}research_round_661_checks.json)、[条件账]({p}unified_physics_condition_ledger_661.md)。'
         '原标量记录、共同Hamiltonian、一般规范场及量子GR仍开放。')
order=('**当前执行顺序（661后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[662共同局部来源、关系区域与引力约束]({p}round662_drafts/STATUS.md)，'
       '回查原lapse歧义和已有来源，检验实际约束接口；目标不改。')
"""+publish[b:]
publish=publish.replace('原手征权重、局部配对与观测变换的共同表示','原物理Weyl观测与辅助测度的共同反射泛函')
publish=publish.replace('原权重与局部正表示的观测接口','原物理观测与辅助测度的共同内积')
publish=publish.replace('权重相同仍须保留观测与来源映射','原物理局部代数与共同背景必须同时保留')
write('publish_round661.py',publish)
write('postcheck_round661.py',replace((HERE/'postcheck_round660.py').read_text('utf8'),mapping))
write('round661_drafts/research_note_661_draft.md',(HERE/'research_note_661.md').read_text('utf8'))
review="""661 primary-agent review only; no independent agent review.
Previous goal turn completed660 and saved actual661 probe: progress.
Current navigation,660 note/results,661 entry and original624/625/647 read.
No live Python process at entry. Missing guessed657 filename resolved through
actual file search; no reliance on that failed read. No goals or tasks changed.
Mature2010 Weyl RP applies only to specified local spin components.
Original index-zero free integral factors after saturating mirror variables.
All departments share D, time reflection, cut and source. No new fitted state.
Schur product positivity is inherited, not a novel general theorem.
Auxiliary bosonic variables commute, so no extra graded tensor-product sign.
Positive half chosen t>=L for fermions; swapping657 scalar halves conjugates
its positive Gram, preserving positivity. Original normalization uses658.
OS Hilbert tensor identity is per fixed finite box, not an all-time H theorem.
Analytic support obstruction: nx1 nt2 D10=-gamma4/2, and mapped local w is
Jminus^dagger(I-D)psi, giving exact cross-half operator norm1/2.
Entry nx2 nt4 matrix mismatch is conditional on E; averaging without sources
may remove particular anomalous pairs. No unproved average difference claim.
Full S9xS9 check uses radial Beta9/2 weight and original2-time Pfaffian,
closed-chain exponent16 (not8). Z=(9/2)_16/(9)_16, feature eigenvalue16/25.
704-dimensional joint Gram retains all physical16 internal indices.
Spatial finite kernel sample does not replace full sphere integration proof.
1024 Nambu decomposable bilinear source Pfaffian checks original complex phase.
Common-source holonomy uses original16 representation charge-4. Initial code
mistakenly looked for vector charge-2, failed before results, was corrected
from authoritative charge list, and complete run then passed. Not a theorem.
Seven-point radial Legendre rule includes Beta3,2 density and r insertion;
background source follows physical Weyl insertion plus auxiliary normalized
mean. Free RP is not inferred for the charged holonomy diagnostic.
Original scalar s instrument, Gauss/ADM constraints, general gauge, same Gibbs,
infinite time and continuum remain open. E never identified with physical s.
2 substantive groups,14 equations. Entry probe rerun not another test group.
Return next to common local sources,585 lapse gap and existing regions; no
further cognitive device design, target reduction or arbitrary clock addition.
"""
for name in ('research_note_661.md','joint_physical_auxiliary_state.py',
             'joint_physical_auxiliary_state_results.json','unified_physics_condition_ledger_661.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round661_drafts/final_review.txt',review)
print('661 note and review frozen; verification/publication prepared')
