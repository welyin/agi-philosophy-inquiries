"""Freeze 652 common reference forms; preserve prior evidence and navigation."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        pat=re.escape(key)
        if key.isdecimal():pat=r'(?<!\d)'+pat+r'(?!\d)'
        pats.append(pat)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_material_boundary_transport':'joint_quantum_reference_forms',
    '650':'651','651':'652','652':'653','3188':'3191','3191':'3194',
    '1262':'1265','1265':'1268','2165':'2175','2175':'2186',
    'material_wall_joint_entry.md':'common_quantum_reference_entry.md',
    'material_wall_signature_probe':'curved_reference_velocity_probe'}
verify=replace((HERE/'verify_round651.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original causal material wall and full boundary-geometry transport.',
                      'Verify original quantum reference forms and common energy representation.')
verify=verify.replace("'round652_drafts/curved_reference_velocity_probe_results.json')",
                      "'round652_drafts/curved_reference_velocity_probe_results.json',\n           'round652_drafts/integration_before_design.md')")
verify=verify.replace("'unified_physics_condition_ledger_652.md','round653_drafts/STATUS.md'):",
                      "'unified_physics_condition_ledger_652.md','round653_drafts/STATUS.md',\n                 'round652_drafts/integration_before_design.md','round652_drafts/common_quantum_reference_entry.md'):")
write('verify_round652.py',verify)
publish=replace((HERE/'publish_round651.py').read_text('utf8'),mapping|{'## 297.':'## 298.','## 202.':'## 203.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第652轮完成：** [原物质量子参考、共同能量与区域表示的联合连接]({p}research_note_652.md)'
         '原光滑物质参考及速度平方共享自伴表示、原总能量矩控制和规范区域映射；'
         '同一谱子空间以正形式保留平方漏出项及几何偏导。'
         '三组、十六式通过，最新652／3194，1268份编号科学文件、2186份保护证据。'
         '[核验]({p}research_round_652_checks.json)、[条件账]({p}unified_physics_condition_ledger_652.md)。'
         '限原固定图共同表示；联合量子坐标、连续手征和动态引力仍开放。')
order=('**当前执行顺序（652后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[653完整图物质与连续手征分支]({p}round653_drafts/STATUS.md)，'
       '核同一状态、规范测度和动力学的真实连接，不继续参考装置或精度优化；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物质图册与共同边界输送','原量子参考与同一能量表示')
publish=publish.replace('物质标签须与因果边界共同相容','参考量子表示须与原过程共同相容')
publish=publish.replace('原物质参考、因果区域与完整边界变分的共同连接','原物质量子参考、共同能量与区域表示的联合连接')
write('publish_round652.py',publish)
write('postcheck_round652.py',replace((HERE/'postcheck_round651.py').read_text('utf8'),mapping))
write('round652_drafts/research_note_652_draft.md',(HERE/'research_note_652.md').read_text('utf8'))
review="""652 primary-agent review; no independent agent review.
Previous user-steering turn recorded execution order but no new science.
This goal turn revalidated current651 evidence and performed the next concrete
joint representation check; no scientific round credited to administrative work.
Read root README, research_direction, RESEARCH_STATE,651 note and receipt,
647/651 ledger,652 drafts and latest integration-before-design instruction.
CIM process inventory was denied; Get-Process fallback returned no Python
process. No prior live handle or652 completed artifacts existed on entry.
Read554,574-575,588,598,603,617,625,649. Searched prior notes for complete flows,
self-adjoint references and Galerkin/Ritz.588 already has cutoff flow tools;
560 already has square-before-compression. Neither is claimed as new theory.
Read primary Bandara-Saratchandran3.4. Our first-order field is not elliptic:
prove complete bounded flow directly, do not invoke their elliptic shortcut.
T=h^2/2 is575's old smooth material variable. Full Cartesian measure retained.
Direct divergence gives drift -F phi/(3M), Delta T=F(4-h^2/(2M)). Global
gradient/divergence bounds make the ORIGINAL fields complete, no cutoff added.
Flow parameter is a mathematical group, not added physical Hamiltonian.
Half-density generator has a smooth compact core, commutes Gauss; singleton
theta geometry is preserved, not duplicated. h=0 and F=0 receive no new wall.
Classical chart determinant h^3 only on old h>0 patch. No global time or
operator-valued inverse of classical coordinates is claimed.
Same original H (including all CAR potentials/hopping) determines V_f on core.
H+C+1 dominates kinetic plus identity by old full mass form bound. Hence
finite energy controls V^2 and first absolute Q moment; no higher moment
or universal continuum energy bound asserted.
W_f is an EXPLICIT bounded gauge-invariant finite-difference reference
prescription on original graph; not silently replacing geodesic H interaction.
Q=W-V^2 self-adjoint by bounded perturbation; R=b-Q positive form has D(V).
Energy eigenvectors need only lie in D(V), not assumed in D(V^2).
Union P_R is V-graph dense: smooth core in H form domain, then energy spectral
approximation and the proved form bound. Ritz resolvents converge by positive
variational projection. Finite ordered bounded words converge, not arbitrary
joint spectral probabilities. Same625 state map used, no rethermalization.
Square leak is positive and fixed by original velocity. It also contributes
to original geometric partial derivative. That derivative is not total thermal
response nor a new stress tensor. Source derivatives of arbitrary histories
are NOT inherited from zero-order strong convergence.
617/649 cut J is identity on matter/theta/CAR: intertwines flows, form domains,
same spectral projections. Matching-space unitary equivalence is not independent
regional preparation, ADM constraint solution or gravitational entropy.
Numerics: full-target divergence and chart Jacobian; same625 64D original
radial diagnostic spectrum, actual discrete commutator. Matrix tests ONLY
normal W=0 sector, no full Gauss spectrum, no infinite numerical error bound.
The finite ordered-word errors are not monotone; report actual rank24 increase.
Only variational energy error is proved/tested monotone. All original tests
passed on first run with stated tolerances. No image tests or goal changes.
16 display formulas; scope and shared/extra inputs audited. Next switch to
other cross-branch conditions, no detector/precision/architecture design.
"""
for name in ('research_note_652.md','joint_quantum_reference_forms.py',
             'joint_quantum_reference_forms_results.json','unified_physics_condition_ledger_652.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round652_drafts/final_review.txt',review)
print('652 artifacts frozen; verification/publication prepared')
