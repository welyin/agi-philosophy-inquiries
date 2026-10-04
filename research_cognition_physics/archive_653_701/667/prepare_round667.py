"""Prepare667, preserving prior numbered and development evidence."""
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


mapping={'joint_reference_mass_stabilizer':'joint_reference_spatial_process',
    '665':'666','666':'667','667':'668','3237':'3239','3239':'3241',
    '1307':'1310','1310':'1313','2316':'2329','2329':'2342'}
verify=replace((HERE/'verify_round666.py').read_text('utf8'),mapping)
verify=verify.replace('Verify reference/mass common stabilizer and transported actual auxiliary weights.',
    'Verify local reference transport, actual spatial masses and inherited geometry scope.')
verify=verify.replace('import joint_contact_gauss_history as prior_model','import joint_reference_mass_stabilizer as prior_model')
verify=verify.replace('particle_hole_contact_probe','reference_scale_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_667.md','round667_drafts/research_note_667_draft.md',
           'round667_drafts/final_review.txt','round668_drafts/STATUS.md',
           'round667_drafts/reference_scale_entry.md','round667_drafts/reference_scale_probe.py',
           'round667_drafts/reference_scale_probe_results.json',
           'round667_drafts/joint_reference_spatial_process_initial.py',
           'round667_drafts/joint_reference_spatial_process_initial_results.json',
           'round667_drafts/initial_result_retained_before_curved_check.json')
"""+verify[b:]
write('verify_round667.py',verify)
publish=replace((HERE/'publish_round666.py').read_text('utf8'),mapping|{'## 312.':'## 313.','## 217.':'## 218.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第667轮完成：** [局部参考、空间传播与共同量子过程的相容条件]({p}research_note_667.md)'
         '局部参考字典保区域和原过程，空间传播、质量及几何来源须共同输送；无限空Fock实现并非必要。'
         '两组、十六式通过，最新667／3241，1313份编号科学文件、2342份保护证据。'
         '[核验]({p}research_round_667_checks.json)、[全条件账]({p}unified_physics_condition_ledger_667.md)。'
         '明确复用382—386、425、522—523；旧条件性空间成果保留，实际共同尺度映射与量子GR仍开放。')
order=('**当前执行顺序（667后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[668原非零质量与手征辅助测度]({p}round668_drafts/STATUS.md)，'
       '核同一物理观测、权重、反射和来源，复用既有空间证明；目标不改。')
"""+publish[b:]
publish=publish.replace('原参考排序、质量与内部对称性的共同条件','局部参考、空间传播与共同量子过程的相容条件')
publish=publish.replace('原参考、质量与同一内部对称性','局部参考与原空间传播的共同连接')
publish=publish.replace('参考字典须与状态、权重和来源共同输送','复用旧维数与坐标证明，核实际物质和尺度映射')
write('publish_round667.py',publish)
write('postcheck_round667.py',replace((HERE/'postcheck_round666.py').read_text('utf8'),mapping))
write('round667_drafts/research_note_667_draft.md',(HERE/'research_note_667.md').read_text('utf8'))
review="""667 primary-agent proof/code audit; no independent agent review.
Inherited666 complete, actual667 entry saved. Latest user reminder keeps goal
and task unchanged. No new app task, automation, agents or visual checks.
Read382-384,386,425,522-523 reports/checks: qualitative384 lower bound no longer
requires extra Lipschitz/Hausdorff regularity.383 uses actual free involution.
386 finite strict margin does not require readout limit/absolute visibility.
386 strong dilatation and425 half/cost are alternative bridges, not cumulative.
425 needs no half homomorphism or global contraction.522 absolute all-pair
thermal criterion differs from local/angular accuracy.523 already jointly
realizes actual direction instrument and thermal dimension under its inputs.
Do not reopen these old proofs or count them as new; remaining task is their
same-object embedding in current matter/state/scale process. Existing algebraic
thermal states522-523 are acknowledged. Adding nodes is not spatial refinement.
Infinite empty Fock normality obstruction is proven by decreasing projections
and normality, not inferred from finite trace tests. Local CAR star automorphism
exists;16 transformed modes/site allow even implementers. Full original gauge
determinant character cancels. No ordinary tensor splitting of Gauss space.
Finite complete H domains/records/sources transported by constant unitary.
Infinite dynamics statement conditional on actual point-norm convergence in
appropriate observable algebra; not all B(H) dynamics or unbounded boson limit.
State and finite CP records transport; derivative convergence requires control.
Actual614 map used on all96/128 original modes, all original complex masses and
nontrivial SM links. Real scalar hopping becomes t Lambda R_left. Odd neutral
loop blocks the specific naive rewrite, not valid nonbipartite physical models.
Bipartite rephasing requires simultaneous Dirac sign changes; Majorana unchanged.
Actual574 curved target edge, not flat numerator, enters final boson check.
Initial code/results preserved as development evidence, not a failed experiment.
BdG thermal constant trace retained; an exact6-mode neutral invariant sector
independently checks it with64-dimensional Fock. Source finite difference checked.
Numerics are frozen-background quadratic tests, not full Gauss/boson integrals
or contact interacting trace. Geometry derivative restricted to stated family.
653 charge obstruction remains. Conditional certificates unchanged by alpha.
604 existing kinetic/doubling analysis is reused, not proposed as new668 work.
Next668 is actual nonzero mass in660-661 physical auxiliary observation dictionary.
Two fresh groups,16 equations. Four physics branches remain distinct; no final
dimension, unique SM, infinite dynamics or quantum GR derivation claimed.
"""
for name in ('research_note_667.md','joint_reference_spatial_process.py',
             'joint_reference_spatial_process_results.json','unified_physics_condition_ledger_667.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round667_drafts/final_review.txt',review)
print('667 report/review frozen; verification and publication prepared')
