"""Freeze the actual boundary mapping; preserve all previous evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    parts=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        parts.append(p)
    return re.sub('|'.join(parts),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_auxiliary_normalization_support':'joint_domain_wall_boundary_mapping',
    '657':'658','658':'659','659':'660','3214':'3217','3217':'3221',
    '1283':'1286','1286':'1289','2239':'2246','2246':'2256'}
verify=replace((HERE/'verify_round658.py').read_text('utf8'),mapping)
verify=verify.replace('Verify exact zero-sector support and strict original auxiliary normalization.',
    'Verify actual domain-wall boundary, signed weight, limits and source normalization.')
verify=verify.replace('(3,0,0)','(4,0,0)').replace('run=3','run=4')
verify=verify.replace("    result=core.read(model.TARGET)","""    import importlib.util
    probe_path=HERE/'round659_drafts/temporal_memory_probe.py'
    spec=importlib.util.spec_from_file_location('round659_probe',probe_path)
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    assert probe.run()==core.read(probe_path.with_name('temporal_memory_probe_results.json'))
    result=core.read(model.TARGET)""")
verify=verify.replace("'round659_drafts/final_review.txt','round660_drafts/STATUS.md')",
    """'round659_drafts/final_review.txt','round660_drafts/STATUS.md',
           'round659_drafts/temporal_process_entry.md','round659_drafts/temporal_memory_probe.py',
           'round659_drafts/temporal_memory_probe_results.json')""")
write('verify_round659.py',verify)
publish=replace((HERE/'publish_round658.py').read_text('utf8'),mapping|{'## 304.':'## 305.','## 209.':'## 210.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第659轮完成：** [原手征权重、局部边界表示与共同来源的连接]({p}research_note_659.md)'
         '原带相位辅助权重接到指定domain-wall真实边界；受控双极限保留整个球面平均。'
         '边界变换与体减除来源不可另选或删去。'
         '四组、十六式通过，最新659／3221，1289份编号科学文件、2256份保护证据。'
         '[核验]({p}research_round_659_checks.json)、[全条件账]({p}unified_physics_condition_ledger_659.md)。'
         '共同物理时间、原CAR、一般规范场及量子GR仍开放。')
order=('**当前执行顺序（659后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[660原测度与共同物理时间]({p}round660_drafts/STATUS.md)，'
       '核成熟反射／Hamiltonian结果的实际映射；不优化第五方向层数或装置，目标不改。')
"""+publish[b:]
publish=publish.replace('原球面平均的零模支撑与严格归一','原手征权重、局部边界表示与共同来源的连接')
publish=publish.replace('原球面平均的零模支撑','原测度的真实边界连接')
publish=publish.replace('原平均、端点支撑与归一反射内积共同成立','原测度、边界与来源须共用同一映射')
write('publish_round659.py',publish)
write('postcheck_round659.py',replace((HERE/'postcheck_round658.py').read_text('utf8'),mapping))
write('round659_drafts/research_note_659_draft.md',(HERE/'research_note_659.md').read_text('utf8'))
review="""659 primary-agent review only; no independent-agent review.
Previous turn was orientation/priority discussion, not scientific progress.
Current turn reads authoritative658, full ledger647 and659 actual entry.
No live Python process on entry; only current verified sessions used.
No goal edits, new tasks, automations or image checks.
Original finite free Wilson X, B, all16 T matrices and S9 kept.
Mature source is Kikukawa1710.11618 section6.4 eq234, not its eq239.
Modified last diagonal Qplus A removes precisely last Qminus bar rows.
Finite chiral transfer follows original Dprime zero equations, u0=0.
T positive by block factorization, but this is auxiliary fifth direction.
(2+aX)M=I+T and gamma5 aX M=I-T fix H_a and actual boundary map.
R_infinity=(2+aX)^-1 range(Pnegative(H_a)); dropping map is incorrect.
Even AP free physical time ensures sine vector nonzero and invertible
initial chiral overlap, for each fixed finite box. No uniform volume claim.
Finite filter error uses two spectral blocks and cond(M) K q^(L-1).
Scaled H_a perturbation bound and resolvent sign projector bound valid only
when delta<gap. Field-map projector gap <=a||H||/2, code uses weaker bound.
Both a->0 and aL->infinity needed for inherited P; fixed a has bias.
Free antiunitary B symmetry keeps reference restricted pairing unitary.
Grassmann saturation: nonzero bar modes force all nonzero psi directions;
only right null pairing survives. Boundary frame Jacobians cancel E ratio.
Full1024 Nambu Pfaffian check retains original internal16 and signed phase.
Fixed-volume Pfaffian contraction derivative bound gives uniform S9 limit;
658 Z>0 permits normalized bounded auxiliary insertions. No finite-a state
positivity or general geometric derivative-limit exchange claimed.
Reference/bulk determinants retained for source kappa: nonzero residual
source verified, rather than deleting a source-dependent normalization.
Analytic result applies to fixed finite free branch, not interacting gauge,
physical32CAR Hamiltonian, common positive physical time or quantum GR.
All27 conditions refreshed; different model branches still distinct.
612/613 already cited2026 Kikukawa abstract; next660 reads full material
and must not repeat old CHN naive axial-charge filtering counterexample.
4 substantive groups,16 equations. Entry probe reproduced but not counted
as an extra round/group. Historical incomplete entry retained as history.
"""
for name in ('research_note_659.md','joint_domain_wall_boundary_mapping.py',
             'joint_domain_wall_boundary_mapping_results.json','unified_physics_condition_ledger_659.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round659_drafts/final_review.txt',review)
print('659 report and code frozen; verification/publication prepared')
