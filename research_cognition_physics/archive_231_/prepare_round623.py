"""Freeze 623 science and prepare history-preserving publication."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        pat=re.escape(key)
        if key.isdecimal():pat=r'(?<!\d)'+pat+r'(?!\d)'
        patterns.append(pat)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as stream:stream.write(text)


mapping={'joint_causal_source_functional':'joint_operator_domain_completion',
         '621':'622','622':'623','623':'624',
         '3106':'3109','3109':'3112','1175':'1178','1178':'1181',
         '1944':'1951','1951':'1958'}
verify=replace((HERE/'verify_round622.py').read_text('utf8'),mapping)
verify=verify.replace('Verify common quadratic CTP sources and causal coefficient limits.',
                     'Verify original operator domain and actual thermal process expansion.')
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==18")
write('verify_round623.py',verify)
publish=replace((HERE/'publish_round622.py').read_text('utf8'),
                mapping|{'## 268.':'## 269.','## 173.':'## 174.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第623轮完成：** [原完整量子过程的共同算符域与真实热响应]({p}research_note_623.md)'
         '原束缚势共同控制动能、边耦合与CAR质量；指定参数族在同一Gauss热态中具有实际二阶过程展开，'
         '接通622的噪声、迟致及接触项。'
         '三组、十八式通过，最新623／3112，1181份编号科学文件、1958份保护证据。'
         '[核验]({p}research_round_623_checks.json)、[条件账]({p}unified_physics_condition_ledger_623.md)。'
         '限固定图与正外部几何；连续、动态引力与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（623后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[624原实际记录与共同量子来源的接口]({p}round624_drafts/STATUS.md)，'
       '先回查原记录及尺度合同，再补真实共同过程的缺口；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('共同二阶过程、时间排序与虚跃迁','共同算符域、原物质与真实热过程')
publish=publish.replace('二阶因果匹配不构成维数选择','原完整过程展开不构成维数选择')
publish=publish.replace('同一量子过程的有限窗噪声、迟致核与接触项','原完整量子过程的共同算符域与真实热响应')
assert 'hashlib.sha256' in publish
write('publish_round623.py',publish)
write('round623_drafts/research_note_623_draft.md',(HERE/'research_note_623.md').read_text('utf8'))
review='''623 primary-agent review. No independent agent review.
622 is verified scientific progress; the intervening conversational priority
answer was not a new scientific round. Current turn revalidated files and
took the available next research action without changing the active goal.
Read current README/direction/state/622 note/results/623 entry and original
574/589/598/603 definitions. CIM process query was denied; Get-Process
fallback confirmed no Python processes started on or after 2026-10-01.
Older visible processes were left untouched; no new app task or timer.
Original polynomial U is exact in hyperboloid coordinates. Global coercive
and Laplacian estimates are proved, not inferred from 25 sampled points.
Configuration dimension five is not a physical spacetime assertion.
Shubin's scalar theorem applies only to T+W on a complete manifold.
Integration-by-parts plus closedness proves operator-domain intersection;
common form comparability alone is not used to assert operator domains.
Central compact Casimir commutes with all original invariant momenta,
including direction-mixed kinetic terms; node Laplacians have their joint
spectral decomposition. Derivative bounds follow on the same fibers.
Original geodesic edges are subordinate by triangle inequality through
the fixed target origin. Full CAR mass bound includes complex Dirac and
Majorana pairing. Magnetic potential and fixed-graph hopping are bounded.
The perturbation is infinitesimally operator-bounded. Resolvent argument
gives the same self-adjoint operator as the existing form construction.
Parameters keep target M,L,u fixed, all geometry uniformly positive, and
vary original gauge-invariant coefficients smoothly; no dynamic geometry.
Kato common-domain hypotheses checked against Schmid-Griesemer primary
paper. Duhamel used on D vectors, with strong derivative on one side and
backward derivative on the other for weak second Peano expansion.
No common D(H^2), source essential self-adjointness, strong second U
derivative, or all-source Frechet C2 conclusion is claimed.
Thermal trace exchange uses the summable original A2 moment and an
explicit pairwise graph-norm bound, not the convergence of coefficients
alone. Each branch shares the same initial state. Hessian contact retained.
Instantaneous correlations are HS Gram pairings; unproved products of
unbounded operators are not used. Abstract 621/622 counterexamples remain.
Three numerical groups passed on their first run. Finite-radius boundary
is nonzero and retained; numerical factor checks do not replace global
proof or compute the full Gauss spectrum. No numerical result was retuned.
Proof claims fixed-history scalar second expansion; ordinary strong C2
propagator and graph-uniform or UV claims are explicitly withheld.
Next: audit original record insertion and scale interfaces, no new cognitive
hardware. All original science and frozen evidence preserved.
'''
for name in ('research_note_623.md','joint_operator_domain_completion.py',
             'joint_operator_domain_completion_results.json','unified_physics_condition_ledger_623.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round623_drafts/final_review.txt',review)
print('623 verifier, publication script, frozen draft and review prepared')
