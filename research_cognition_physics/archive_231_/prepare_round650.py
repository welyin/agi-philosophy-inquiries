"""Freeze650: actual constraints and common regional boundary flux."""
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
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_area_reduction_process':'joint_gravity_boundary_flux',
    '648':'649','649':'650','650':'651',
    '3182':'3185','3185':'3188','1256':'1259','1259':'1262',
    '2144':'2155','2155':'2165'}
verify=replace((HERE/'verify_round649.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original area quotient, process obstruction and common geometry extension.',
    'Verify original gravity constraint completion and common boundary flux.')
verify=verify.replace('import verify_round649 as previous\n','')
verify=verify.replace('base=previous.verify()',"base=core.read(HERE/'research_round_649_checks.json')\n    assert base['all_reported_checks_passed']")
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_650.md','round650_drafts/research_note_650_draft.md',
           'round650_drafts/final_review.txt','round651_drafts/STATUS.md',
           'round650_drafts/boundary_flux_entry.md','round650_drafts/original_boundary_flux_probe.py',
           'round650_drafts/original_boundary_flux_probe_results.json')
"""+verify[end:]
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==16")
verify=verify.replace('previous_results_unchanged=True,full_historical_science_rerun=False,',
    'previous_results_unchanged=True,full_historical_science_rerun=False,\n        previous_round_receipt_and_all_evidence_hashes_checked=True,')
write('verify_round650.py',verify)
publish=replace((HERE/'publish_round649.py').read_text('utf8'),mapping|{'## 295.':'## 296.','## 200.':'## 201.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第650轮完成：** [原引力约束、面积匹配与物质边界通量的联合条件]({p}research_note_650.md)'
         '原物质模有满足全部经典初始约束的小数据解族；两侧面积匹配仍有单侧通量，'
         '完整物质—几何动量拼接使两侧通量相消。'
         '三组、十六式通过，最新650／3188，1262份编号科学文件、2165份保护证据。'
         '[核验]({p}research_round_650_checks.json)、[条件账]({p}unified_physics_condition_ledger_650.md)。'
         '限给定领先Einstein作用、固定类时切口与经典分支；关系区域及全量子统一仍开放。')
order=('**当前执行顺序（650后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[651物质关系区域、动态嵌入与共同边界配对]({p}round651_drafts/STATUS.md)，'
       '核原参考的实际输送及完整边界条件，不继续面积或控制器参数扫描；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('面积约化与同一量子几何的共同过程','面积匹配与实际物质—引力通量')
publish=publish.replace('正物理内积还须容纳原演化','区域演化还须保持完整边界配对')
publish=publish.replace('面积匹配的物理内积、动力学下降与共同量子几何','原引力约束、面积匹配与物质边界通量的联合条件')
write('publish_round650.py',publish)
write('round650_drafts/research_note_650_draft.md',(HERE/'research_note_650.md').read_text('utf8'))
review="""650 primary-agent review; no independent agent review.
Latest649 navigation, completed science and2155 evidence files checked.
The immediately previous explanatory turn was no progress, not a blocking
condition:650 unfinished probe and a concrete next proof remained available.
No live Python process was found. No goal or automation changed.
Read573,600,618,619,647-649 and650 saved entry/probe/results; searched old notes
for AF completion and boundary symplectic flux to avoid duplicate claims.
Iyer-Wald78-82 and Harlow-Wu1-6 and3.2 read. Their tools are inputs; no novel
general covariant-phase-space theorem claimed.
The same original U radial Hessian and K produce the two physical mass modes.
Zero gauge and fermions with fixed real Higgs radial orientation form a
consistent classical sector; no linear or nonlinear angular gauge current.
The AF compact initial scalar data are followed by an actual constraint solve:
Newton fixed point of -8 Delta psi=B psi+2U psi^5, p=0 and spatial k=0.
The continuous kernel norm R^2/16 and global B,U bounds prove contraction,
positivity and local uniqueness. Numerical samples are not sup certificates.
Smooth parameter dependence follows invertibility of I-DT; all linear metric
perturbations vanish because the original vacuum stress starts quadratically.
Local hyperbolic development reuses573's principal-symbol and constraint
propagation mapping. No future completeness or T3 linearization claim.
The cutoff AF waves are not assigned the uncut plane-wave time formula.
At t=0 their compact normal flux is positive and integrates to0.609908597826.
Each side area only has zero FIRST variation. Exact equality is between two
restrictions of the SAME solution; no exactly fixed-area solution family is
claimed. This suffices against area-matching-alone implies single-region
zero-flux. It does not refute other allowed boundary conditions.
At the vacuum both metric perturbations vanish. Full Jordan normal curvature
variation is retained, producing precisely the same scalar pairing. Freezing
it is an expressly inconsistent comparison, not a competing physical frame.
Field-independent xi, fixed region, standard unextended Omega and unrestricted
specified variations are part of the integrability test. Exterior derivative
of alpha cancels mixed delta Q. Nonzero curl excludes only this Hamiltonian
one-form; boundary symplectic extensions/flux charges remain possible.
The same leading action's full momenta match on the transmission subspace.
Taking exterior derivative there cancels flux; standalone regions need not
be autonomous. General-jet FD checks are kinematic, AF proof is on shell.
This background's original four material coordinates are degenerate; no
relation-region result or full quantum gravity is inferred.651 addresses
that genuine object mismatch. No hardware or resource optimization.
Three meaningful numerical groups,16 formulas. A missing qquad slash was
fixed before freezing. No visual checks. Old frozen evidence retained.
Previous receipts plus all historical hashes checked; only650 and its actual
imported/probe dependencies rerun, not every old numerical experiment.
"""
for name in ('research_note_650.md','joint_gravity_boundary_flux.py',
             'joint_gravity_boundary_flux_results.json','unified_physics_condition_ledger_650.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round650_drafts/final_review.txt',review)
print('650 frozen; verifier and publisher prepared')
