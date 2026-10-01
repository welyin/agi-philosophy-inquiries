"""Freeze655 full graph temporal process and geometric normalization interface."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        patterns.append(p)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_mass_state_time_limit':'joint_full_graph_transfer_sources',
    '653':'654','654':'655','655':'656','3198':'3202','3202':'3206',
    '1271':'1274','1274':'1277','2196':'2206','2206':'2218',
    'spatial_mass_common_entry.md':'normalization_geometry_entry.md',
    'neutral_mass_transfer_probe':'lapse_extension_probe'}
verify=replace((HERE/'verify_round654.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original mass, common state and auxiliary temporal limit.',
                     'Verify original full graph transfer, lapse sources and geometric normalization.')
verify=verify.replace("'round655_drafts/lapse_extension_probe_results.json')",
    "'round655_drafts/lapse_extension_probe_results.json',\n           'round655_drafts/round654_formula15_erratum.md','round655_drafts/development_diagnostic.md')")
verify=verify.replace("text['display_formulas']==20","text['display_formulas']==22")
verify=verify.replace("(?:qquad|quad)","(?:qquad|quad|longrightarrow|operatorname|stackrel)")
verify=verify.replace("    result=core.read(model.TARGET)","""    import importlib.util
    entry_spec=importlib.util.spec_from_file_location('entry655',HERE/'round655_drafts/lapse_extension_probe.py')
    entry=importlib.util.module_from_spec(entry_spec);entry_spec.loader.exec_module(entry)
    assert entry.run()==core.read(HERE/'round655_drafts/lapse_extension_probe_results.json')
    result=core.read(model.TARGET)""")
write('verify_round655.py',verify)
publish=replace((HERE/'publish_round654.py').read_text('utf8'),mapping|{'## 300.':'## 301.','## 205.':'## 206.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第655轮完成：** [原完整图的质量转移、共同时间来源与几何归一边界]({p}research_note_655.md)'
         '正质量分割接回原全图、Gauss与固定原态记录，正时间常lapse来源共同收敛。'
         '辅助欧氏投影不等于真实遗忘；固定计数归一不可直接认作协变真空体积项。'
         '四组、二十二式通过，最新655／3206，1277份编号科学文件、2218份保护证据。'
         '[核验]({p}research_round_655_checks.json)、[条件账]({p}unified_physics_condition_ledger_655.md)、'
         '[654式15勘误]({p}round655_drafts/round654_formula15_erratum.md)。'
         '限原有限图；实际空间手征测度、连续与量子引力仍开放。')
order=('**当前执行顺序（655后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[656实际空间overlap测度与正图过程]({p}round656_drafts/STATUS.md)，'
       '恢复原辅助Pfaffian的真实空间依赖，核共同映射；不继续单独优化时间步，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原质量态与共同时间极限','原全图过程与几何归一')
publish=publish.replace('质量、状态与辅助测度须共享同一极限','全图演化、时间来源与体积归一须共同签收')
publish=publish.replace('原完整质量、共同量子态与手征辅助测度的同一时间极限',
                        '原完整图的质量转移、共同时间来源与几何归一边界')
write('publish_round655.py',publish)
write('postcheck_round655.py',replace((HERE/'postcheck_round654.py').read_text('utf8'),mapping))
write('round655_drafts/research_note_655_draft.md',(HERE/'research_note_655.md').read_text('utf8'))
review="""655 primary-agent review; no independent agent review.
Previous goal turn was progress:654 completed/published/postchecked,655 actual
lapse-extension entry and reproducible probe saved. Current turn re-read root
and research READMEs/direction/state/latest654/ledger/results; reran654 verifier.
All2206 existing evidence hashes intact. No live Python process found.
No goal change, new thread, automation, subagent or image checks.
Read585 local lapse ambiguity,598 original graph/CAR/confinement,623 operator
core and domain,640 normalization limit,643 full heat source,645 volume/lapse,
646 positive-time/zero-time boundary. Searched prior notes for Chernoff/Vitali
and pair splitting; not repeating the585 general ambiguity as a new result.
Primary Chernoff reference and Arendt-Nikolski Theorem2.1 checked;3+1 source
definitions checked against Gourgoulhon. Existing theorems not claimed as new.
S_a retains original K,V,J and full original C. Full M_a, not its unbounded
creation factor alone, bounded by original quartic confinement. Complex-time
bound checked on Re z>0 compact sets. Tangency proof places M_a-I on the
original compact core, not on a heat-spread wavefunction.623 supplies the core.
Chernoff gives strong convergence. Banach-valued Vitali/Cauchy give positive-
time constant-lapse derivatives, not generic metric derivatives or zero-time
energy convergence. All constants may depend on fixed graph.
Original fixed trace-class rho used; no inference of Gibbs-trace convergence
from strong convergence alone. H_a=-logS_a/a well-defined since S_a has no
kernel; common lower bound and R_a=(I-S_a)/a yield strong-resolvent convergence.
Real-time bounded original record probabilities follow; realtime unbounded
source derivatives not asserted. Original joint Gauss action retained.
Auxiliary sitewise copy is an EXTRA candidate, not actual spatial overlap
Pfaffian. Euclidean excited suppression not unitary forgetting. Explicit two
subsequences show no general real-time auxiliary-state limit; this counterstate
has divergent energy and does not exclude uniformly finite-energy branches.
Raw count energy's lapse and psi derivatives compared using original volume
w=eps^3 psi^6. Matching scalar value at a reference geometry does not match
spatial derivative. No dust/GR interpretation inferred from only this block.
Volume-only scalar additivity fixes exponential form under declared restrictive
inputs, not Lambda coefficient. Subtract-count/add-volume is a stated matching
choice, not a prediction. A geometry-dependent scalar is not a global identity
when geometry is quantized. No hidden vacuum deletion.
Full64 original conditional spatial factor checked, actual complex masses,
nontrivial link, independent endpoint gauge changes. All five sectors retained.
Four-mode independent original neutral Fock comparison and original fixed-state
effect probe. First development run's invalid O(a^2) effect expectation failed;
parameters unchanged, corrected to actual O(a) convergence and preserved log.
654 formula15 missing backslash documented as frozen erratum; old hashes intact.
4 groups and22 formulas. Next656 returns to actual spatial auxiliary measure;
does not optimize time grid or design cognition devices.
"""
for name in ('research_note_655.md','joint_full_graph_transfer_sources.py',
             'joint_full_graph_transfer_sources_results.json','unified_physics_condition_ledger_655.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round655_drafts/final_review.txt',review)
print('655 frozen; verification/publication prepared')
