"""Freeze 633: original continuous records, sources, and causal boundary."""
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
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_background_contact_matching':'joint_continuum_record_sources',
         '631':'632','632':'633','633':'634',
         '3133':'3136','3136':'3140','1205':'1208','1208':'1211',
         '2015':'2022','2022':'2029'}
verify=replace((HERE/'verify_round632.py').read_text('utf8'),mapping)
verify=verify.replace('Verify shared vacuum potential, local Weyl contacts and background conditions.',
                      'Verify original continuous record states, sources and causal counterexample.')
verify=verify.replace("['display_formulas']==16","['display_formulas']==18")
verify=verify.replace("==(3,0,0)","==(4,0,0)").replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=4,')
write('verify_round633.py',verify)
publish=replace((HERE/'publish_round632.py').read_text('utf8'),
                mapping|{'## 278.':'## 279.','## 183.':'## 184.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第633轮完成：** [连续记录、态正则性与同一物质的几何来源]({p}research_note_633.md)'
         '原混合中性场的光滑记录后态保短距结构与有限能源，并共同给质量和应力来源；'
         '同一场两分离模排除瞬时分布式读取。'
         '四组、十八式通过，最新633／3140，1211份编号科学文件、2029份保护证据。'
         '[核验]({p}research_round_633_checks.json)、[条件账]({p}unified_physics_condition_ledger_633.md)。'
         '限自由固定背景；真实仪器、图映射和GR仍开放。')
order=('**当前执行顺序（633后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[634记录条件化、来源与几何约束]({p}round634_drafts/STATUS.md)，'
       '核同一过程的状态边界与守恒；不继续波包或装置优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原真空势与几何接触项的共同匹配','原连续记录与物质几何来源的共同连接')
publish=publish.replace('共同背景匹配不选择空间维数','记录条件的整合不选择空间维数')
publish=publish.replace('同一真空势、局部接触项与背景方程的共同匹配','连续记录、态正则性与同一物质的几何来源')
write('publish_round633.py',publish)
write('round633_drafts/research_note_633_draft.md',(HERE/'research_note_633.md').read_text('utf8'))
review="""633 primary-agent review. No independent agent review.
Previous conversational goal turn was no progress: it restated priorities.
This round advances the frozen 633 entry; no prior science/status file changed.
Latest README, direction, state, 632 note and 633 entry were reread.
Old Python processes remain live; no current Oct-1 research worker was seen.
No old process was terminated or assumed to have failed.
524 FV scalar probe, 592 finite-graph record regularity and 624 conditional
source response were reviewed. Their general results are inherited.
Sanders 0911.1304 propositions 4.8 and 4.12 were read from the primary paper;
FV 1810.06512 section 3 was checked for causal process vs algebra locality.
Hadamard closure and CAR/CP identities are mature tools, not new general claims.
Original Y, F and 632 tree background are retained; the sterile field is not
replaced by a mass eigenfield or an independent massless probe.
The neutral invariant sector retains Dirac and Majorana mixing and both spins.
All charged species remain spectators in the fixed free continuum branch.
Nambu c-dagger(-k) requires +K(k).T in the lower block. Particle-hole,
anticommutator and full neutral mass dispersions are independently checked.
The original field spectral weights sum to one; no Dirac/Majorana counting
factor is appended to the expectation for one normalized measured mode.
Smooth compact Cauchy data guarantee all weighted momentum moments;
finite polynomial insertions give a smooth finite-rank two-point difference.
This proves Hadamard and finite energy for these particular operations,
not for arbitrary bounded local algebra elements or arbitrary QFT states.
Finite Fock checks test covariance/source normalization, not the continuum
regularity theorem. The continuum packet is truly compact, not Gaussian.
Vacuum energy uses a fixed renormalization; only state differences are finite.
Record-dependent resetting of 632 vacuum matching is explicitly disallowed.
Conditional energy formula uses Hvac Omega=0; it is not copied onto arbitrary
conditional mass or stress observables. Nonselective shared source identities
are derived separately from the same covariance difference.
Pressure denotes one third of integrated spatial stress trace, not necessarily
an isotropic local fluid or a solved metric. All values are diagnostic units.
Numerical Fourier tail is bounded analytically by an H2 norm; its coefficient
was evaluated numerically, not interval-certified. Grid agreement is not a
certified total error. No numerical cutoff is used in the existence proof.
Causal witness uses the original sterile field's two separated smooth modes
and a different initial state, not the vacuum injection sample.
A nonzero vacuum empty-mode projection plus creation polynomial supplies
a Hadamard finite-energy extension, justified by strict Fourier Schwarz.
It excludes instantaneous whole-region distributed implementation, not every
finite-duration local measurement or all possible cognition architectures.
Algebraic locality, energy regularity, and physical implementation are not conflated.
632 local coefficients, dimension, background, Y and physical apparatus remain
inputs/open conditions. No full Gauss continuum mapping, interacting SM,
dynamical metric, or independent GR derivation is claimed.
Four groups passed first execution. Eighteen formulas checked without images.
Next 634 integrates state boundaries, record conditioning and common conservation,
not deeper record hardware design. Broad active goal remains unchanged.
"""
for name in ('research_note_633.md','joint_continuum_record_sources.py',
             'joint_continuum_record_sources_results.json','unified_physics_condition_ledger_633.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round633_drafts/final_review.txt',review)
print('633 science frozen; verifier and publisher prepared')

