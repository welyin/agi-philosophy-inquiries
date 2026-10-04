"""Freeze the actual weight/observable/RP interface and prepare publication."""
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


mapping={'joint_domain_wall_boundary_mapping':'joint_local_mirror_process',
    '658':'659','659':'660','660':'661','3217':'3221','3221':'3225',
    '1286':'1289','1289':'1292','2246':'2256','2256':'2269'}
verify=replace((HERE/'verify_round659.py').read_text('utf8'),mapping)
verify=verify.replace('Verify actual domain-wall boundary, signed weight, limits and source normalization.',
    'Verify full onsite mirror completion, original observations and source/RP interface.')
verify=verify.replace('temporal_memory_probe','mirror_schur_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_660.md','round660_drafts/research_note_660_draft.md',
           'round660_drafts/final_review.txt','round661_drafts/STATUS.md',
           'round660_drafts/physical_process_entry.md','round660_drafts/mirror_schur_probe.py',
           'round660_drafts/mirror_schur_probe_results.json','round660_drafts/source_receipt.json',
           'round660_drafts/talk_20260629_kikukawa.pdf','round660_drafts/talk_20260629_kikukawa.txt')
    receipt=core.read(HERE/'round660_drafts/source_receipt.json')
    assert core.digest(HERE/'round660_drafts/talk_20260629_kikukawa.pdf')==receipt['sha256']
    assert core.digest(HERE/'round660_drafts/talk_20260629_kikukawa.txt')==receipt['text_sha256']
"""+verify[b:]
write('verify_round660.py',verify)
publish=replace((HERE/'publish_round659.py').read_text('utf8'),mapping|{'## 305.':'## 306.','## 210.':'## 211.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第660轮完成：** [原手征权重、局部配对与观测变换的共同表示]({p}research_note_660.md)'
         '同一局部配对精确保留原带相位权重及变换后的观测；自由有限盒候选的反射正性与严格归一成立。'
         '原观测的时间支撑尚未接通，帧Jacobian及来源必须同时输送。'
         '四组、十六式通过，最新660／3225，1292份编号科学文件、2269份保护证据。'
         '[核验]({p}research_round_660_checks.json)、[条件账]({p}unified_physics_condition_ledger_660.md)。'
         '原共同Hamiltonian、一般规范场及量子GR仍开放。')
order=('**当前执行顺序（660后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[661原观测代数与共同物理时间]({p}round661_drafts/STATUS.md)，'
       '核实际观测的时间支撑、反射及演化映射，限定接口后回到其它共同条件；目标不改。')
"""+publish[b:]
publish=publish.replace('原手征权重、局部边界表示与共同来源的连接','原手征权重、局部配对与观测变换的共同表示')
publish=publish.replace('原测度的真实边界连接','原权重与局部正表示的观测接口')
publish=publish.replace('原测度、边界与来源须共用同一映射','权重相同仍须保留观测与来源映射')
write('publish_round660.py',publish)
write('postcheck_round660.py',replace((HERE/'postcheck_round659.py').read_text('utf8'),mapping))
write('round660_drafts/research_note_660_draft.md',(HERE/'research_note_660.md').read_text('utf8'))
review="""660 primary-agent review only; no independent-agent review.
Previous goal turn was an ordering/status answer, not scientific progress.
Current turn reread current659 navigation, ledger, results,660 actual entry.
Get-CimInstance was unavailable; fallback Get-Process showed no Python.
Only actual current command handles used. No new app tasks, goal changes,
automations, image inspection or independent-agent delegation.
Original D,B,16 T,S9, four-dimensional lattice are inherited inputs.
New completion uses common E in the two onsite pairings, exp(+xi N xi/2).
2017 projected mirror and unprojected eq206 are distinct; its pair coefficient
has no1/2, so numeric couplings cannot be copied without conversion.
All b,c,d,e mixed blocks retained before exact determinant-one shear.
Sign in eprime contains minus half Kl^-T F c; all cross blocks checked.
Pfaffian phase and detS kept; no positive square root substituted.
Weight identity can extend to Kl zeros within same index-zero patch;
inverse source map cannot. Other topological index sectors not covered.
Original bar-E observable algebra is not identical to common-E completion.
Bare e-e mismatch is conditional on E; no claim of a nonzero no-source
S9-averaged pair. E-dependent test observations witness full algebra mismatch.
Wick kernel is -N^-1 for positive exponent; source transport fixes two/four
point insertions. Frame-rephasing test generates spurious18i if detS dropped.
Full holonomy integral uses exact degree10 radial quadrature, not Monte Carlo.
RP relies on mature free-overlap positive cross decomposition, not symmetry
or sampled eigenvalues. Onsite even Vplus+theta(Vplus), compact positive S9
product preserve the cone. Candidate action is not literally2010 eq51.
Strict normalization: constant E0 pairing scaled lambda, compressed positive
lambda Qplus+lambda^-1 Qminus is invertible via free antiunitary B symmetry;
Pf remains real and nonzero from positive free detD. Then658 average applies.
Free-RP input restricted to unit gauge background and physical even AP time.
No RP assertion for tested charged holonomy, no general gauge reconstruction.
Positive-half support of dressed observations is a separate needed map.
Original655 CAR/H,612 spectrum restriction, general gauge and continuum GR
remain open. No fifth-layer limit required here;659 bulk source stays intact
for its separate representation. Not a derivation of3+1, group or couplings.
4 substantive test groups,16 equations. Initial probe rerun not a new group.
661 entry records a specific observable-interface test, not cognitive design.
All frozen historical science must be hash-checked before publication.
"""
for name in ('research_note_660.md','joint_local_mirror_process.py',
             'joint_local_mirror_process_results.json','unified_physics_condition_ledger_660.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round660_drafts/final_review.txt',review)
print('660 report/code review frozen; verification/publication prepared')
