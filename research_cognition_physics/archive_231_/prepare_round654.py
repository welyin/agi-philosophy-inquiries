"""Freeze 654 joint mass/state/time limit and prepare reviewed publication."""
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


mapping={'joint_chiral_character_state':'joint_mass_state_time_limit',
    '652':'653','653':'654','654':'655','3194':'3198','3198':'3202',
    '1268':'1271','1271':'1274','2186':'2196','2196':'2206',
    'common_chiral_state_entry.md':'spatial_mass_common_entry.md',
    'holonomy_charge_probe':'neutral_mass_transfer_probe'}
verify=replace((HERE/'verify_round653.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original subgroup auxiliary character and actual temporal transfer.',
                      'Verify original mass, common state and auxiliary temporal limit.')
write('verify_round654.py',verify)
publish=replace((HERE/'publish_round653.py').read_text('utf8'),mapping|{'## 299.':'## 300.','## 204.':'## 205.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第654轮完成：** [原完整质量、共同量子态与手征辅助测度的同一时间极限]({p}research_note_654.md)'
         '原全质量的候选Pfaffian给同一正转移、态及来源；有限步状态有一阶修正。'
         '原固定辅助核在共同连续时间极限只留下规范单态，其余能隙发散。'
         '四组、二十式通过，最新654／3202，1274份编号科学文件、2206份保护证据。'
         '[核验]({p}research_round_654_checks.json)、[条件账]({p}unified_physics_condition_ledger_654.md)。'
         '限空间平凡及紧标量背景；完整手征、空间与引力仍开放。')
order=('**当前执行顺序（654后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[655原时间归一、几何来源与空间传播]({p}round655_drafts/STATUS.md)，'
       '审查原真空标量与同一几何来源，接回空间过程；不继续单独优化时间步，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原辅助测度与规范时间拼接','原质量态与共同时间极限')
publish=publish.replace('手征状态须与测度及时间拼接共同相容','质量、状态与辅助测度须共享同一极限')
publish=publish.replace('原手征辅助测度、规范状态与时间拼接的共同表示','原完整质量、共同量子态与手征辅助测度的同一时间极限')
write('publish_round654.py',publish)
write('postcheck_round654.py',replace((HERE/'postcheck_round653.py').read_text('utf8'),mapping))
write('round654_drafts/research_note_654_draft.md',(HERE/'research_note_654.md').read_text('utf8'))
review="""654 primary-agent review; no independent agent review.
User's priority is integration of remaining conditions before cognitive design.
No goal edit, new app task, automation or subagent. No image checks.
Read latest navigation/state/direction,653 ledger/receipt and654 pending entry.
No live Python process found. Existing654 STATUS is frozen under653 and kept.
Reused actual598 mass coefficients,614 fixed particle-hole dictionary,643 CAR
functions,653 frozen auxiliary spectrum. Read623 full operator-domain scope.
Searched numbered notes for fixed auxiliary ratio, asinh and pair transfer;
did not repeat612 branch-cut obstruction or646 free time-window result.
Kikukawa5.3 and Klich paired-CAR trace source checked. Nonzero physical mass
insertion explicitly chosen here, not attributed to a completed overlap Yukawa
prescription. Original Y, all32 modes, subgroup and fixed dictionary retained.
Pfaffian identity proved through the exact block recurrence U and CAR lift;
polynomial identity fixes phase globally, no principal determinant square root.
Positive per-step transfer does NOT make arbitrary conditional path weights
real or positive. Three-time original history actually has a complex trace.
BCH sign audited with explicit neutral density: H_a=H+a[Cdag,C]/2+O(a^2).
Vacuum constant retained together with the normal number term. Energy error
O(a^2) does not remove state error O(a). Full32 covariance and original five
source derivatives checked; explicit original four-mode density independently
checks Nambu covariance formula. Energy-only refitting fails finite-step state.
All common-limit assertions restricted to compact scalar backgrounds inside
F>0 and finite CAR; no exchange with full unbounded configuration integral.
Same auxiliary kernel ratio8/17 gives exact trace-norm convergence to unique
constant harmonic. Excited energies diverge1/a. No strongly continuous limit
on the full extended Hilbert space; the surviving subspace is distinguished.
Auxiliary holonomy derivative vanishes with same tail; original physical source
does not disappear. Fixed nonzero charged Higgs is not Gauss invariant; joint
covariance requires transforming/quantizing the scalar. No illegal claim of
commutation of fixed-background mass and full Gauss projection.
Original 2^(-32N) and lambda0^N explicitly restored in vacuum energy; removing
them from normalized matter states does not settle geometry/lapse or gravity.
Four tests passed,20 formulas checked. New result is the common interface and
the fixed-kernel limit conflict, not a new general theorem about all fermions.
Next655 prioritizes original vacuum normalization/geometry and spatial entry,
not auxiliary-device design or more isolated time-step optimizations.
"""
for name in ('research_note_654.md','joint_mass_state_time_limit.py',
             'joint_mass_state_time_limit_results.json','unified_physics_condition_ledger_654.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round654_drafts/final_review.txt',review)
print('654 frozen; verification/publication prepared')
