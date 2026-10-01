"""Freeze 625 and prepare history-preserving publication."""
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
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as stream:stream.write(text)


mapping={'joint_record_source_response':'joint_source_preserving_compression',
         '623':'624','624':'625','625':'626',
         '3112':'3115','3115':'3118','1181':'1184','1184':'1187',
         '1958':'1966','1966':'1973'}
verify=replace((HERE/'verify_round624.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original records, conditional source process and score constraints.',
                      'Verify joint finite CP compression of records and source coefficients.')
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==20")
verify=verify.replace(",'round625_drafts/development_diagnostic.md'","")
assert 'development_diagnostic' not in verify
anchor="    links=0\n"
assert verify.count(anchor)==1
verify=verify.replace(anchor,"    assert all(c>=32 or c in (10,9) for c in (HERE/'research_note_625.md').read_bytes())\n"+anchor)
verify='\n'.join(
    r"    assert not re.search(r'(?<!\\)\b(?:qquad|quad)\b', (HERE/'research_note_625.md').read_text('utf8'))"
    if line.strip().startswith('assert not re.search') else line
    for line in verify.split('\n'))
write('verify_round625.py',verify)
publish=replace((HERE/'publish_round624.py').read_text('utf8'),
                mapping|{'## 270.':'## 271.','## 175.':'## 176.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第625轮完成：** [记录与二阶来源共用一个有限尺度过程]({p}research_note_625.md)'
         '同一显式有限CP族共同逼近原有限历史的记录及零、一、二阶来源系数；'
         '初态、仪器尾项与Hessian必须共同输送。'
         '三组、二十式通过，最新625／3118，1187份编号科学文件、1973份保护证据。'
         '[核验]({p}research_round_625_checks.json)、[条件账]({p}unified_physics_condition_ledger_625.md)。'
         '限固定图与有限矩、有限来源菜单；数值仅为64维诊断，局域连续及GR仍开放。')
order=('**当前执行顺序（625后，优先于下方历史安排）：** 先尽量整合剩余共同条件，认知设计后置。'
       '接[626共同量子过程与连续作用的条件]({p}round626_drafts/STATUS.md)，'
       '按条件账核共同对象、态、来源与尺度，再处理相容性卡点；不继续谱精度优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('实际记录、条件选择与共同来源','同一有限过程的记录与来源匹配')
publish=publish.replace('来源分辨性不构成维数选择','有限过程共同匹配不选择空间维数')
publish=publish.replace('实际记录、条件响应与共同来源过程','记录与二阶来源共用一个有限尺度过程')
write('publish_round625.py',publish)
write('round625_drafts/research_note_625_draft.md',(HERE/'research_note_625.md').read_text('utf8'))
review='''625 primary-agent review. No independent agent review.
624 verified and published before this dependent unit. 592 already proves
CP finite-history value approximation; 623 gives the original operator
domain and real source expansion; 624 supplies actual recorded histories.
New connection is their common finite CP model, including its own compressed
state, instruments, and actual time evolution, through weak order two.
The original full graph, Gauss, CAR content and source family are unchanged.
P_R is fixed at the base parameter; no moving-projector derivative omitted.
Ground choice need not be unique or pointwise positive after adding CAR.
The tail channel preserves arbitrary passive reference marginals and has
finite Kraus representation on the retained finite-dimensional system.
The instrument tail is positive and essential to trace preservation.
Naive system-only energy control fails for the Stinespring adjoint.
An explicit fixed enlarged space and positive auxiliary energy labels give
exact forward and backward intertwining and graph convergence. These are
proof weights, not a physical free reset or an autonomous apparatus.
Source derivatives have uniform A-relative bounds on the common domain.
Two insertions are placed on opposite inner-product sides; no unjustified
operator product on D(A), strong second derivative, or D(H^2) is asserted.
Telescoping is restricted to bounded words or compatible graph-space maps.
Finite A2 initial moment permits dominated trace sums. Finite histories and
fixed finite source menus are explicit; no uniform cutoff rate or arbitrary
observable, all-times, continuum, or full source-space Frechet claim.
Original effects >= I/4 preserve conditional probability denominators.
Monitored source coefficients are not relabelled as unmonitored kernels.
Numerics use original radial potential, measure, inverse metric and readouts
in a declared 64-dimensional finite-difference quadratic-form diagnostic.
This is neither the original full graph spectrum nor a certified continuum
discretization. Intermediate ranks are visibly inaccurate; full-rank equality
is an endpoint algebra check, not an infinite-dimensional convergence proof.
Finite-difference step-halving error is reported. Linear real log artifact
is not claimed as physical attenuation. Wrong Gibbs, absent Hessian and
missing instrument-tail substitutions are checked with unchanged parameters.
All three saved groups passed on first run; verifier reruns actual outputs.
Nine carriage-return TeX escapes and three bare qquad tokens were repaired
in the unpublished note, before draft freeze. No historical files changed.
Condition ledger adds an execution grouping, not another research round or
claim that related conditions imply each other. Design is deferred, while
resource requirements remain. Next entry returns to common scale, action
and source compatibility, without additional rank optimization.
Unified physics goal remains active and unchanged.
'''
for name in ('research_note_625.md','joint_source_preserving_compression.py',
             'joint_source_preserving_compression_results.json','unified_physics_condition_ledger_625.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round625_drafts/final_review.txt',review)
print('625 draft, review, verifier and publisher prepared')
