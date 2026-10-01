"""Freeze 626 and prepare verified history-preserving publication."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    keys=[]
    for key in sorted(mapping,key=len,reverse=True):
        pattern=re.escape(key)
        if key.isdecimal():pattern=r'(?<!\d)'+pattern+r'(?!\d)'
        keys.append(pattern)
    return re.sub('|'.join(keys),lambda m:mapping[m.group()],text)


def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_source_preserving_compression':'joint_equilibrium_memory_matching',
         '624':'625','625':'626','626':'627',
         '3115':'3118','3118':'3121','1184':'1187','1187':'1190',
         '1966':'1973','1973':'1980'}
verify=replace((HERE/'verify_round625.py').read_text('utf8'),mapping)
verify=verify.replace('Verify joint finite CP compression of records and source coefficients.',
                      'Verify shared thermal memory and source-only locality obstruction.')
verify=verify.replace("text['display_formulas']==20","text['display_formulas']==14")
write('verify_round626.py',verify)
publish=replace((HERE/'publish_round625.py').read_text('utf8'),
                mapping|{'## 271.':'## 272.','## 176.':'## 177.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第626轮完成：** [共同参考态对局部有效作用的记忆限制]({p}research_note_626.md)'
         '原完整Gauss热态的均匀来源保跨时间关联；'
         '消去全部物质后只保源局部二次项，不能匹配两分离脉冲。'
         '同一初始能谱精确补全该来源部门。'
         '三组、十四式通过，最新626／3121，1190份编号科学文件、1980份保护证据。'
         '[核验]({p}research_round_626_checks.json)、[条件账]({p}unified_physics_condition_ledger_626.md)。'
         '数值限原32模式条件物质；不排除保留物质的局部EFT，连续与GR仍开放。')
order=('**当前执行顺序（626后，优先于下方历史安排）：** 继续整合共同条件，认知设计后置。'
       '接[627共同有效作用的保留与消去部门]({p}round627_drafts/STATUS.md)，'
       '核原规范、物质及几何来源的共同分工和态；不继续脉冲或记忆资源优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('同一有限过程的记录与来源匹配','共同参考态与局部作用的记忆限制')
publish=publish.replace('有限过程共同匹配不选择空间维数','时间来源记忆不选择空间维数')
publish=publish.replace('记录与二阶来源共用一个有限尺度过程','共同参考态对局部有效作用的记忆限制')
write('publish_round626.py',publish)
write('round626_drafts/research_note_626_draft.md',(HERE/'research_note_626.md').read_text('utf8'))
review='''626 primary-agent review. No independent agent review.
Previous goal turn made verified progress: 625 published and its state,
code, results and navigation receipts inspected before this continuation.
Current process inspection found no Python job started on the current date.
625 protected files and 626 frozen entry were not changed.
Historical search covered 281,305,558,585,590-591,600-603,612,620-625 and
the earlier CTP bridge. Conserved variance and static/retarded difference
are inherited 602-603 results; CTP and Gaussian spectral tools are mature.
New joint restriction concerns source-only temporal locality after eliminating
all retained memory, not the locality of a microscopic or Wilsonian action.
The actual full fixed-graph Gauss thermal state has strictly positive energy
variance by faithfulness, infinite-dimensional physical space and compact
resolvent. All moments exist. Uniform lapse is a declared original source.
There are no intermediate instruments in this diagnostic. A matching menu
containing this subcase must pass it; no claim made for all recorded histories.
H is not the full ADM constraint; positivity and Gibbs assumptions are not
silently moved to a total constrained quantum-gravity theory.
Exact CTP depends on the two elapsed-time integrals. Finite second moment
gives its actual Peano coefficient without source-limit exchange.
Any diagonal-supported local quadratic kernel has zero cross term for
disjoint smooth compact pulses. Plus/minus cases give the same local value
but exact variances 4Va2 and zero. Minimax coefficient bound follows directly.
This is a quadratic-coefficient bound, not a finite-amplitude uniform bound.
Boundary jets vanish; a boundary functional of the total source integral
already retains nonlocal initial data and is not in the excluded class.
Slow-pulse claim explicitly fixes bounded local coefficients and broadens
the time window. The fixed-window two-pulse bound does not need this limit.
Endpoint-fixed reparametrization keeps total elapsed time. Nonzero integral
branch differences are not called pure gauge. No local diffeomorphism or
gravitational Ward identity is claimed from this global source alone.
One common energy spectral measure restores the commuting source sector
exactly. It is not a new bath, an actual projective energy measurement, a
hidden-variable replacement of all quantum observables, or a memory device.
Gaussian constant noise is only a quadratic representation; third cumulant
and failure of independently redrawing energy are explicitly checked.
Original CAR, complex Yukawa/Majorana constants, scalar vacuum and beta
are unchanged. 32 positive BdG energies encode 32, not 64, physical modes.
Neutral Fock matrix characteristic function is independently checked.
Numerical variance is the frozen mass-sector value, not full Gauss variance.
Mixed conserved D is a prior result checked for continued use, not a new
fundamental theorem. No pulse optimization or autonomous hardware pursued.
Source literature was inspected: Hu-Verdaguer influence/noise distinction
and Crossley-Glorioso-Liu retained conserved modes. Their geometric or fluid
models and local KMS assumptions are not imported as proved project facts.
Three scientific groups passed on first run. A document patch had a malformed
add-file line and made no mutation; corrected before note freeze.
Next work returns to common retained/eliminated matter and geometry sources.
Unification goal is active, unchanged, and incomplete.
'''
for name in ('research_note_626.md','joint_equilibrium_memory_matching.py',
             'joint_equilibrium_memory_matching_results.json','unified_physics_condition_ledger_626.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round626_drafts/final_review.txt',review)
print('626 draft, review, verifier and publisher prepared')
