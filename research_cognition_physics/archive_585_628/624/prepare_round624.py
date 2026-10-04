"""Freeze 624 and prepare verified history-preserving publication."""
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


mapping={'joint_operator_domain_completion':'joint_record_source_response',
         '622':'623','623':'624','624':'625',
         '3109':'3112','3112':'3115','1178':'1181','1181':'1184',
         '1951':'1958','1958':'1966'}
verify=replace((HERE/'verify_round623.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original operator domain and actual thermal process expansion.',
                     'Verify original records, conditional source process and score constraints.')
verify=verify.replace("text['display_formulas']==18","text['display_formulas']==14")
verify=verify.replace("'round625_drafts/STATUS.md')\n    preserved=",
                      "'round625_drafts/STATUS.md','round624_drafts/development_diagnostic.md')\n    preserved=")
assert "'round624_drafts/development_diagnostic.md'" in verify
write('verify_round624.py',verify)
publish=replace((HERE/'publish_round623.py').read_text('utf8'),
                mapping|{'## 269.':'## 270.','## 174.':'## 175.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第624轮完成：** [实际记录、条件响应与共同来源过程]({p}research_note_624.md)'
         '原仪器与完整CAR来源共享有限历史二阶过程；条件筛选须保记录概率导数，'
         '实际来源得分限制记录合并后的分辨性。'
         '三组、十四式通过，最新624／3115，1184份编号科学文件、1966份保护证据。'
         '[核验]({p}research_round_624_checks.json)、[条件账]({p}unified_physics_condition_ledger_624.md)。'
         '限固定图与有限矩历史；连续、动态引力与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（624后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[625同一有限尺度映射的记录与来源合同]({p}round625_drafts/STATUS.md)，'
       '回查原谱近似和虚跃迁，核记录、噪声、迟致及接触项的共同匹配；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('共同算符域、原物质与真实热过程','实际记录、条件选择与共同来源')
publish=publish.replace('原完整过程展开不构成维数选择','来源分辨性不构成维数选择')
publish=publish.replace('原完整量子过程的共同算符域与真实热响应','实际记录、条件响应与共同来源过程')
assert 'hashlib.sha256' in publish
write('publish_round624.py',publish)
write('round624_drafts/research_note_624_draft.md',(HERE/'research_note_624.md').read_text('utf8'))
review='''624 primary-agent review. No independent agent review.
623 verified and published before beginning this dependent unit. 624 entry
and 586/592/593/598/617 records audited; 592 already contains graph-domain
and H2 budgets, and 586 already contains Fisher conditional decomposition.
Neither is reclassified as a new theorem. Added connection is the actual
full CAR common-source history with original measurement insertions.
CAR potential commutes with primitive L(s). Original lower potential bound
controls the unchanged kinetic commutator by shifted A. Both L and adjoint
preserve common D, enabling weak two-insertion argument for finite words.
Finite A2 suffices; Gibbs commutation is not required in this extension.
Original strictly positive effects give p_history >= 4^-k, not a uniform
infinite-history bound. No extra instrument or autonomous hardware chosen.
Same-source branch Z diagonal is record probability, not one. Summed branch
Z for unequal sources is generally the monitored process, not unmonitored
623. Conditional normalization keeps parameter-dependent probabilities.
Bounded observable alone is not enough for second weak expansion; D and
adjoint-D preservation or equivalent control explicitly required.
Future effect and past state distinction checked against primary Gammelmark
et al paper. Source insertion is a Duhamel functional derivative, not an
unproved independently executable exp(-i source) or backward physical time.
HS Gram positivity and normalized curvature refer to the fixed Kraus
representation/purification, not accessible mixed-state QFI or spacetime.
Numerics evaluate true original full-H instantaneous f(s) commutator, for
which all matrix potentials and spectators cancel analytically. They do not
simulate a replacement or full graph time evolution. Singlet current is
not a spacetime propagation velocity. Selected read occurs after flow.
Initial nonzero-magnitude diagnostic failed for the smaller branch at
8.84e-5 versus an arbitrary 1e-4 threshold. Diagnostic preserved separately;
no parameter changed. Final check compares correction to floating spacing.
Three final groups reproduce. Coarse score identity is inherited algebra
applied to actual source probabilities, not a newly invented theorem or
replacement of coarse reports by a direct square-root instrument.
Next unit must integrate original scale and source contracts, no endless
precision tuning or new cognitive hardware. Full unification remains open.
'''
for name in ('research_note_624.md','joint_record_source_response.py',
             'joint_record_source_response_results.json','unified_physics_condition_ledger_624.md',
             'round624_drafts/development_diagnostic.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round624_drafts/final_review.txt',review)
print('624 draft, review, verifier and publisher prepared')
