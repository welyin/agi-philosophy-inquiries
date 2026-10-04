"""Prepare 772 working-entry publication, with unchanged scientific counters."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
mapping = {'771': '772', '770': '771', '3496': '3499', '3757': '3773',
           'local_contact_entry': 'first_response_entry'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda text: pattern.sub(lambda m: mapping[m[0]], text)
source = convert((HERE/'publish_round771_entry.py').read_text('utf8'))
source = source.replace('report_text.count("$$") == 12', 'report_text.count("$$") == 14')
source = source.replace('text_checks["display_formulas"] == 6', 'text_checks["display_formulas"] == 7')
source = source.replace('entry_calibrations_reproduced=2', 'bookkeeping_audit_reproduced=True, entry_calibrations_reproduced=0')
source = source.replace('background_split_anomaly_computed=False, all_sector_Ward_claimed=False',
                        'one_loop_single_source_Ward_inherited=True, interacting_BRST_lift_proven=False')
start, stop = source.index('summary = '), source.index('planned = ')
summary = '**772首阶响应入口已核，正式仍771／3499：** [工作报告]({p}round772_drafts/research_note_772_working.md)把771完整来源接入旧732／754约束和Cauchy响应；这是直接应用，不新增轮次。保原自由态，下一项核有限阶共同插入与物理态/记录提升，不以全阶背景输送拖延首阶均值。[条件增量]({p}round772_drafts/entry_condition_ledger.md)、[入口核验]({p}round772_drafts/first_response_entry_checks.json)。目标及旧限定结论保持。'
source = source[:start]+'summary = '+repr(summary)+'\n'+source[stop:]
source = source.replace('**772局部接触工作已保存', '**772首阶响应入口已核')
post = convert((HERE/'postcheck_round771_entry.py').read_text('utf8'))
post = post.replace('**772局部接触工作已保存', '**772首阶响应入口已核')
for name, content in [('publish_round772_entry.py', source), ('postcheck_round772_entry.py', post)]:
    with (HERE/name).open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)

review = dict(
    working_round=772, latest_completed_round=771, new_scientific_tests=0,
    primary_agent_review=True, independent_agent=False, visual_inspection=False,
    direct_reuse=[601, 623, 624, 625, 704, 732, 734, 735, 754, 756, 763, 765, 767, 768, 769, 770, 771],
    spatial_results_not_reproved=[382, 383, 384, 386, 425, 522, 523],
    analytic_review=[
        '771 conserved source and its original density pairing meet the old smooth Cauchy-source contract.',
        '754 supplies initial constraint data including color mean; homogeneous initial correction remains a preparation choice.',
        'No compact-time support of the actual source or zero initial correction is silently imposed.',
        'Fixed-source first-order mean identity is a direct application, not a new theorem or a finite-strength bound.',
        'Off-shell Hessian-generator identity uses density coordinates and does not mix adjoints from changing pairings.',
        'An O(epsilon) off-shell defect is not promoted to an obstruction to the already valid first-order mean on B0.',
        'Operator BRST-state recursion is explicitly conditional, with no interacting charge or domain claimed.',
        'Connected one-loop vertex count is finite; no-vertex two-point function and explicit counterterms are separate.',
        'Formal in-in retarded response is not required to be a symmetric single-history Hessian.',
    ],
    numerical_scope='No new numerical physical experiment. Exact integer graph bookkeeping and provenance audit only.',
    literature=[
        dict(url='https://arxiv.org/html/hep-th/9807078', location='Section 4, equation (4.8), Theorem 4; section 5',
             verified='Formal BRST physical-state lifting and positivity require the specified representation, hermitian nilpotent deformation, and free kernel positivity with null vectors in the image.',
             limitation='Explicit model implementation is QED; original mixed coupled model hypotheses are not established here.'),
        dict(url='https://arxiv.org/html/1306.1058', location='Section 5',
             verified='The quantum-gravity state construction is formal and conditional on free representation and interacting Ward data.',
             limitation='No convergent finite-strength interacting state or full original matter model. Positivity symbols checked against original 9807078 rather than copied from this HTML.'),
        dict(url='https://arxiv.org/html/gr-qc/0209075', location='Section III, equations (19)-(24)',
             verified='Causal in-in response involves a retarded stress commutator and local contact terms.',
             limitation='Matter semiclassical/large-N context does not supply all coupled quantum-gravity sectors.'),
        dict(url='https://arxiv.org/html/1804.07640', location='Theorem 3.13 and section 4',
             verified='Background comparison on physical cohomology requires an additional anomaly insertion condition; gravitational possibilities are left open.',
             limitation='This is inherited 769 scope, not a newly established failure.'),
    ],
    one_source_Ward_inherited=True, absolute_formal_response_is_direct_application=True,
    interacting_physical_state_proven=False, original_record_map_proven=False,
    finite_strength_self_consistency_proven=False, active_goal_unchanged=True,
    next='Keep original fixed background and examine the minimal finite-order common causal insertions, BRST-state/observable and record lift; do not reopen the conserved-source Cauchy problem.')
with (HERE/'round772_drafts/entry_scope_review.json').open('x', encoding='utf8') as stream:
    json.dump(review, stream, ensure_ascii=False, indent=2)
print('Prepared 772 entry publication; no completed-round counter change.')
