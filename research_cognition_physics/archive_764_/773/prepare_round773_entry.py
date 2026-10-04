"""Prepare a scope-only working entry; no new scientific result or counters."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
mapping = {'772': '773', '771': '772', '3499': '3502', '3773': '3789',
           'first_response_entry': 'local_insertion_entry'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda text: pattern.sub(lambda m: mapping[m[0]], text)
source = convert((HERE/'publish_round772_entry.py').read_text('utf8'))
start, stop = source.index('artifacts = ('), source.index('report_text = ')
source = source[:start] + '''artifacts = (
    "round773_drafts/research_note_773_working.md",
    "round773_drafts/entry_condition_ledger.md",
    "round773_drafts/entry_scope_review.json",
)
''' + source[stop:]
source = source.replace('report_text.count("$$") == 14', 'report_text.count("$$") == 6')
source = source.replace('text_checks["display_formulas"] == 7', 'text_checks["display_formulas"] == 3')
source = source.replace('for extra in artifacts[4:]:', 'for extra in artifacts[1:2]:')
source = source.replace('bookkeeping_audit_reproduced=True', 'primary_source_scope_review_completed=True')
summary = '**773共同插入入口已核，正式仍772／3502：** [工作报告]({p}round773_drafts/research_note_773_working.md)区分自由态收缩与局部反项，明确当前实际作用及观测插入的障碍验收；不额外要求全H¹消失，也不把771一次来源直接扩大为共同相互作用处方。[条件增量]({p}round773_drafts/entry_condition_ledger.md)、[入口核验]({p}round773_drafts/local_insertion_entry_checks.json)。本项不增加完成轮次，目标保持。'
start, stop = source.index('summary = '), source.index('planned = ')
source = source[:start]+'summary = '+repr(summary)+'\n'+source[stop:]
source = source.replace('**773首阶响应入口已核', '**773共同插入入口已核')
post = convert((HERE/'postcheck_round772_entry.py').read_text('utf8'))
post = post.replace('**773首阶响应入口已核', '**773共同插入入口已核')

review = dict(
    working_round=773, latest_completed_round=772, cumulative_scientific_tests=3502,
    new_scientific_tests=0, primary_agent_review=True, independent_agent=False,
    visual_inspection=False, scope='Primary-source applicability and next-insertion audit only.',
    preserved_results=[732, 735, 754, 767, 768, 769, 770, 771, 772],
    direct_reuse=['772 free representation and mean/source/initial-data transport',
                 '771 single-source local Ward normalization',
                 '770/772 finite ordinary external-leg vertex counting',
                 '769/772 distinction between background transport and fixed-background first mean'],
    primary_source_checks=[
        dict(url='https://arxiv.org/pdf/1803.10235',
             locations=['Theorems 3, 10, 11, 12', 'Equations (219)-(223) and adjoining qualifications'],
             scope='Action anomaly removal precedes interacting differential; observable/contact lifting has additional local-functional hypotheses. No automatic promotion of a source-only scheme.'),
        dict(url='https://arxiv.org/html/1705.03480', locations=['Sections 2.1-2.4'],
             scope='Background-field counterterm structure assumes the specified gauge algebra, absence of measure anomaly and local remaining divergences. No original real-time state constructed.'),
    ],
    reviewed_distinctions=[
        'Algebraic Fock contraction is not a support-preserving local-jet contraction.',
        'Actual anomaly coefficients vanishing is not full ghost-number-one cohomology vanishing.',
        'Task-specific first-order observable lifting needs the actual obstruction to be exact in the allowed class.',
        'Local counterterm consistency includes field derivatives/subinsertions at fixed free background.',
        'Ordinary external-leg counting is not the full BRST/antifield insertion menu.',
        'An isolated nonintegrable contact does not by itself disprove a full retarded source prescription.',
    ],
    interacting_action_anomaly_computed=False, local_primitive_constructed=False,
    interacting_state_or_instrument_proven=False, original_Q_continuum_map_proven=False,
    original_single_source_Ward_preserved=True,
    free_Fock_contraction_used_as_local_counterterm=False,
    stronger_global_cohomology_axiom_added=False, active_goal_unchanged=True)

targets = [('publish_round773_entry.py', source), ('postcheck_round773_entry.py', post),
           ('round773_drafts/entry_scope_review.json', json.dumps(review, ensure_ascii=False, indent=2))]
assert all(not (HERE/name).exists() for name, _ in targets)
for name, content in targets:
    with (HERE/name).open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)
print('Prepared 773 scope-only entry; scientific counters unchanged.')
