"""Prepare a verified working entry without changing scientific counters."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
mapping = {'770': '771', '769': '770', '3494': '3496', '3742': '3757',
           'finite_valence_entry': 'local_contact_entry'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda s: pattern.sub(lambda m: mapping[m[0]], s)
source = convert((HERE/'publish_round770_entry.py').read_text('utf8'))
start, stop = source.index('artifacts = ('), source.index('spec = ')
source = source[:start]+'''artifacts = (
    "round771_drafts/research_note_771_working.md",
    "round771_drafts/local_contact_entry.py",
    "round771_drafts/local_contact_entry_results.json",
    "round771_drafts/entry_scope_review.json",
    "round771_drafts/entry_condition_ledger.md",
)
'''+source[stop:]
source = source.replace('report_text.count("$$") == 4', 'report_text.count("$$") == 12')
source = source.replace('text_checks["display_formulas"] == 2', 'text_checks["display_formulas"] == 6')
source = source.replace('entry_calibrations_reproduced=1', 'entry_calibrations_reproduced=2')
start, stop = source.index('summary = '), source.index('planned = ')
summary = '**771局部接触工作已保存，正式仍770／3496：** [研究报告]({p}round771_drafts/research_note_771_working.md)明确密度配对变分的必要局部接触，剩余Ward改写为减除核的交换余项；成熟标量修复未直接覆盖原耦合对象。[条件增量]({p}round771_drafts/entry_condition_ledger.md)、[入口核验]({p}round771_drafts/local_contact_entry_checks.json)。同一实际态保持，绝对来源修复仍待证；不新增完成轮次，目标保持。'
source = source[:start]+'summary = '+repr(summary)+'\n'+source[stop:]
source = source.replace('**771固定背景来源入口已推进', '**771局部接触工作已保存')
# Verify every added report's local links, not just the lead report.
anchor = 'text_checks = core.text_checks(HERE/artifacts[0])'
source = source.replace(anchor, '''for extra in artifacts[4:]:
    for link in core.link_parser()((HERE/extra).read_text("utf8")):
        assert ((HERE/extra).parent/link).resolve().exists(), link
'''+anchor)
post = convert((HERE/'postcheck_round770_entry.py').read_text('utf8'))
post = post.replace('**771固定背景来源入口已推进', '**771局部接触工作已保存')
for name, content in [('publish_round771_entry.py', source), ('postcheck_round771_entry.py', post)]:
    with (HERE/name).open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)
review = dict(
    working_round=771, latest_completed_round=770, new_scientific_tests=0,
    primary_agent_review=True, independent_agent=False,
    analytic_checks=[
        'Density factor M is zeroth order and nondegenerate, not required positive.',
        'Product rule gives delta D=M^-1 delta L - (M^-1 delta M)D.',
        'Cyclic integration is used only with smooth state-minus-subtraction kernels and compact tests.',
        'Both bisolution identities give the sign in the remaining commutator trace.',
        'The full auxiliary operator is not asserted to be normally hyperbolic.',
        'Only the required combined kernels have coincidence limits; background jet order is not bounded by three.',
        'M contact and residual repair must satisfy the same permitted normalization category before acceptance.',
    ],
    numerical_checks=[
        'Full original 111 principal matrix and scalar fiber pairing derivative calibrate algebra only.',
        'W=0 in the finite calibration is not called a physical state; no actual two-point function is replaced.',
        'Two entry calibrations reproduce saved results; no new scientific test count.',
    ],
    literature=[
        dict(url='https://arxiv.org/html/gr-qc/0109048',
             location='Lemma 2.1 and Theorem 2.1',
             verified='Scalar eta_D=D/[2(D+2)], including external-potential exchange source.',
             limitation='No coupled graviton/YM/H5/ghost theorem supplied.'),
        dict(url='https://arxiv.org/html/1202.5107',
             location='Section 2, Theorem 6 and Section 6',
             verified='Scalar framework; zeta-to-actual-state correspondence assumes static setting and positive Euclidean operator.',
             limitation='No arbitrary dynamic coupled background/state correspondence supplied.'),
        dict(url='https://arxiv.org/html/hep-th/0306138',
             location='Sections 2.1, 2.3 and 4.1',
             verified='General local Laplace-type coefficients are usable after the object/sign map.',
             limitation='Their diagonal traces alone do not prove the required same-state local Ward repair.'),
    ],
    direct_reuse=[734, 735, 765, 768, 769, 770],
    spatial_results_not_reproved=[382, 383, 384, 386, 425, 522, 523],
    all_sector_Ward_repair_proven=False, counterterm_category_acceptance_proven=False,
    absolute_conserved_backreaction_proven=False, active_goal_unchanged=True,
    visual_inspection=False,
    next='Calculate original joint equation-remainder jets or prove an applicable variational local construction, then verify the complete permitted repair.')
with (HERE/'round771_drafts/entry_scope_review.json').open('x', encoding='utf8') as stream:
    json.dump(review, stream, ensure_ascii=False, indent=2)
print('Prepared 771 working entry publication.')
