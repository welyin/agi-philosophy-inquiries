"""Prepare working-entry publication, without incrementing scientific counts."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
mapping = {'769':'770', '768':'769', '3491':'3494', '3726':'3742',
           'absolute_source_entry':'finite_valence_entry'}
pattern = re.compile('|'.join(map(re.escape, sorted(mapping, key=len, reverse=True))))
convert = lambda s: pattern.sub(lambda m: mapping[m[0]], s)

source = convert((HERE/'publish_round769_entry.py').read_text('utf8'))
start, stop = source.index('artifacts = ('), source.index('spec = ')
source = source[:start]+'''artifacts = (
    "round770_drafts/research_note_770_working.md",
    "round770_drafts/finite_valence_entry.py",
    "round770_drafts/finite_valence_entry_results.json",
    "round770_drafts/entry_scope_review.json",
)
'''+source[stop:]
source = source.replace('report_text.count("$$") == 14', 'report_text.count("$$") == 4')
source = source.replace('text_checks["display_formulas"] == 7', 'text_checks["display_formulas"] == 2')
source = source.replace('entry_calibrations_reproduced=2', 'entry_calibrations_reproduced=1')
source = source.replace('inherited_767_769_free_state_preserved', 'inherited_767_768_free_state_preserved')
start, stop = source.index('summary = '), source.index('planned = ')
summary = '**770固定背景来源入口已推进，正式仍769／3494：** [工作报告]({p}round770_drafts/research_note_770_working.md)区分背景展开与涨落展开：完整背景可保在传播子中，首阶一外腿来源只需原三阶顶点；无需先证明无限背景级数收敛。局部规范化和Ward修复仍须实核，未签收绝对来源。[入口核验]({p}round770_drafts/finite_valence_entry_checks.json)。计数校准不新增科学轮次，目标保持。'
source = source[:start]+'summary = '+repr(summary)+'\n'+source[stop:]
source = source.replace('**770绝对来源工作报告已保存', '**770固定背景来源入口已推进')
post = convert((HERE/'postcheck_round769_entry.py').read_text('utf8'))
post = post.replace('**770绝对来源工作报告已保存', '**770固定背景来源入口已推进')
for name, content in [('publish_round770_entry.py', source), ('postcheck_round770_entry.py', post)]:
    with (HERE/name).open('x', encoding='utf8', newline='\n') as stream:
        stream.write(content)
review = dict(
    working_round=770, latest_completed_round=769, new_scientific_tests=0,
    reviewed_by_primary_agent=True, independent_agent=False,
    result='One-loop fixed-background tadpole uses exact third fluctuation derivatives; background amplitudes need not be Taylor truncated.',
    assumptions=['On-shell background; full quadratic operator resummed.',
                 'Centered inherited free state and formal first-order quantum expansion.',
                 'Remaining interacting bulk vertices have fluctuation valence at least three.'],
    not_conclusions=['Gauge anomaly cancellation or local counterterm primitive existence.',
                     'Finite degree in actual background fields or convergence of full interacting perturbation theory.',
                     'Physical realization of every enumerated valence pattern.'],
    primary_sources=['https://arxiv.org/html/1501.07014', 'https://arxiv.org/pdf/1803.10235'],
    scope='An entry simplification, not a new completed scientific round.',
    next='Actual local all-sector source subtraction and its on-shell Ward defect; keep exact background coefficients.',
    active_goal_unchanged=True)
with (HERE/'round770_drafts/entry_scope_review.json').open('x', encoding='utf8') as stream:
    json.dump(review, stream, ensure_ascii=False, indent=2)
print('Prepared 770 working entry publication.')
