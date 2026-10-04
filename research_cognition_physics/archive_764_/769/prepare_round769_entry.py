"""Prepare publishers for the 769 working report, without science counters."""
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
DRAFT = HERE / "round769_drafts"
mapping = {"768": "769", "767": "768", "3488": "3491", "3714": "3726",
           "joint_source_entry": "absolute_source_entry"}
pattern = re.compile("|".join(map(re.escape, sorted(mapping, key=len, reverse=True))))


def replace(source):
    return pattern.sub(lambda match: mapping[match[0]], source)


source = replace((HERE / "publish_round768_entry.py").read_text("utf8"))
source = source.replace("import ast\n", "import ast\nimport importlib.util\n")
source = source.replace('artifacts = ("round769_drafts/absolute_source_entry.md",)', '''artifacts = (
    "round769_drafts/research_note_769_working.md",
    "round769_drafts/bv_source_entry_calibration.py",
    "round769_drafts/bv_source_entry_calibration_results.json",
    "round769_drafts/entry_condition_ledger.md",
    "round769_drafts/entry_scope_review.json",
)
spec = importlib.util.spec_from_file_location("entry769calibration", HERE/artifacts[1])
calibration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calibration)
assert calibration.run() == core.read(HERE/artifacts[2])
ast.parse((HERE/artifacts[1]).read_text("utf8"))''')
source = source.replace('report_text.count("$$") == 6', 'report_text.count("$$") == 14')
source = source.replace('text_checks["display_formulas"] == 3', 'text_checks["display_formulas"] == 7')
a, b = source.index("summary = "), source.index("planned = ")
source = source[:a] + '''summary = "**769绝对来源工作报告已保存，正式仍768／3491：** [研究报告]({p}round769_drafts/research_note_769_working.md)核成熟EFT/BV定理的原对象映射，区分规范反常、来源接触及背景比较；两组原势/规范固定局部校准已复算，未宣称实际量子反常或全部门Ward完成。下一项先接固定背景首阶绝对来源，不先要求全阶背景独立。[条件增量]({p}round769_drafts/entry_condition_ledger.md)、[入口核验]({p}round769_drafts/absolute_source_entry_checks.json)。不新增完成轮次，目标保持。"
''' + source[b:]
source = source.replace("**769完整来源入口已推进", "**769绝对来源工作报告已保存")
source = source.replace("new_scientific_tests=0, entry_text_checks=text_checks,",
                        "new_scientific_tests=0, entry_calibrations_reproduced=2, entry_text_checks=text_checks,")
source = source.replace("inherited_768_free_Hadamard_state_preserved=True,",
                        "inherited_767_768_free_state_preserved=True, background_split_anomaly_computed=False,")
with (HERE / "publish_round769_entry.py").open("x", encoding="utf8", newline="\n") as stream:
    stream.write(source)

post = replace((HERE / "postcheck_round768_entry.py").read_text("utf8"))
post = post.replace("**769完整来源入口已推进", "**769绝对来源工作报告已保存")
with (HERE / "postcheck_round769_entry.py").open("x", encoding="utf8", newline="\n") as stream:
    stream.write(post)

review = dict(
    working_round=769, completed_round=768, numbered_scientific_tests_added=0,
    primary_agent_review_completed=True, independent_agent_review=False,
    visual_checks=False, old_scope_preserved=True, active_goal_unchanged=True,
    results_scope="Two off-shell finite-jet identity checks, not computed anomalies or loop sources.",
    conditional_derivation="Assuming a common renormalized BV vertex functional, its one-loop source obeys R*J+R1*E=0. Off-shell response includes the generator contact. In-in source connection is still required.",
    original_dedup=[580, 600, 628, 629, 735, 765, 766, 767, 768],
    sources=[
        dict(url="https://arxiv.org/html/1501.07014", sections="2.3 and 8",
             use="Formal nonrenormalizable EFT anomaly theorem, with explicit completeness and one-loop hypotheses.",
             limitation="No blanket signature of the original coupled real-time state/source or its background transport."),
        dict(url="https://arxiv.org/html/1306.1058v5", sections="3.4, equations 64-67",
             use="Local BV master equation and pure-gravity cohomology argument.",
             limitation="Pure-gravity anomaly reasoning does not cover original U(1) and chiral matter directly."),
        dict(url="https://arxiv.org/html/1804.07640", sections="3.3.4, 4, remark 4.1",
             use="Background-dependent gauge-fixing insertion is an additional check for quantum background comparison.",
             limitation="Potential obstructions are not proven occurrences. Its full hypotheses remain to be checked for original coupled model."),
        dict(url="https://link.springer.com/article/10.1007/s00023-023-01368-0", sections="Theorem 1",
             use="Locally covariant Adler-Bardeen bridge for renormalizable YM.",
             limitation="Not a direct theorem for dynamical gravity, H5 sigma model and the original nonzero matter background."),
    ],
    next="Construct a same-state all-sector local source at fixed background and formal first order; do not delay this by requiring all-order background independence.",
)
with (DRAFT / "entry_scope_review.json").open("x", encoding="utf8") as stream:
    json.dump(review, stream, ensure_ascii=False, indent=2)
print("Prepared 769 entry publishers and scope review.")
