"""Prepare 615 review and publication without changing frozen evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replace(text,mapping):
    parts=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r"(?<!\d)"+p+r"(?!\d)"
        parts.append(p)
    return re.sub("|".join(parts),lambda m:mapping[m.group()],text)

def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open("x",encoding="utf8",newline="\n") as f:f.write(text)

mapping={"joint_spinor_subgroup_mass":"joint_subgroup_measure_source",
         "613":"614","614":"615","615":"616",
         "3084":"3087","3087":"3090","1151":"1154","1154":"1157",
         "1884":"1891","1891":"1900"}
verify=replace((HERE/"verify_round614.py").read_text("utf8"),mapping)
verify=verify.replace("Verify the original subgroup and CAR mass dictionary.",
                      "Verify the exact finite subgroup measure and common source.")
verify=verify.replace("'round616_drafts/STATUS.md')",
    "'round616_drafts/STATUS.md','round615_drafts/measure_source_entry.md','round615_drafts/probe_auxiliary.py')")
verify=verify.replace("text['display_formulas']==12","text['display_formulas']==14")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round615.py",verify)
publish=replace((HERE/"publish_round614.py").read_text("utf8"),
                mapping|{"## 260.":"## 261.","## 165.":"## 166."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第615轮完成：** [原子群的精确辅助测度与共同响应]({p}research_note_615.md)'
         '同一原子群平坦族的辅助Pfaffian积分已精确完成且严格为正；'
         '它与物理行列式共同产生来源，删除辅助项会漏掉非零响应。'
         '三组、十四式通过，最新615／3090，1157份编号科学文件、1900份保护证据。'
         '[核验]({p}research_round_615_checks.json)、[条件账]({p}unified_physics_condition_ledger_615.md)。'
         '限明确有限欧氏族；一般局域性、非零质量、重建及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（615后，优先于下方历史安排）：** 继续先合并共同条件，认知设计后置。'
       '接[616原物质、全局几何与量子态]({p}round616_drafts/STATUS.md)，'
       '核spin存在／选择、原群及状态来源，不继续孤立优化本轮积分；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("原群与原质量的共同表示字典","同一测度、物理权重与共同来源")
publish=publish.replace("表示容器不等于物理规范群扩大","有限测度正性与全局物理前提")
publish=publish.replace("原物质、全局商群与成熟手征表示的共同字典","原子群的精确辅助测度与共同响应")
assert "hashlib.sha256" in publish
write("publish_round615.py",publish)
write("round615_drafts/research_note_615_draft.md",
      (HERE/"research_note_615.md").read_text("utf8"))
review="""615 primary-agent review. No new independent agent review.
The prior goal response only restated the research ordering: no-progress.
This unit changes the evidence with an exact integral and source calculation.
Use Kikukawa's specified measure, not an asserted completed general chiral
theory. D is half the previous convention, so GW has coefficient 2.
One site in each of four Euclidean directions, three identity spatial links,
original hypercharge temporal holonomy, antiperiodic temporal fermions.
Four dimensions, the boundary choice and m0=1 are explicit inputs.
Scalar background phi=0 means the old mass matrix vanishes; the original
Yukawa parameters are not re-fitted, and scalar fluctuations are not integrated.
No claim of solving the original massive dynamical bosonic model.
The Wilson gap is exactly 1. The chirality bases are actual rotating spinors.
Their full Jacobian is independent of theta; a fixed phase sets normalization.
Pairing uses i gamma5 C_D and the 614 Clifford matrices, not naive C_D.
Pfaffian reduction is analytic; pivoted skew elimination checks full matrices.
The uniform S9 integral has Beta(3,2) radial measure and exact rational
coefficients. Six Gauss nodes suffice for its degree-11 integrand.
a and b cannot both vanish. Positivity is proved for the whole chosen family,
not inferred from samples. This does not give all-volume locality.
At pi/3 the pure weak auxiliary configuration is singular; the integrated
measure remains positive. An individual zero is not failure of the integral.
The physical determinant has zeros outside the declared small-angle window;
global positivity of the auxiliary factor is not global nonzero full Z.
theta is observable holonomy variation, not a pure gauge Ward direction.
The source includes both physical and auxiliary terms. The induced action
has the opposite sign to log Z, consistent with the source convention.
Alternative fixed-r measure preserves the subgroup and second moments but
is not proven a full valid chiral theory. S9 remains fixed in the main candidate.
Original subgroup has no vector singlet, so the Spin9 sufficient proof cannot
be applied by inclusion. Its failure is not a no-go; this family is positive.
Conditional pullback does not identify a geometric lambda or calculate GR.
Fixed-link Euclidean Z is not a normal thermal trace or time reconstruction.
Final scope keeps spin/geometry, full mass, locality, states, scales and GR open.
The existing entry and exploratory script are preserved as support, not new tests.
"""
for n in ("research_note_615.md","joint_subgroup_measure_source.py",
          "joint_subgroup_measure_source_results.json","unified_physics_condition_ledger_615.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round615_drafts/final_review.txt",review)
print("615 review and publication prepared")
