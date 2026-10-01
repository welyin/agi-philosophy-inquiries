"""Prepare 619 primary review and append-only publication."""
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

mapping={"joint_geometric_boundary_matching":"joint_chiral_boundary_gluing",
         "617":"618","618":"619","619":"620",
         "3095":"3098","3098":"3101","1163":"1166","1166":"1169",
         "1915":"1922","1922":"1930"}
verify=replace((HERE/"verify_round618.py").read_text("utf8"),mapping)
verify=verify.replace("Verify original nonminimal boundary momenta and joint sources.",
                      "Verify original Nambu boundary ranks, transmission, and frame sources.")
verify=verify.replace("'round620_drafts/STATUS.md')",
                      "'round620_drafts/STATUS.md','round619_drafts/boundary_current_entry.md')")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round619.py",verify)
publish=replace((HERE/"publish_round618.py").read_text("utf8"),
                mapping|{"## 264.":"## 265.","## 169.":"## 170."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第619轮完成：** [原手征物种、边界电流与共同区域拼接]({p}research_note_619.md)'
         '原Nambu物种的指定对称局部反射需32维酉匹配而交织秩至多2；'
         '整体透射可保原质量与来源，同一F还固定旋量通量和质量权重。'
         '三组、十四式通过，最新619／3101，1169份编号科学文件、1930份保护证据。'
         '[核验]({p}research_round_619_checks.json)、[条件账]({p}unified_physics_condition_ledger_619.md)。'
         '限明示连续分支与边界规则；重建、尺度及量子GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（619后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[620区域对象、共同量子态与几何来源]({p}round620_drafts/STATUS.md)，'
       '回查旧态、时间谱和来源结果，补真实联合接口；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("原F、共同边界动量与来源","手征区域、法向流与共同来源")
publish=publish.replace("无源接合的条件性匹配不选择维数","连续旋量边界要求不推出维数")
publish=publish.replace("原非最小耦合、几何边界与标量来源的共同匹配","原手征物种、边界电流与共同区域拼接")
assert "hashlib.sha256" in publish
write("publish_round619.py",publish)
write("round619_drafts/research_note_619_draft.md",
      (HERE/"research_note_619.md").read_text("utf8"))
review="""619 primary-agent review. No independent agent review.
Previous 618 is verified progress. This is common-condition integration.
Lorentz canonical Weyl principal symbol, four dimensions and the spin frame
are explicit continuous-branch inputs, not the completed 598 graph limit.
Original Q,u,d,L,e,nu modules, charges, handedness and complex Y are retained.
Nambu basis is (c_p,c^dagger_-p); hole principal normal block is +A^T, not
-A^T, because its kinetic matrix is -h_pr(-p)^T. Gauge and spin hole
generators are -q and -J^T. Code includes both signs and the original pairing.
Maximal null boundary subspaces of inertia (32,32) are graphs of unitaries.
This is a necessary local condition; vanishing a single diagonal current
or imposing all zero traces does not prove a self-adjoint domain.
At homogeneous zero Higgs and constant singlet, original gauge symmetry and
axial spin rotation are unbroken. Local linear derivative-free boundaries
preserving both require intertwining equal (q,j) weights. The actual table
has only two allowed entries, both neutral nu. Rank <=2 cannot reach32.
Allowing all complex Nambu unitaries is an enlargement before particle-hole
compatibility, so the obstruction is not based on excluding Majorana or
antilinear particle boundary conditions. Nonlinear boundaries are not covered.
Nonzero Higgs-only, derivative, nonlocal, symmetry breaking and additional
boundary-state branches are explicitly outside the no-go. No general claim
about Standard Model boundaries or all interaction mechanisms is made.
Transmission uses combined normal form diag(A,-A); its graph has dimension
64 in128. The numerical normalized graph checks linear algebra, not cloning.
Actual L2 restriction on disjoint regions is unitary. For flat full-space
canonical principal symbol and bounded Hermitian lower-order terms, domain
is H1 and self-adjoint. Matching traces characterize global H1; a constant
internal gauge change gives the conjugate cut domain. Variable transitions
require the corresponding connection and are not silently omitted.
Original masses and their parameter derivatives transform with fixed T;
no field-dependent derivative of T is discarded. No normal thermal trace or
interacting chiral measure is inferred from one-particle self-adjointness.
The charged Gaussian is an exact characteristic in the stipulated continuum
background, with Higgs zero but original neutral Majorana nonzero. It is not
a prepared Gauss-invariant full state or an Einstein background solution.
Spinor F^(-3/4) and mass F^(-1/2) accompany g_E=F g_J. Surface pairing and
spacetime mass density retain all factors. Source differentiation compares
the same pulled-back family, not incompatible fixed frames.
Nambu overall half factors are common. Pairing conservation is not ordinary
Majorana particle-number conservation. Coefficient bilinear tests do not
establish a quantum measure, anomaly cancellation, or source renormalization.
617 finite-graph,618 classical geometry,619 continuum fibre remain distinct
until a common state/scale/physical-source map is supplied. Full goal stays open.
"""
for n in ("research_note_619.md","joint_chiral_boundary_gluing.py",
          "joint_chiral_boundary_gluing_results.json","unified_physics_condition_ledger_619.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round619_drafts/final_review.txt",review)
print("619 review and publication prepared")

