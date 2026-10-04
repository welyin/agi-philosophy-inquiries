"""Prepare 618 primary review and append-only publication."""
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

mapping={"joint_region_energy_gluing":"joint_geometric_boundary_matching",
         "616":"617","617":"618","618":"619",
         "3092":"3095","3095":"3098","1160":"1163","1163":"1166",
         "1907":"1915","1915":"1922"}
verify=replace((HERE/"verify_round617.py").read_text("utf8"),mapping)
verify=verify.replace("Verify original compact-group matching, full energy, and sources.",
                      "Verify original nonminimal boundary momenta and joint sources.")
verify=verify.replace(",'round618_drafts/region_source_entry.md'","")
verify=verify.replace("text['display_formulas']==12","text['display_formulas']==14")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round618.py",verify)
publish=replace((HERE/"publish_round617.py").read_text("utf8"),
                mapping|{"## 263.":"## 264.","## 168.":"## 169."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第618轮完成：** [原非最小耦合、几何边界与标量来源的共同匹配]({p}research_note_618.md)'
         '原F同时固定几何与五标量边界动量，混合Schur量恒M；'
         '无源类时拼接须同时匹配法向数据，框架变换须保边界标量来源。'
         '三组、十四式通过，最新618／3098，1166份编号科学文件、1922份保护证据。'
         '[核验]({p}research_round_618_checks.json)、[条件账]({p}unified_physics_condition_ledger_618.md)。'
         '限给定二导数经典作用；全物质边界、尺度及量子GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（618后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[619原物质边界电流与共同拼接]({p}round619_drafts/STATUS.md)，'
       '回查旋量与手征范围，再核原规范／费米边界及来源；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("规范匹配、非对角能源与共同来源","原F、共同边界动量与来源")
publish=publish.replace("表示切分不等于真实空间细化","无源接合的条件性匹配不选择维数")
publish=publish.replace("原规范区域拼接、完整演化与共同能源来源","原非最小耦合、几何边界与标量来源的共同匹配")
assert "hashlib.sha256" in publish
write("publish_round618.py",publish)
write("round618_drafts/research_note_618_draft.md",
      (HERE/"research_note_618.md").read_text("utf8"))
review="""618 primary-agent review. No independent agent review.
Prior 617 is verified progress. This round keeps condition integration first.
The original five real fields, F=M-phi^2/6, M=2, potential and target metric
are retained. The Lorentz Einstein/Jordan two-derivative action is an input.
Boundary is smooth timelike, outward normal spacelike with n^2=+1; curvature
K_ab=h_a^mu h_b^nu nabla_mu n_nu. Covariant induced metric is varied.
No null, corner, higher-curvature or quantum-gravity result is claimed.
The action includes +F K on that boundary. Scalar variation yields f K-v;
metric covector yields F(K h^ab-K^ab)+h^ab nF with the stated factor of 1/2.
This matches the primary appendix with its f_literature=F/2 convention.
Two sides share h and phi. Sums of outward momenta must cancel a specified
surface action. Holding both artificial boundaries fixed is not gluing.
Trace/scalar block has Schur F+3|grad F|^2/2=M; its inverse is explicit.
Traceless inverse divides F, so F=0 is excluded even though Schur remains M.
No-source matching removes first normal derivative jumps in adapted local
coordinates. This is an interface condition, not a bulk/global existence
theorem or a quantum constraint closure result. Degenerate and null cases
from mature literature are not excluded in general.
The same F gives the old nonflat Einstein target. Full boundary cotangent
transformation includes Pi_E h_J grad F in scalar momentum. Omitting it
is a change of boundary variation, not a harmless convention.
Test tension is expressly added only as a diagnostic. In Einstein variables
its action is -tau sqrt|h_E| F^(-3/2); scalar source must be varied from it.
Local junction data cancel both metric and scalar surface sources. No global
thin-shell spacetime satisfying all bulk constraints has been constructed.
The numeric action uses nonperiodic polynomial profiles and the original V.
Jordan and Einstein curvatures plus their own GHY terms are independently
integrated. Einstein lapse sqrt F is retained. Variation checks include
nonzero bulk Euler terms and nonzero endpoints; no false on-shell assumption.
The matrix Hessian has determinant -6M and is indefinite; invertibility is
not stability or positivity of a gravity quantum measure.
Gauge and fermion fields are zero for this scalar/gravity variation; arbitrary
Higgs profiles are off shell, not asserted a consistent full gauge truncation.
All matter boundary terms, full EFT, chirality, continuum maps and GR remain
open. Mature junction theory is not claimed new; only original-model mapping
and joint source/variational checks are the project increment.
"""
for n in ("research_note_618.md","joint_geometric_boundary_matching.py",
          "joint_geometric_boundary_matching_results.json","unified_physics_condition_ledger_618.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round618_drafts/final_review.txt",review)
print("618 review and publication prepared")

