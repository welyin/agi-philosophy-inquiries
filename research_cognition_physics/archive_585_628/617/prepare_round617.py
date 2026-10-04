"""Prepare 617 primary review and append-only publication."""
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

mapping={"joint_spin_thermal_source":"joint_region_energy_gluing",
         "615":"616","616":"617","617":"618",
         "3090":"3092","3092":"3095","1157":"1160","1160":"1163",
         "1900":"1907","1907":"1915"}
verify=replace((HERE/"verify_round616.py").read_text("utf8"),mapping)
verify=verify.replace("Verify original neutral spin obstruction and thermal source contract.",
                      "Verify original compact-group matching, full energy, and sources.")
verify=verify.replace("'round618_drafts/STATUS.md')",
                      "'round618_drafts/STATUS.md','round617_drafts/region_source_entry.md')")
verify=verify.replace("==(2,0,0)","==(3,0,0)").replace("fresh_tests=dict(run=2,","fresh_tests=dict(run=3,")
verify=verify.replace("text['display_formulas']==11","text['display_formulas']==12")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round617.py",verify)
publish=replace((HERE/"publish_round616.py").read_text("utf8"),
                mapping|{"## 262.":"## 263.","## 167.":"## 168."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第617轮完成：** [原规范区域拼接、完整演化与共同能源来源]({p}research_note_617.md)'
         '原紧商群的匹配切分保持原Gauss过程、CAR、热态和几何来源；'
         '非对角电动能须共同分配，重复计能或删剪切均有明确差额。'
         '三组、十二式通过，最新617／3095，1163份编号科学文件、1915份保护证据。'
         '[核验]({p}research_round_617_checks.json)、[条件账]({p}unified_physics_condition_ledger_617.md)。'
         '仅表示切分，非独立区域或真实细化；连续与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（617后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[618区域拼接、共同作用与几何边界]({p}round618_drafts/STATUS.md)，'
       '先回查旧边界和引力变分结果，再核同一来源的缺口；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("原物种、spin存在与热态来源","规范匹配、非对角能源与共同来源")
publish=publish.replace("条件性spin存在不选择维数或空间边界","表示切分不等于真实空间细化")
publish=publish.replace("原中性费米子、旋量结构与热态边界的共同条件","原规范区域拼接、完整演化与共同能源来源")
assert "hashlib.sha256" in publish
write("publish_round617.py",publish)
write("round617_drafts/research_note_617_draft.md",
      (HERE/"research_note_617.md").read_text("utf8"))
review="""617 primary-agent review. No new independent agent review.
This is condition integration, not a cognitive device design.
The compact group is the original Z6 quotient, with normalized Haar measure.
Multiplication and cutoff gauge are well-defined on that quotient. In (U,A)
coordinates, cut invariance is independence of A, establishing onto, not only
an isometric embedding. Peter-Weyl uses actual admissible representations.
Whole-graph map is identity on scalars and original CAR modes. Endpoint Gauss
intertwines; cut Gauss gives matching, not independent regional preparability.
Momentum on B must be transported by Ad(A). Coefficients depend on A, so B
derivatives do not differentiate them. Both lifted derivatives preserve the
matching subspace. Non-Abelian and U1 directions follow the same rule.
The actual 589 coefficient is a positive matrix across outgoing directions,
not a collection of diagonal Casimirs. Positive split KA+KB=K is sufficient;
it is a bookkeeping freedom on the matched space, not a physical parameter.
B lift contains A; no independent regional Hamiltonian or new local control
permission is asserted. An unchanged edge uses its original derivative.
Original geometry weights and lattice spacing are not recomputed at the cut.
All original multiplication potentials, CAR hopping and masses pull back
through AB. There is no new boundary fermion or new tensor sign convention.
Closed form domain is J D(q); matching Hilbert space is the unitary image.
Full unbounded-domain and thermal results follow by form transport. They
are not inferred from a finite PW matrix or claimed on all extended states.
Source differentiation uses a fixed background Hilbert identification,
geometry-independent J and the inherited differentiable common form domain.
No bounded stress operator on every state, quantum geometry, or Einstein
constraint follows. Geometry-dependent splitting must also be differentiated.
Ordinary thermal trace is over the physical matched space. Parity sectors
are retained according to 616; no hidden supertrace or discarded sector.
Instrument histories and arbitrary reference are intertwined but their
physical controllability/locality is not furnished by an abstract encoding.
Numerical electric fixture has three original weak PW sectors, fixed right
indices, before endpoint Gauss. The 4096x8 map tests genuine off-diagonal
geometry. Full physical equivalence is analytic, not this fixture's spectrum.
Numerical interaction fixture uses actual L representation in one fixed
spin component and four endpoint modes; it is not a full fermion spectrum.
Doubled energy and omitted shear are two specific failed prescriptions.
The proof is a connection to mature gauge gluing, not a novel general theorem
about gauge factorization, nor a claim of physical continuum refinement.
Independent fields, full chirality, spacetime dimension, GR, and predictions
remain open. Next entry requires historical audit before any new round.
"""
for n in ("research_note_617.md","joint_region_energy_gluing.py",
          "joint_region_energy_gluing_results.json","unified_physics_condition_ledger_617.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round617_drafts/final_review.txt",review)
print("617 review and publication prepared")

