"""Prepare 614 primary review, verification and append-only publication."""
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

mapping={"joint_hamiltonian_gauge_contract":"joint_spinor_subgroup_mass",
         "612":"613","613":"614","614":"615",
         "3081":"3084","3084":"3087",
         "1148":"1151","1151":"1154",
         "1877":"1884","1884":"1891"}
verify=replace((HERE/"verify_round613.py").read_text("utf8"),mapping)
verify=verify.replace("Verify spatial Hamiltonian, original gauge algebra and common sources.",
                      "Verify the original subgroup and CAR mass dictionary.")
begin=verify.index("    ledger=")
end=verify.index("    links=0",begin)
verify=verify[:begin]+verify[end:]
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round614.py",verify)
publish=replace((HERE/"publish_round613.py").read_text("utf8"),
                mapping|{"## 259.":"## 260.","## 164.":"## 165."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第614轮完成：** [原物质、全局商群与成熟手征表示的共同字典]({p}research_note_614.md)'
         '原16内部表示、Z₆商、32模式CAR质量及来源可共同映入候选容器的子群，'
         '无需扩大物理规范群；若扩大为全Spin(10)，原实singlet Majorana接法失败。'
         '三组、十二式通过，最新614／3087，1154份编号科学文件、1891份保护证据。'
         '[核验]({p}research_round_614_checks.json)、[条件账]({p}unified_physics_condition_ledger_614.md)。'
         '有限质量字典不是完整手征测度或时间重建，局域性、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（614后，优先于下方历史安排）：** 继续整合原群、量子过程和共同来源，认知设计后置。'
       '接[615原子群的Weyl测度与共同来源]({p}round615_drafts/STATUS.md)，'
       '核配置限制、测度变化与同一来源的真实连接，保留成熟候选的局域性和重建边界；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("同一Hamiltonian、规范及几何来源","原群与原质量的共同表示字典")
publish=publish.replace("手征、规范和来源不能分别签收","表示容器不等于物理规范群扩大")
publish=publish.replace("同一Hamiltonian的手征、规范与几何来源条件",
                        "原物质、全局商群与成熟手征表示的共同字典")
assert "hashlib.sha256" in publish
write("publish_round614.py",publish)
write("round614_drafts/research_note_614_draft.md",
      (HERE/"research_note_614.md").read_text("utf8"))
review="""614 primary-agent review. No new independent agent review.
Source 1710.11618 explicitly permits subgroup restriction and discusses SM
with right-handed neutrinos; its measure-locality issue remains within the
paper's stated scope. No claim of solving that open problem or that no later
results exist. The inaccessible 2026 slides are not used as evidence.
Original right-handed entries are charge-conjugated to all-left convention.
The SU5 map has precisely the inherited Z6 kernel, not just a Lie-algebra match.
The wedge^4 module is faithful on SU5, so the representation adds no kernel.
The explicit signed basis dictionary includes the weak epsilon intertwiner.
The basis Hodge-complement map is linear as defined, not an implicit assertion
that a complex Hodge star on arbitrary vectors is complex-linear.
Internal 16 and the five auxiliary Clifford oscillators are not spacetime
dimensions or extra physical modes. 32 original CAR modes remain 32.
The finite particle-hole map preserves CAR and the actual original BdG
matrix. All original complex Y and F dependence remain, including the
conjugation/sign of the neutral Majorana pair. Zero normal-order constant is
justified for this off-diagonal mass block; not generalized to full hopping.
Same fixed transform preserves the scalar-source derivatives. Moving Weyl
fibres, gauge measure and actual full time evolution are not reconstructed.
Full Spin10 gauging is a separate branch. With a full-group singlet scalar,
the Majorana coefficient fails the extra Cartan identity. Its infinitesimal
stabilizer is su5 (24 directions), not a claimed complete global stabilizer.
The 126 dimension is standard representation theory, computed exactly by
the Weyl formula. A scalar completion requirement is limited to unbroken
full Spin10 with linear elementary fields and renormalizable Yukawa terms.
No 126 or extra gauge bosons are required merely by the subgroup carrier.
Claims do not cover auxiliary-measure resource cost, locality, positive
Hamiltonian reconstruction, original interacting Gauss thermal state or GR.
"""
for n in ("research_note_614.md","joint_spinor_subgroup_mass.py",
          "joint_spinor_subgroup_mass_results.json","unified_physics_condition_ledger_614.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round614_drafts/final_review.txt",review)
print("614 review and publication prepared")
