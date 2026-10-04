"""Prepare 616 primary review and append-only publication."""
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

mapping={"joint_subgroup_measure_source":"joint_spin_thermal_source",
         "614":"615","615":"616","616":"617",
         "3087":"3090","3090":"3092","1154":"1157","1157":"1160",
         "1891":"1900","1900":"1907"}
verify=replace((HERE/"verify_round615.py").read_text("utf8"),mapping)
verify=verify.replace("Verify the exact finite subgroup measure and common source.",
                      "Verify original neutral spin obstruction and thermal source contract.")
verify=verify.replace(",'round616_drafts/measure_source_entry.md','round616_drafts/probe_auxiliary.py'","")
verify=verify.replace("==(3,0,0)","==(2,0,0)").replace("fresh_tests=dict(run=3,","fresh_tests=dict(run=2,")
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==11")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round616.py",verify)
publish=replace((HERE/"publish_round615.py").read_text("utf8"),
                mapping|{"## 261.":"## 262.","## 166.":"## 167."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第616轮完成：** [原中性费米子、旋量结构与热态边界的共同条件]({p}research_note_616.md)'
         '给定定向全局双曲3+1分支可合并spin存在条件；原中性费米子阻止规范补偿。'
         '普通Gauss热态固定时间闭合及来源，双宇称周期超迹不能替代正态。'
         '两组、十一式通过，最新616／3092，1160份编号科学文件、1907份保护证据。'
         '[核验]({p}research_round_616_checks.json)、[条件账]({p}unified_physics_condition_ledger_616.md)。'
         '维数、全局几何、空间spin选择及图到连续映射仍开放，无新增独立代理审查。')
order=('**当前执行顺序（616后，优先于下方历史安排）：** 继续整合共同条件，认知设计后置。'
       '接[617原费米区域、Gauss拼接与共同态]({p}round617_drafts/STATUS.md)，'
       '回查旧区域结果，核同一群、态与来源的真实拼接；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("同一测度、物理权重与共同来源","原物种、spin存在与热态来源")
publish=publish.replace("有限测度正性与全局物理前提","条件性spin存在不选择维数或空间边界")
publish=publish.replace("原子群的精确辅助测度与共同响应","原中性费米子、旋量结构与热态边界的共同条件")
assert "hashlib.sha256" in publish
write("publish_round616.py",publish)
write("round616_drafts/research_note_616_draft.md",
      (HERE/"research_note_616.md").read_text("utf8"))
review="""616 primary-agent review. No new independent agent review.
Geometry: smooth oriented time-oriented boundaryless globally hyperbolic
four-dimensional spacetime is explicitly assumed. Bernal-Sanchez splitting
and orientable 3-manifold parallelizability imply spin existence, not the
choice of dimension, topology, preferred framing, connection or spin lift.
Spin choices form a torsor, not canonically the group H1 itself.
Eight choices on the marked fixed T3 do not classify large-diffeomorphism orbits.
The original nu is a gauge singlet both before and after the 614 dictionary.
No element of the original quotient group acts as minus identity on it.
Thus a diagonal spin-gauge central quotient cannot compensate its spin sign.
No claim excludes altered matter, composites or enlarged gauge groups.
The 598 full Hamiltonian preserves parity; the 603 physical thermal trace
is finite. Both physical parity sectors are nonzero: gauge-invariant bosonic
states times vacuum, and a neutral nu creation. Heat is strictly positive on
each nonzero sector. Therefore both Z_even and Z_odd are strictly positive.
If the supertrace is nonzero, one positive parity projector has negative
normalized supertrace. Zero supertrace is not normalizable. Parity projectors
are even observables; superselection alone does not discard a sector.
An explicitly fixed-parity reference is a different legal choice and uses
the AP/P combination. That exception is stated, not silently ruled out.
Coherent-state trace gluing applies to the specified Hamiltonian; the reverse
implication from an arbitrary AP Euclidean measure to positive reconstruction
is not claimed. Thermal time circle is not a Lorentzian closed time curve.
The numerical mass block comes directly from old mass_matrices at phi=(0,0,0,0,s),
retaining old complex Y_s and F=M0^2-s^2/6. No kinetic reduction is asserted.
Four-state spectrum and source are exact for that mass block. Full Gauss
nonpositivity relies on the analytic argument, not the small fixture.
No source counterterm or zero-energy subtraction was changed. s=0 is handled
without taking the logarithm of zero. s>0 is stated for mass derivative.
Global geometry and fixed graph thermal branches still need a common scale
mapping; no generic nonstationary background is assigned an equilibrium.
General measure, full chiral dynamics, GR and predictions remain unproved.
"""
for n in ("research_note_616.md","joint_spin_thermal_source.py",
          "joint_spin_thermal_source_results.json","unified_physics_condition_ledger_616.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round616_drafts/final_review.txt",review)
print("616 review and publication prepared")
