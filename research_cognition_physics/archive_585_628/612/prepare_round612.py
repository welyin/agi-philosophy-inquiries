"""Prepare current primary review and guarded publication; preserve old rounds."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():
            p=r"(?<!\d)"+p+r"(?!\d)"
        patterns.append(p)
    return re.sub("|".join(patterns),lambda m:mapping[m.group()],text)

def write(name,text):
    path=HERE/name
    path.parent.mkdir(exist_ok=True)
    with path.open("x",encoding="utf8",newline="\n") as f:
        f.write(text)

mapping={
    "joint_original_mass_spinor_bridge":"joint_overlap_time_spectrum",
    "610":"611","611":"612","612":"613",
    "3075":"3078","3078":"3081",
    "1142":"1145","1145":"1148",
    "1863":"1870","1870":"1877",
}
verify=replace((HERE/"verify_round611.py").read_text("utf8"),mapping)
verify=verify.replace("text['display_formulas']==12","text['display_formulas']==14")
verify=verify.replace("Verify projected process completion, full histories and common sources.",
                      "Verify actual overlap time spectrum and the joint thermal contract.")
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round612.py",verify)

publish=replace((HERE/"publish_round611.py").read_text("utf8"),
                mapping|{"## 257.":"## 258.","## 162.":"## 163."})
publish=re.sub(r"navigation_before_round612_\d+",
               "navigation_before_round612_20261001",publish)
start=publish.index("summary=")
end=publish.index("planned={}",start)
header="""summary=('**第612轮完成：** [完整时间关联、正演化与热态条件的共同匹配]({p}research_note_612.md)'
         '继承自由overlap核的非零动量关联含正连续谱；可正演化重建，'
         '却不能同时精确全部时间匹配原固定图的紧谱与有限热迹。'
         '三组、十四式通过，最新612／3081，1148份编号科学文件、1877份保护证据。'
         '[核验]({p}research_round_612_checks.json)、[条件账]({p}unified_physics_condition_ledger_612.md)。'
         '仅排除明示自由强匹配分支；有限窗口、完整规范及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（612后，优先于下方历史安排）：** 先整合条件，认知系统设计后置。'
       '接[613共同匹配范围与剩余条件整合]({p}round613_drafts/STATUS.md)，'
       '回填598—612结果，比较量子过程、物质、热态和几何来源的同一尺度合同；'
       '不继续向记忆或辅助硬件设计下钻，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("原物种质量、旋量及来源共同字典","同一时间谱、正演化与热态合同")
publish=publish.replace("静止质量不替代完整实时间过程","正演化与原热态条件须共同匹配")
publish=publish.replace("原物种质量、旋量投影与共同来源字典",
                        "完整时间关联、正演化与热态条件的共同匹配")
assert "hashlib.sha256" in publish and "sha257" not in publish
write("publish_round612.py",publish)
write("round612_drafts/research_note_612_draft.md",
      (HERE/"research_note_612.md").read_text("utf8"))

review="""612 primary-agent final review; no independent agent review.
Contract is explicit: inherited free massless m0=1 kernel, finite periodic space,
fixed cutoff, infinite Euclidean time, fixed normal ground-state observables.
This does not rule out the original nonzero interacting parameter point unless
its proposed dictionary is required to include this free slice.
Gamma conventions, normalization and actual inherited hopping block agree.
Direct inverse and specified free Weyl block both give the tested scalar.
The inverse propagator spectral density is not confused with the action kernel.
Pole residue, branch discontinuity and positivity follow analytically.
The zero spatial momentum exception and equal-time contact limitation are stated.
The second diagnostic Bloch momentum is not claimed to lie on the L=4 grid.
All-order Hankel positivity comes from interval support, not a small matrix test.
Compact-resolvent obstruction uses finite positive moment uniqueness and a
normal ground-state vector, not the incorrect claim that all infinite spaces
have discrete spectra. The earlier 596 theorem is explicitly inherited.
The cyclic scalar construction is not claimed to reconstruct all CAR or gauge
correlations. Exact spectral inclusion, unlike a finite-window fit, enforces
infinite heat trace through infinitely many orthogonal finite-energy vectors.
Field spectral weights do not count states in a partition function.
Finite Euclidean path integrals are not automatically canonical thermal traces.
Quadrature errors and energy-insertion tests are illustrative validation; they
are not the proofs of infinite rank or divergent trace.
2010 free/Weyl and specific non-gauge Yukawa results retain their scope.
The 2026 Spin(10) item has been read only as an abstract, not used as a theorem.
Integration ledger retains quantum, matter, state, scale and source conditions.
Cognitive device design stays deferred. Goal unchanged; gravity still open.
"""
for name in ("research_note_612.md","joint_overlap_time_spectrum.py",
             "joint_overlap_time_spectrum_results.json",
             "unified_physics_condition_ledger_612.md"):
    review+=name+" sha256="+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+"\n"
write("round612_drafts/final_review.txt",review)
print("612 report, review and guarded publication helpers prepared")
