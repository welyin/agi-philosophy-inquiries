"""Prepare a scope-reviewed 613 without rewriting any frozen evidence."""
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

mapping={"joint_overlap_time_spectrum":"joint_hamiltonian_gauge_contract",
         "611":"612","612":"613","613":"614",
         "3078":"3081","3081":"3084",
         "1145":"1148","1148":"1151",
         "1870":"1877","1877":"1884"}
verify=replace((HERE/"verify_round612.py").read_text("utf8"),mapping)
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==12")
verify=verify.replace("Verify actual overlap time spectrum and the joint thermal contract.",
                      "Verify spatial Hamiltonian, original gauge algebra and common sources.")
anchor="    links=0\n"
assert anchor in verify
verify=verify.replace(anchor,
    "    ledger=(HERE/'unified_physics_condition_ledger_613.md').read_text('utf8')\n"
    "    ids=re.findall(r'^\\|C(\\d\\d) ',ledger,re.M)\n"
    "    assert ids==[f'{i:02d}' for i in range(1,28)]\n"+anchor)
assert "cognitive_foundation_bridge_605_navigation.json" in verify
write("verify_round613.py",verify)
publish=replace((HERE/"publish_round612.py").read_text("utf8"),
                mapping|{"## 258.":"## 259.","## 163.":"## 164."})
start=publish.index("summary=");end=publish.index("planned={}",start)
header="""summary=('**第613轮完成：** [同一Hamiltonian的手征、规范与几何来源条件]({p}research_note_613.md)'
         '成熟连续时间候选可保有限热态，但朴素手征筛选破坏原弱代数和U(1)周期；'
         '同谱海能与压力不能仅由单点Λ抵消共同匹配。'
         '三组、十二式通过，最新613／3084，1151份编号科学文件、1884份保护证据。'
         '[核验]({p}research_round_613_checks.json)、[C01—C27更新总账]({p}unified_physics_condition_ledger_613.md)。'
         '仅限指定向量样候选和接法；真实手征、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（613后，优先于下方历史安排）：** 先合并物理成立条件，认知系统设计后置。'
       '接[614原一代物质与成熟手征构造]({p}round614_drafts/STATUS.md)，'
       '核原表示、商群与Yukawa／Majorana能否共用候选输入，逐项记录新增场及来源；'
       '不重复精度优化，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace("同一时间谱、正演化与热态合同","同一Hamiltonian、规范及几何来源")
publish=publish.replace("正演化与原热态条件须共同匹配","手征、规范和来源不能分别签收")
publish=publish.replace("完整时间关联、正演化与热态条件的共同匹配",
                        "同一Hamiltonian的手征、规范与几何来源条件")
assert "hashlib.sha256" in publish
write("publish_round613.py",publish)
write("round613_drafts/research_note_613_draft.md",
      (HERE/"research_note_613.md").read_text("utf8"))
review="""613 primary-agent review. No new independent agent review.
Previous goal turn: progress, verified 612 science and navigation publication.
CHN identities and the known charge-quantization issue are attributed to
primary literature; not counted as discoveries. Full normalization follows
the existing m0=1 Hermitian Euclidean gamma convention.
This is a spatial/continuous-time alternative, not the exact 612 time kernel.
The inherited hopping block, h^2 and q5 identities are checked independently.
Absence of extra zeroes follows analytically from sine=0 and W<0, not merely
the numerical corner scan. One Dirac still contains both chiralities.
P=(1+-q5)/2 is the explicitly tested naive filter. Its nonidempotency gives
the original weak algebra defect. The old U1 periods are inherited; failures
are operator identities, not anomaly coefficients. Some charges/corners are
exceptions, and are not hidden.
Normalizing q5 repairs a static block but lacks a continuous corner extension;
this is not a no-go claim for all chiral regularizations or all nonlocal ones.
The thermal calculation uses 16 vectorlike copies and 64 modes per node,
not a claimed completion of the original 32-mode chiral model.
Negative one-particle eigenvalues do not make finite CAR energy unbounded.
The fixed-grid scale derivative keeps beta and labels fixed. Sea subtraction
and its geometry derivative are handled together.
Pressure of this spatial regulator is not physical Lorentz-invariant vacuum
pressure. The single cosmological counterterm result has a stated scope;
general covariant counterterms and the quantum GR connection remain open.
Old 578-579 source logic and 602 thermodynamic identities are explicitly
reused. Actual interacting Yukawa/Majorana, Gauss thermal state, original
full chiral measure and spacetime generation are not claimed.
Ledger covers C01-C27 without turning audit headings into new axioms.
The large 2026 slides were not obtained as trusted text and are not used.
Current priority remains joint conditions; cognitive devices are deferred.
"""
for n in ("research_note_613.md","joint_hamiltonian_gauge_contract.py",
          "joint_hamiltonian_gauge_contract_results.json",
          "unified_physics_condition_ledger_613.md"):
    review+=n+" sha256="+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+"\n"
write("round613_drafts/final_review.txt",review)
print("613 review and guarded publication prepared")
