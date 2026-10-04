"""Prepare 766 artifacts without replacing previously frozen work."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(name, content):
    p = HERE/name
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("x", encoding="utf8", newline="\n") as f:
        f.write(content)


lines = (HERE/"unified_physics_condition_ledger_765.md").read_text("utf8").splitlines()
lines[0] = "# 联合条件总账：766精确规范关系与Hadamard二点候选"
lines[2] = "2026-10-04。接[765全账](unified_physics_condition_ledger_765.md)，回填[766报告](research_note_766.md)。[结果](joint_microlocal_gauge_projection_results.json)、[核验](research_round_766_checks.json)。同一E自由分支中短距离与精确规范/CCR已相容；正性尚未证明。"
updates = {
    "C01": "766同一完整D1有精确规范与CCR的Hadamard二点候选；未经物理正性不能签收为态",
    "C04": "766全阶频率分解与原D1K=KD0共同给平滑交织缺陷及其精确修正",
    "C09": "766利用753原共同参考及内部配置的无连续稳定子，证明Cauchy规范Gram算符核零",
    "C11": "766新椭圆Cauchy规范逆及频率相容投影成立；线性规范身份不等于全量子约束完成",
    "C19": "766精确规范/对易/短距离可共同实现；剩同一物理商上正性及允许平滑修正",
    "C22": "766尚未正态化或完成全部门二次来源；735单费米圈及原任务映射的范围保持"}
seen = set()
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith("|"+key+" "):
            parts = line.split("|"); parts[2] += "；"+value
            lines[i] = "|".join(parts); seen.add(key)
assert seen == set(updates)
write("unified_physics_condition_ledger_766.md", "\n".join(lines)+"""

## 766短距离与精确规范身份的共同候选

由765全耦合正常双曲算符的一般Hadamard投影，构造频率适配的辅助正配对。完整归一Cauchy规范映射是单射主符号的一级过定椭圆算符；753无连续联合稳定子消去其全局核。Gram逆给B∈Ψ^-1及Π∈Ψ^0，精确去规范并保物理电荷，且与频率投影只差平滑交换子。

成熟Cauchy规范修正公式的所有复合都有明确有限阶，故修正为平滑。修正后同一完整D1二点候选精确保规范、Hermitian及全对易核，也保Hadamard波前。这不是正态存在性：正性必须在同一物理商上检验。有限代数负方差见证只排除“规范修复自动保正”的错误推论。

766没有新增作用、物种、参考场或普遍认知公理；没有把原来源右逆未经证明称为任意拟微分阶。已有无稳定子结果在当前共同构造中实质消除了全局核条件。物理Hadamard正态、全来源、实际记录及原Q连续映射保持开放。

## 下一项：共同物理正性及平滑修正

接767，处理同一候选在N/imK上的正性。必须保精确规范与对易关系，仅做不改变短距离奇性的修正。不得以紧初片自动推断有限秩修复足够；须核物理正频椭圆估计、低频谱及无限秩平滑余项。旧空间、604、649/699和统一目标保持。
""")
write("round766_drafts/research_note_766_draft.md", (HERE/"research_note_766.md").read_text("utf8"))
write("round767_drafts/STATUS.md", """# 第767轮入口：同一物理商上的正性及保规范平滑修正

接[766](../research_note_766.md)、[全账](../unified_physics_condition_ledger_766.md)。

1. 766已给同一完整耦合D1的Hadamard二点候选，精确保Hermitian、规范和完整CCR；未证明物理正性。不要回退到另一份不具短距离证书的Gaussian态。
2. Cauchy规范映射K为Ψ^1，753无连续联合稳定子给核零；频率适配辅助H、Gram逆B、Π和成熟平滑修正均已有映射。不要再扫投影矩阵或主符号秩。
3. 正性准确对象是 ±(Πz,Q1 c1±Πz)，z∈ker K^ddagger。辅助H正不承担它。应先核物理频率子空间的椭圆正性，再控制平滑/低频缺陷。
4. 734/766入口已说明平滑未必有限秩，有限负谱修复需证书。保留允许无限秩、快速衰减的共同平滑修正分支，但必须证明规范、CCR和实结构共同保持。
5. 即使正态构造成立，全部门局部Ward、同阶玻色/费米来源、原Q连续映射及原实际记录/装置仍在总验收单；不得将自由阶段当作全量子引力。
6. 旧空间382—386、425、522—523以及604、649/699边界不改。目标保持，不新增应用任务、调度或图像。
""")
write("round766_drafts/scope_and_dedup_review.md", """# 766范围与去重复核

- 上一轮765及766入口为有效进展。本轮回读导航、最新结果、731/754、753无稳定子证明；无运行中Python实验。
- 734的平滑无限秩结论及766入口校准复用，不增加科学轮次。数值三组只核本轮新公式。
- 原完整D1/D0形式自伴且正常双曲，经相容联络运输和半密度化适配一般不定纤维Hadamard投影；未搬用真空引力正性部分。
- 归一Cauchy资料使K为Ψ^1；空协向量下规范生成符号单射，给过定椭圆性；全核由原753连续联合稳定子零证明消除。数值未负责全局零模。
- sharp为辅助H伴随，ddagger为真实Q伴随；不能混用。新B不是未经阶数审查的旧来源R。
- 频率相容由全阶交织缺陷平滑及有限阶椭圆逆推出。修复采用完整补项，未只压缩导致丢失全CCR。
- 精确规范Hadamard二点候选仍未正。修正c不必投影，物理负方差校准仅反驳方法捷径，不反驳当前模型存在性。
- 未重新选择背景、作用或输入维数/群，没有原Q→E、实际L_r(s)、全来源或全量子约束完成声明。
- 主代理人工检查算符对象、支撑、逆、所有伴随及准确文献前提；无独立代理或图像检查。
""")
write("round766_drafts/literature_scope_audit.json", json.dumps(dict(
    sources=[
        dict(url="https://arxiv.org/html/2311.11043",
             checked="Sections3.1-3.2, formulas3.30-3.33: general indefinite Hermitian normally hyperbolic setting; section3.3.2 gauge-sign smoothing.",
             actual_mapping="Full coupled D_i moved by normal metric-connection transport and density conjugation to constant Hermitian fiber form and scalar elliptic spatial principal part.",
             not_imported="Vacuum Einstein requirement and trace-reversal/TT-synchronous physical positivity proofs in later sections."),
        dict(url="https://arxiv.org/html/1403.7153",
             checked="Theorem3.17, Remark3.18: Cauchy gauge repair, smooth correction, exact CCR and explicit separate positivity assumption.",
             actual_mapping="Compact normalized Cauchy K, newly proved elliptic inverse B and projection Pi. Operators have finite Psi order, ensuring correction smoothing.",
             not_imported="Pure YM Hadamard-state existence theorem or unconditional physical positivity.")],
    project_increment="Original irreducible joint background supplies the global kernel condition and a frequency-compatible elliptic gauge section; all exact gauge/CCR identities and Hadamard singularities have a common candidate.",
    physical_positive_Hadamard_state_proven=False,
    all_sector_Ward_or_original_task_mapping_proven=False), ensure_ascii=False, indent=2)+"\n")
files = ("research_note_766.md", "joint_microlocal_gauge_projection.py",
         "joint_microlocal_gauge_projection_results.json", "unified_physics_condition_ledger_766.md")
write("round766_drafts/final_review.txt",
      "Primary reviewer only. Checked full on-shell operators, Cauchy normalization orders, "
      "injective principal symbol and use of original joint-background stabilizer proof. "
      "Auxiliary positive H and indefinite physical Q adjoints are distinct. "
      "Gram inverse is a new Psi^-2 elliptic inverse on the fixed irreducible background. "
      "All smoothing compositions have established finite orders. "
      "Both compensating terms in the GW repair are retained; exact sum, adjoint and gauge identities "
      "do not assert idempotency or positivity. A finite negative physical variance is only a "
      "scope witness for the invalid positivity shortcut. No physical Hadamard state or all-sector "
      "renormalized source is certified.\n"+
      "\n".join(n+" "+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+"\n")

mapping = {"765": "766", "764": "765", "3679": "3691", "3691": "3705",
           "3479": "3482", "3482": "3485", "1607": "1610",
           "joint_covariant_gauge_complex": "joint_microlocal_gauge_projection"}
pattern = re.compile("|".join(map(re.escape, sorted(mapping, key=len, reverse=True))))


def convert(text):
    return pattern.sub(lambda m: mapping[m[0]], text)


verify = convert((HERE/"verify_round765.py").read_text("utf8"))
a, b = verify.index("    extra = ("), verify.index("    new = ")
verify = verify[:a]+"""    extra = ('unified_physics_condition_ledger_766.md',
             'round766_drafts/research_note_766_draft.md',
             'round766_drafts/final_review.txt',
             'round766_drafts/literature_scope_audit.json',
             'round766_drafts/scope_and_dedup_review.md',
             'round767_drafts/STATUS.md',
             'round766_drafts/physical_hadamard_entry.md',
             'round766_drafts/smooth_positivity_calibration.py',
             'round766_drafts/smooth_positivity_calibration_results.json',
             'round766_drafts/physical_hadamard_entry_checks.json',
             'round766_drafts/entry_postpublication_checks.json')
"""+verify[b:]
verify = verify.replace("checks['display_formulas'] == 14", "checks['display_formulas'] == 16")
a, b = verify.index("                complete_linear_covariant_gauge_complex_constructed="), verify.index("                primary_code_and_note_review_completed=")
verify = verify[:a]+"""                actual_coupled_Cauchy_gauge_ellipticity_and_injectivity_proven=True,
                exact_frequency_compatible_gauge_projection_constructed=True,
                joint_Hadamard_gauge_CCR_candidate_constructed=True,
                physical_positive_Hadamard_state_proven=False,
                numerical_checks_only_finite_algebraic_calibrations=True,
                all_sector_quantum_Ward_completed=False,
                original_Q_to_E_process_equivalence_proven=False,
"""+verify[b:]
write("verify_round766.py", verify)

publish = convert((HERE/"publish_round765.py").read_text("utf8"))
a, b = publish.index("summary = "), publish.index("planned = ")
publish = publish[:a]+"""summary = '**第766轮完成（精确规范与短距离共同候选）：** [频率相容的规范投影]({p}research_note_766.md)利用原背景无连续联合稳定子，构造真实椭圆规范逆；同一完整自由场二点候选精确保规范、CCR及Hadamard奇性。物理正性仍未证，不能称已得量子态。三组、十六式通过，最新766／3485，1610份编号科学文件、3705份保护证据。[核验]({p}research_round_766_checks.json)、[条件账]({p}unified_physics_condition_ledger_766.md)。'
order = '**当前执行顺序（766后，优先于下方历史安排）：** 接[767同一物理商的正性与平滑修正]({p}round767_drafts/STATUS.md)，保已得全部精确身份，核物理椭圆正性与平滑缺陷；不延伸投影参数。[范围审计]({p}round766_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
"""+publish[b:]
publish = publish.replace("第766轮完成（共同因果量子对象）", "第766轮完成（精确规范与短距离共同候选）")
publish = publish.replace("## 411.", "## 412.").replace("## 316.", "## 317.")
publish = publish.replace("完整耦合线性规范复形与因果物理代数", "频率相容的规范投影与短距离共同候选")
publish = publish.replace("旧空间合同保持，共同因果代数与物理短距离准备", "旧空间合同保持，精确规范与物理正性")
publish = publish.replace("完整耦合线性规范复形", "频率相容的规范投影")
publish = publish.replace("next_round=766", "next_round=767")
write("publish_round766.py", publish)
post = convert((HERE/"postcheck_round765.py").read_text("utf8"))
post = post.replace("range(584, 766)", "range(584, 767)")
post = post.replace("第766轮完成（共同因果量子对象）", "第766轮完成（精确规范与短距离共同候选）")
write("postcheck_round766.py", post)
print("766 artifacts and verification/publication tools prepared.")
