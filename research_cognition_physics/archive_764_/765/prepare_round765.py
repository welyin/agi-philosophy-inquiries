"""Exclusive preparation of 765 reports, ledger and publication tools."""
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


lines = (HERE/"unified_physics_condition_ledger_764.md").read_text("utf8").splitlines()
lines[0] = "# 联合条件总账：765完整耦合规范复形与共同因果物理代数"
lines[2] = "2026-10-04。接[764全账](unified_physics_condition_ledger_764.md)，回填[765报告](research_note_765.md)。[结果](joint_covariant_gauge_complex_results.json)、[核验](research_round_765_checks.json)。原E线性量子化分支继续，不等同原Q连续匹配。"
updates = {
    "C01": "765原耦合协变物理代数与764 canonical物理空间辛同构，正则正态可沿此映射",
    "C03": "765紧支规范不变线性观测有共同因果代数；原实际L_r(s)及全非线性关系记录未接通",
    "C04": "765原完整Hessian与完整规范生成元组成正常双曲规范复形",
    "C05": "765在输入的同一g上证明线性合法观测微因果性，未生成g或维数",
    "C09": "765规范条件含实际曲率/标量零阶混合；不得照搬独立de Donder加Lorenz伴随",
    "C11": "765无穷小线性规范商连接成立；全量子约束及大规范识别仍须原范围处理",
    "C19": "765正则态与物理Hadamard准备严格分开，真空GR/独立YM存在定理不直接覆盖原耦合背景",
    "C22": "765同阶玻色/费米来源对象已明确，短距离态及共同局部量子Ward尚未完成"}
seen = set()
for i, line in enumerate(lines):
    for key, value in updates.items():
        if line.startswith("|"+key+" "):
            parts = line.split("|"); parts[2] += "；"+value
            lines[i] = "|".join(parts); seen.add(key)
assert seen == set(updates)
write("unified_physics_condition_ledger_765.md", "\n".join(lines)+"""

## 765共同因果量子对象的线性连接

在原753在壳背景、正F、573 Einstein框架作用及紧T³短时发展上，完整Hessian P和协变化规范生成元K满足PK=0。选按原动能归一的非退化纤维配对，构造D1=P-KK†与D0=-K†K；全部二阶主部为同一g的波算子，低阶物质混合全部保留。先进/推迟Green算符的交织给因果物理可观测商。

原canonical初值与协变解的同一作用辛流、线性规范唯一性和紧初片共同给物理辛同构。764正则正态可拉回；空间样分离的紧支合法线性观测对易。不把非局域截面坐标、h(x)、a(x)或原固定图仪器直接当成这些局域观测。

这是原输入Einstein作用的条件连接，不是维数、共同g或量子引力从认知生成。仍缺当前耦合系统的物理Hadamard准备、同阶全来源与Ward、原Q映射及自主实际记录。正则态、局部parametrix和正物理Hadamard态三者分开。735单费米圈范围不扩大。

## 下一项：同一短距离准备的规范与正性

接766，优先核高频物理投影、精确约束与正性可否同时实现。现有R提供精确初始约化，但任意辅助Gaussian宽度没有正确短距离证书；成熟真空GR/纯YM结果不可直接相乘。旧空间及604、649/699与统一目标保持。
""")
write("round765_drafts/research_note_765_draft.md", (HERE/"research_note_765.md").read_text("utf8"))
write("round766_drafts/STATUS.md", """# 第766轮入口：共同物理Hadamard准备的规范相容与正性

接[765](../research_note_765.md)、[条件账](../unified_physics_condition_ledger_765.md)、[H1—H3](../round755_drafts/cognitive_joint_candidate_working_report.md)。

1. 原完整协变P/K/D1/D0已在E线性分支组成正常双曲规范复形，因果物理代数与764 canonical约化辛同构。不要重复主部、自由度或协向量扫描。
2. 正则正态已存在；短距离parametrix的存在不等于物理正Hadamard态。当前要共同满足方程、正确对易核、规范商良定义、正性及微局部谱条件。
3. 回查Gérard2311.11043的Einstein真空条件、GW1403.7153纯YM背景条件及WZ1407.8079的独立正性/规范前提。优先构造/适配当前全耦合对象，不将无法直接套定理称为不存在。
4. 764 R和P给精确Cauchy约化；须核其与全阶高频分解的相容，不只核首项正性。允许低频有限秩修正，但须证明足够且不改变物理对易/规范或旧任务合同。
5. 同阶完整来源含玻色约束Hessian及费米来源。735只是一费米圈；物理态之后仍须全部门局部Ward和原来源/记录字典。
6. 原Q连续匹配、实际原L_r(s)、注能和自主准备保持联合验收项。旧空间382—386、425、522—523以及604、649/699不改；目标保持。
""")
write("round765_drafts/scope_and_dedup_review.md", """# 765范围与去重复核

- 573已证明原经典短时发展，不重复；764是全初片canonical辛约化，本轮新接同一完整Hessian的协变因果物理代数。
- 采用原Einstein框架、二导数最小动能、正F和在壳背景；不直接用真空引力Hessian，不加入752反例中的梯度耦合。
- h不迹反转；V度规为DeWitt配对/(4κ)、w g逆和K_AB，W为g/(2κ)与w。伴随零阶混合按完整K逐项计算，负导数符号明确。
- T=-K映射Hack—Schenkel一般框架；R_HS=-D0、Q_HS=D0均双曲。本文P与764截面P不同，报告以文字区分。
- 闭T³上所有光滑规范参数空间紧支；原乘子/数据规范匹配及作用Green恒等式给canonical桥。不同全局支撑类、边界区域因子及大规范识别不被顺带证明。
- 数值校准完整二阶时空主部，但未离散低阶全Hessian；原背景零阶检查只含全部标量块和颜色磁曲率，明确不含全部电场/电弱曲率。
- 真空GR、独立YM、单KG文献只给方法/比较，不代替当前全耦合物理Hadamard构造；BRST正性亦单列。
- 规范商线性观测不是h(x)或实际原仪器。原Q→E、全来源、相互作用和自主记录仍未签收。
- 初次程序失败是对非零块不必要地要求>0.1，改为非零容差；修复前代码保留。764发布后出现历史字内多空格，保留观察版本后精确恢复冻结SHA；无科学内容更改。
- 三组新校准、十四式；未做独立代理或图像检查。目标、空间旧结论、604及649/699保持。
""")
sources = [
    dict(url="https://arxiv.org/html/1205.3484",
         checked="Definition3.4, Theorem3.12, Theorem5.2 and Green/Cauchy pairing in section5.",
         mapping="Full coupled P,K, T=-K; normally hyperbolic D1,D0 checked in note765. Compact T3 removes the two spatial-support gauge-image distinction.",
         not_imported="Vacuum gravity example, universal nondegeneracy or Hadamard positivity."),
    dict(url="https://arxiv.org/html/2311.11043",
         checked="Section1.1 and section2 require Ric=Lambda g and compact Cauchy surfaces.",
         use="Pure-gravity comparison/method only; actual matter-supported background not automatically covered."),
    dict(url="https://arxiv.org/html/1403.7153",
         checked="Theorems1.1/1.2; Hypothesis2.3; conditions pos/g.i./microlocal.",
         use="Independent linear Yang-Mills about an on-shell pure YM connection. Not a full quantum metric/Higgs theorem."),
    dict(url="https://arxiv.org/html/1403.3957",
         checked="Section2 quantizes coupled Einstein-Klein-Gordon on arbitrary on-shell background.",
         use="Method comparison; not a ready-made Hadamard existence theorem for the full project model."),
    dict(url="https://arxiv.org/html/1407.8079",
         checked="Section3, Lemma3.2, Remark3.3 and Proposition3.4.",
         use="Short-distance, gauge compatibility and physical positivity remain separate; modulo-gauge relations insufficient for automatically extending all nonlinear observables.")]
write("round765_drafts/literature_scope_audit.json",
      json.dumps(dict(sources=sources, Hadamard_existence_for_current_model_claimed=False,
                      new_project_result="Complete coupled covariant gauge object and canonical/causal physical-algebra matching.",
                      source_text_copied=False), ensure_ascii=False, indent=2)+"\n")
files = ("research_note_765.md", "joint_covariant_gauge_complex.py",
         "joint_covariant_gauge_complex_results.json", "unified_physics_condition_ledger_765.md")
write("round765_drafts/final_review.txt",
      "Primary review only. Checked original action normalization and real derivative formal-adjoint signs. "
      "K includes covariant diffeomorphism on gauge and scalar fields, including curvature and Higgs zero-order blocks. "
      "D1=P-KKadj and D0=-KadjK satisfy exact intertwiners by on-shell Noether identity, with scalar wave principal symbol. "
      "Canonical matching uses the same Hessian boundary form and compact Cauchy support. "
      "Numerics do not simulate complete lower-order Hessian or propagators. "
      "Positive regular physical state, singular parametrix and positive physical Hadamard state are distinguished. "
      "No all-sector quantum Ward, original record realization or continuum Q-to-E theorem is claimed.\n"+
      "\n".join(n+" "+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+"\n")

mapping = {"764": "765", "763": "764", "3666": "3679", "3679": "3691",
           "3476": "3479", "3479": "3482", "1604": "1607",
           "joint_linear_physical_phase": "joint_covariant_gauge_complex"}
pattern = re.compile("|".join(map(re.escape, sorted(mapping, key=len, reverse=True))))


def convert(s):
    return pattern.sub(lambda m: mapping[m[0]], s)


verify = convert((HERE/"verify_round764.py").read_text("utf8"))
a, b = verify.index("    extra = ("), verify.index("    new = ")
verify = verify[:a]+"""    extra = ('unified_physics_condition_ledger_765.md',
             'round765_drafts/research_note_765_draft.md',
             'round765_drafts/final_review.txt',
             'round765_drafts/literature_scope_audit.json',
             'round765_drafts/scope_and_dedup_review.md',
             'round766_drafts/STATUS.md',
             'round765_drafts/joint_covariant_gauge_complex_before_check_repair.py',
             'round764_drafts/round755_working_report_observed_whitespace_change.md',
             'round764_drafts/historical_whitespace_repair.json')
"""+verify[b:]
verify = verify.replace("checks['display_formulas'] == 15", "checks['display_formulas'] == 14")
a, b = verify.index("                full_classical_linear_reduction_proven="), verify.index("                primary_code_and_note_review_completed=")
verify = verify[:a]+"""                complete_linear_covariant_gauge_complex_constructed=True,
                canonical_and_causal_physical_algebras_matched=True,
                spacelike_linear_observables_commute=True,
                numerical_checks_only_principal_symbols_and_specified_gauge_blocks=True,
                complete_lower_order_Hessian_or_propagator_simulated=False,
                Hadamard_bosonic_state_certified=False,
                all_sector_quantum_Ward_completed=False,
                original_Q_to_E_process_equivalence_proven=False,
"""+verify[b:]
write("verify_round765.py", verify)

publish = convert((HERE/"publish_round764.py").read_text("utf8"))
a, b = publish.index("summary = "), publish.index("planned = ")
publish = publish[:a]+"""summary = '**第765轮完成（共同因果量子对象）：** [完整耦合线性规范复形]({p}research_note_765.md)把原几何、规范与五标量Hessian接成共同双曲规范复形；因果物理代数与764初片空间辛同构，正则态及合法线性微因果性接通。物理Hadamard、全来源与原记录映射仍开放。三组、十四式通过，最新765／3482，1607份编号科学文件、3691份保护证据。[核验]({p}research_round_765_checks.json)、[条件账]({p}unified_physics_condition_ledger_765.md)。'
order = '**当前执行顺序（765后，优先于下方历史安排）：** 接[766规范相容的物理短距离准备]({p}round766_drafts/STATUS.md)，核同一高频物理结构、精确约束及正性；不继续主部/秩扫描。[范围审计]({p}round765_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及目标保持。'
"""+publish[b:]
publish = publish.replace("第765轮完成（领先物理量子分支）", "第765轮完成（共同因果量子对象）")
publish = publish.replace("## 410.", "## 411.").replace("## 315.", "## 316.")
publish = publish.replace("共同线性物理相空间与正态", "完整耦合线性规范复形与因果物理代数")
publish = publish.replace("旧空间合同保持，共同量子相空间与同阶来源", "旧空间合同保持，共同因果代数与物理短距离准备")
publish = publish.replace("共同线性物理相空间与内禀量子态", "完整耦合线性规范复形")
publish = publish.replace("next_round=765", "next_round=766")
write("publish_round765.py", publish)

post = convert((HERE/"postcheck_round764.py").read_text("utf8"))
post = post.replace("range(584, 765)", "range(584, 766)")
post = post.replace("第765轮完成（领先物理量子分支）", "第765轮完成（共同因果量子对象）")
write("postcheck_round765.py", post)
print("Prepared 765 reports, ledger and verification/publication scripts.")
