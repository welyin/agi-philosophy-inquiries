"""Publish the reviewed quantum-reference interface and preserve navigation history."""
from pathlib import Path
import hashlib
import json
import os
import verify_quantum_reference_readout_review as evidence
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
checks = evidence.verify()
assert checks["total_protected_including_this_review"] == 1028
paths = [ROOT/"README.md", RESEARCH/"README.md", RESEARCH/"research_direction.md",
         RESEARCH/"RESEARCH_STATE.md", HERE/"README.md",
         HERE/"spatial_premise_closure_audit.md", HERE/"three_dimensional_four_conditions_review.md"]
raws = {p: p.read_bytes() for p in paths}
folder = HERE/"navigation_before_quantum_reference_readout_20260930"
folder.mkdir(exist_ok=False)
summary = (
    "**参考钟控制量子仪器接口（未编号）：** [审查]({p}quantum_reference_readout_review.md)"
    "把已给FLRW参考钟零模与被测qubit、探针保留为同一联合量子演化；固定有限窗口内，"
    "完整仪器对任意未知输入／旧参考的半diamond误差有O(V⁻¹ᐟ²)上界，并保留钟色散、回冲和背景功账。"
    "增加精度需要增加参考资源；本项读取系统Z，尚未记录空间参考值，也未闭合Einstein回馈。"
    "6项诊断、12式及独立终审通过；[结果]({p}quantum_reference_readout_results.json)、"
    "[核验]({p}quantum_reference_readout_checks.json)。编号仍520／2552，871份编号科学文件；保护证据1028份。")
next_steps = """

### 当前接续：把空间参考值变成实际记录

参考控制仪器已完成核验。它解决固定背景、正常钟零模和有限窗口下的真实联合测量近似；不会自动产生空间参考数值，也未给出新的自洽Einstein解。旧材料保持，最新冻结入口为verify_quantum_reference_readout_review.py。

**下一研究：** 在明确的局域参考场—探针模型中，计算探针实际记录与参考场数值的关系，核验不同位置能否在有限资源下区分，以及相同参考下多个记录是否可形成相容坐标。优先检查成熟局部测量接口及可精确计算的联合Gaussian模型；给定几何、边界、准备、读取和校准均须单列。需要实际传播／读取与反作用共同成立，不能只重写一般Kraus公式。

待核候选是有限区域内具有正联合二次型的参考场和探针场，保留局域耦合、探针源、全协方差及能量账。它可作为逆向物理输入的条件接口；原动态图到该有效模型的映射、连续极限和三维选择另待证明。若仅验证给定三维格点，不算推出三维。

**连接限制：** 将量子钟f(Q)乘到场—探针双线性耦合后，一般不再是精确Gaussian模型；场算符无界，不能直接继承本项对全部未知输入的统一误差。后续接合须给能量／矩约束、截断泄漏界或其它适用的受控近似。几何回馈须包含完整协变应力；源修正小不自动证明度规误差小。

目标继续采用用户确认的双向研究与逻辑统一标准，手动启动，不设置定时任务。宏观三维、GR生成及完整统一仍未完成。
"""
planned, manifest = {}, {}
for p, raw in raws.items():
    key = ("root_" if p.parent == ROOT else "research_" if p.parent == RESEARCH else "archive_") + p.name
    with (folder/key).open("xb") as f:
        f.write(raw)
    manifest[str(p.relative_to(ROOT))] = dict(snapshot=key, sha256=hashlib.sha256(raw).hexdigest())
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
    nl = "\r\n" if b"\r\n" in raw else "\n"
    body = raw.decode(enc).replace("\r\n", "\n")
    assert "参考钟控制量子仪器接口（未编号）" not in body
    prefix = "research_cognition_physics/archive_231_/" if p.parent == ROOT else "archive_231_/" if p.parent == RESEARCH else ""
    block = summary.format(p=prefix)
    if p in paths[:5]:
        title, tail = body.split("\n\n", 1)
        body = title + "\n\n" + block + "\n\n" + tail
    elif p.name == "spatial_premise_closure_audit.md":
        body += "\n\n## 162. 参考钟控制仪器与空间数值读取的区别（未编号）\n\n" + block + next_steps
    else:
        body += "\n\n## 67. 完整量子仪器误差已有界，空间坐标读取继续开放\n\n" + block + "\n\n下一步核实际参考场读数及跨记录相容性。qubit Z仪器、钟参考资源和三维空间不能相互替代。\n"
    if p.name == "RESEARCH_STATE.md":
        body += next_steps
    if p.name == "research_direction.md":
        body += "\n\n**下一数学接口更新：** 参考控制量子仪器已核验，后继直接检验空间参考值的实际记录及校准坐标；同时保留参考反作用。量子钟控制与精确Gaussian场模型不能未经证明拼接，全部输入保证也不能直接推广到无界场算符。详见RESEARCH_STATE.md最新接续。\n"
    planned[p] = body.replace("\n", nl).encode(enc)
with (folder/"manifest.json").open("x", encoding="utf8") as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)
count = 0
for p, body in planned.items():
    for link in core.link_parser()(body.decode("utf-8-sig")):
        assert (p.parent/link).resolve().exists(), (p, link)
        count += 1
for p, raw in raws.items():
    assert p.read_bytes() == raw, ("concurrent change", p)
for p, body in planned.items():
    temp = p.with_name(p.name+".quantum-reference.tmp")
    with temp.open("xb") as f:
        f.write(body)
    os.replace(temp, p)
assert evidence.verify() == checks
report = dict(date="2026-09-30", scientific_base_through_round=520,
    unchanged_numbered_scientific_tests=2552, numbered_scientific_files=871,
    protected_evidence=1028, navigation_files=7, navigation_links=count,
    broken_links=0, previous_navigation_snapshots_preserved=True,
    independent_final_review_completed=True, manual_goal_start_preserved=True,
    scientific_evidence_unchanged_after_navigation=True,
    stage_complete=False, all_checks_passed=True)
with (HERE/"quantum_reference_readout_integration_checks.json").open("x", encoding="utf8", newline="\n") as f:
    f.write(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
print(json.dumps(report, ensure_ascii=False))
