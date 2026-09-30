"""Publish the reviewed material-reference interface; preserve all prior evidence."""
from pathlib import Path
import hashlib
import json
import os
import verify_material_reference_geometry_review as evidence
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checks=evidence.verify()
assert checks["total_protected_including_this_review"]==1021
paths=[ROOT/"README.md",RESEARCH/"README.md",RESEARCH/"research_direction.md",
       RESEARCH/"RESEARCH_STATE.md",HERE/"README.md",
       HERE/"spatial_premise_closure_audit.md",HERE/"three_dimensional_four_conditions_review.md"]
raws={p:p.read_bytes() for p in paths}
folder=HERE/"navigation_before_material_reference_20260930"
folder.mkdir(exist_ok=False)
manual="**运行方式（用户2026-09-30确认）：** 当前目标已由用户更新并手动恢复；后续由用户手动启动或继续，不设置定时研究。"
old_schedule="既有自动研究已同步这一阶段目标，保留每天07:00运行及只在有意义结果、完成、失败或需要决定时通知的设置。"
new_schedule="用户随后明确改为手动继续目标；原定时研究已删除，目标正文已由用户更新并手动恢复。"
summary=("**物理参考场与关系几何接口（未编号）：** [审查]({p}material_reference_geometry_review.md)"
 "把成熟物理参考场接到同一Einstein—标量解：参考读数可定义坐标，参考应力同时进入几何，"
 "并得到同一零方向与有限读数精度的条件性资源账。所选钟在无限固有时间内趋于有限读数，"
 "同一图的精度会恶化；这是参考方案限制。构造适用于所有d≥2，未选择三维，"
 "未把经典场、实际量子读取或原图连续极限混同。"
 "6项诊断与11式通过；[结果]({p}material_reference_geometry_results.json)、"
 "[核验]({p}material_reference_geometry_checks.json)。编号仍520／2552，871份编号科学文件；保护证据1021份。")
next_steps="""

### 当前接续：从物理参考图到实际量子读取

物理参考场的相容构造已完成未编号核验。它替换了“必须先把所有钟尺都从原微模型无额外输入造出”这一唯一推进顺序，但尚未建立原图模型到有效场的映射。旧轮冻结与全部计数不变。

已核读Fewster–Verch《Quantum fields and local measurements》v3的§1、§3和§4（https://arxiv.org/html/1810.06512v3）：在固定、可弯曲、整体双曲背景上，给定局部耦合、探针准备和末端读取，能由散射映射定义系统的量子仪器；因果分解是组合一致性的条件。其范围不包括完整动态量子引力，也没有自动产生探针控制。

**下一候选：** 用物理参考读数指定系统—探针交互的区域，并核实际量子仪器、参考反作用与几何账能否在同一有限精度模型中接续。一个可检验的有效作用量选择是质量正的两场Φ、Ψ及势V=(mΦ²Φ²+mΨ²Ψ²)/2+λ f(X)ΦΨ；取光滑有界f及|λ f|<mΦ mΨ，可使该二次物质势正定。f的形状和作用量仍为模型选择。

把原来外部指定的开关c(x)改为λ f(X)，会使参考场的方程出现λ(∂f/∂X^A)ΦΨ源；不能继续把X当作无扰动的自由参考解。相关局部Noether能量账直接复用325，不为重新求散度增加轮次。Fewster–Verch的固定开关线性模型也不能直接替这个非线性联合模型担保：后继须说明采用的量子化/有效近似、读取精度、反作用余项、实际准备及旧参考关联的范围。

如果只是在固定背景重算一般Kraus公式或复制记录，则不编号。新编号应关闭一个明确的连接，例如给出同一过程下可验证的读数误差与参考扰动控制。三维仍为近期主线之一；不同维数的参考模型怎样受到同一组合和测量要求约束，保持开放。

最新冻结入口为verify_material_reference_geometry_review.py；1021份保护证据、520／2552保持。未宣告空间、GR生成或统一阶段完成。
"""
planned={}
manifest={}
for p,raw in raws.items():
    key=("root_" if p.parent==ROOT else "research_" if p.parent==RESEARCH else "archive_")+p.name
    (folder/key).write_bytes(raw)
    manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    enc="utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
    nl="\r\n" if b"\r\n" in raw else "\n"
    body=raw.decode(enc).replace("\r\n","\n")
    assert "物理参考场与关系几何接口（未编号）" not in body
    body=body.replace(old_schedule,new_schedule)
    prefix="research_cognition_physics/archive_231_/" if p.parent==ROOT else "archive_231_/" if p.parent==RESEARCH else ""
    block=summary.format(p=prefix)
    if p in paths[:5]:
        parts=body.split("\n\n",1)
        body=parts[0]+"\n\n"+manual+"\n\n"+block+"\n\n"+parts[1]
    elif p.name=="spatial_premise_closure_audit.md":
        body+="\n\n## 161. 物理参考场与实际读取的双向接口（未编号）\n\n"+block+next_steps
    else:
        body+="\n\n## 66. 物理参考图可以成立，但参考数量不选择维数\n\n"+block+"\n\n"+manual+"\n\n同一模型的d≥2族说明：d+1个参考读数来自所选维数，不能反过来作为三维来源。下一步将参考可读性与测量反作用联合核算。\n"
    if p.name=="RESEARCH_STATE.md":
        body+=next_steps
    if p.name=="research_direction.md":
        body+="\n\n**下一数学接口：** 接物理参考场与局域量子测量，核同一过程的参考反作用、读数精度及组合一致性。成熟Fewster–Verch接口的固定背景/固定开关范围须保留；详见RESEARCH_STATE.md当前接续。维数选择与原图连续映射仍开放。\n"
    planned[p]=body.replace("\n",nl).encode(enc)
with (folder/"manifest.json").open("x",encoding="utf8") as f:
    json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in raws.items():
    assert p.read_bytes()==raw,("concurrent change",p)
for p,new in planned.items():
    temp=p.with_name(p.name+".material-reference.tmp")
    with temp.open("xb") as f:f.write(new)
    os.replace(temp,p)
count=0
for p in paths:
    text=p.read_text("utf-8-sig")
    assert old_schedule not in text
    assert "物理参考场与关系几何接口（未编号）" in text
    for link in core.link_parser()(text):
        assert (p.parent/link).resolve().exists(),(p,link)
        count+=1
after=evidence.verify()
assert after==checks
report=dict(date="2026-09-30",scientific_base_through_round=520,
            unchanged_numbered_scientific_tests=2552,numbered_scientific_files=871,
            protected_evidence=1021,navigation_files=7,navigation_links=count,
            broken_links=0,manual_goal_start_recorded=True,
            stale_daily_schedule_statement_removed=True,
            previous_navigation_snapshots_preserved=True,
            new_goal_is_active_by_user_request=True,
            independent_final_review_completed=True,
            next_quantum_measurement_interface_sources_checked=True,
            scientific_evidence_unchanged_after_navigation=True,
            stage_complete=False,all_checks_passed=True)
with (HERE/"material_reference_geometry_integration_checks.json").open("x",encoding="utf8",newline="\n") as f:
    f.write(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(report,ensure_ascii=False))
