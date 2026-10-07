"""Publish finite coordinate records and the positive-alignment follow-up."""
from pathlib import Path
import hashlib
import json
import os
import verify_spatial_reference_interfaces as evidence
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
before=evidence.verify()
assert before["total_protected_including_this_review"]==1038
paths=[ROOT/"README.md",RESEARCH/"README.md",RESEARCH/"research_direction.md",
       RESEARCH/"RESEARCH_STATE.md",HERE/"README.md",
       HERE/"spatial_premise_closure_audit.md",HERE/"three_dimensional_four_conditions_review.md"]
raws={p:p.read_bytes() for p in paths}
folder=HERE/"navigation_before_spatial_reference_interfaces_20260930"
folder.mkdir(exist_ok=False)
summary=("**空间参考的实际记录与稳定性接续（未编号）：** "
    "[局域读取]({p}local_reference_probe_review.md)在给定有限格点和参考源上，以同一联合量子演化、"
    "正噪声仪器和有限符号记录，证明单次同时识别全部格点坐标的置信保证；大参考能量与末读反作用均计入。"
    "原裸耦合有扩大区域失稳限制。[共同参考正势]({p}reference_alignment_completion_review.md)"
    "给出明确规则修订：共同模自由、相对模稳定；可精确嵌入旧经典Einstein参考解，量子读取仍含相对模噪声。"
    "两个接口的准备不同，尚未拼成完整量子—引力模型，也未选择三维。共12组诊断、21式与交叉复核通过；"
    "[核验]({p}spatial_reference_interfaces_checks.json)。编号仍520／2552、871份编号科学文件；保护证据1038份。")
next_steps="""

### 当前接续：同一正势规则中的参考形成与坐标读取

已经闭合的有限接口：

1. 给定有限格点、参考梯度和Gaussian准备，真实局部参考—探针作用可产生坐标读数；用正噪声仪器与有限分箱，一次读取同时正确识别所有格点的概率至少1−α。证书覆盖同一联合记录，不使用格点独立噪声或重复读取的独立性假设。示例参考能量很大，不是廉价定位证明。
2. 原裸双线性势不能以固定非零耦合无限扩大。新增完成平方势使共同模无质量、相对模有正质量；全部有限Dirichlet区域的二次型保持正，但能隙随区域扩大趋零。
3. 正势模型的共同参考分支可精确嵌入原经典Einstein—标量解。它从一开始便在两套场中共享参考，不继承空白探针准备；量子相对模的涨落、读取加能和相关回冲仍不可删除。

**下一可检验问题：** 在同一正势局域规则下，从明确的独立探针准备开始，实际参考—探针演化能否同时保留坐标区分证书、有限记录和规模扩展？须按修改后的势重新计算双方源、初始边界及均值；不可只给原矩阵加稳定项却沿用旧强迫项。若采用预先共享参考分支，应明确其形成机制及准备关联，不能把它当作空白探针自然得到的结果。

保留同一标准：给定几何是逆向模型输入，格点数、场数量或内部qubit坐标不自动选择三维；经典Einstein解不自动成为量子协方差应力的解。后续同时审查有限区域到有效几何的适用条件、源反作用与维数选择；无需重算一般辛矩阵恒等式、裸势失稳、旧高斯尾界或只优化b。

最新冻结入口verify_spatial_reference_interfaces.py，继承1038份保护证据；编号520／2552保持。目标手动继续，不设置定时任务。空间及完整统一阶段未结项。
"""
manifest,planned={},{}
for p,raw in raws.items():
    key=("root_" if p.parent==ROOT else "research_" if p.parent==RESEARCH else "archive_")+p.name
    with (folder/key).open("xb") as f:f.write(raw)
    manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    enc="utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf8"
    nl="\r\n" if b"\r\n" in raw else "\n"
    body=raw.decode(enc).replace("\r\n","\n")
    assert "空间参考的实际记录与稳定性接续（未编号）" not in body
    prefix="research_cognition_physics/archive_231_/" if p.parent==ROOT else "archive_231_/" if p.parent==RESEARCH else ""
    block=summary.format(p=prefix)
    if p in paths[:5]:
        title,tail=body.split("\n\n",1)
        body=title+"\n\n"+block+"\n\n"+tail
    elif p.name=="spatial_premise_closure_audit.md":
        body+="\n\n## 163. 从局域场演化到有限坐标记录（未编号）\n\n"+block
        body+="\n\n## 164. 共同参考正势与同一准备的待接合条件（未编号）\n"+next_steps
    else:
        body+="\n\n## 68. 实际坐标记录可以成立，所用场数仍是输入\n\n"+block
        body+="\n\n## 69. 正势稳定性没有选择空间维数\n\n完成平方后的共同模可用于所有声明维数。条件性参考读取与经典Einstein嵌入已分别核验，三维选择仍须在同一框架中额外论证。\n"
    if p.name=="RESEARCH_STATE.md":body+=next_steps
    if p.name=="research_direction.md":
        body+="\n\n**下一接口更新：** 同一正势局域规则中的参考形成、空白探针读取和规模扩展。新增势需重算全部源与准备；既有共享参考分支不能替代形成过程。详见RESEARCH_STATE.md最新接续。维数选择及底层到有效场映射继续开放。\n"
    planned[p]=body.replace("\n",nl).encode(enc)
with (folder/"manifest.json").open("x",encoding="utf8") as f:
    json.dump(manifest,f,ensure_ascii=False,indent=2)
count=0
for p,body in planned.items():
    for link in core.link_parser()(body.decode("utf-8-sig")):
        assert (p.parent/link).resolve().exists(),(p,link)
        count+=1
for p,raw in raws.items():assert p.read_bytes()==raw,("concurrent change",p)
for p,body in planned.items():
    temp=p.with_name(p.name+".spatial-reference.tmp")
    with temp.open("xb") as f:f.write(body)
    os.replace(temp,p)
assert evidence.verify()==before
report=dict(date="2026-09-30",scientific_base_through_round=520,
    unchanged_numbered_scientific_tests=2552,numbered_scientific_files=871,
    protected_evidence=1038,navigation_files=7,navigation_links=count,broken_links=0,
    previous_navigation_snapshots_preserved=True,independent_cross_review_completed=True,
    scientific_evidence_unchanged_after_navigation=True,manual_goal_start_preserved=True,
    stage_complete=False,all_checks_passed=True)
with (HERE/"spatial_reference_interfaces_integration_checks.json").open("x",encoding="utf8",newline="\n") as f:
    f.write(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(report,ensure_ascii=False))
