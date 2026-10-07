"""Publish round 964; historical scientific artifacts remain frozen."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths}
updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 964：同一材料的膨胀、记录与热账"
    old="## 963：共同钟标定的有限证书与EPS接口"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"
      +f"[964报告]({pre}research_note_964.md)把958原材料接入明示的均匀半经典作用；"
      +"同一尺度增大一倍时光／材料频比减半，原记录概率差保持，平均来源与压力功闭合。"
      +"热模熵不变、材料记录关联可回落，膨胀不自动提供时间箭头。"
      +f"[结果]({pre}964/material_cosmology_results.json) · [核验]({pre}964/research_round_964_checks.json)。"
      +"正式964／累计3749，整体目标未完成。\n\n"
      +f"[共同恢复v0.7]({pre}964/drafts/common_recovery_v0_7.md)明确均匀／单极／期望值输入，"
      +"没有把自由光子、实际光学交互及量子引力涨落混为一项完成。"
      +f"停止背景和器件优化；接[965]({pre}965/drafts/STATUS.md)审原生物质与真实传播交互的共同窗口。"
      +"应用目标及任务设置保持。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—963轮共963份","001—964轮共964份")
    if p==paths[4]:
        s=s.replace("当前正式963／累计3748，963已结项","当前正式964／累计3749，964已结项")
    if p==paths[5]:
        s=s.replace("# 231—963轮阶段成果总览","# 231—964轮阶段成果总览",1)
        s=s.replace("231—963的733份","231—964的734份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至963）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至964）",1)
        rows={
        "C09":"|C09 钟尺、参考态与可访问性|1、3、6；522—530、941、962—964|964原材料谱与共形光子谱在同一平均背景有不同缩放，频比可变|单极质量和几何仍输入；未实现发射／吸收／比较仪器，963合同仍条件性|",
        "C19":"|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957、964|964同一有限光子菜单的Gibbs态温度随a降低而熵保持|不是完整黑体或与材料的热平衡；初态和省略模式匹配仍输入|",
        "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—964|964同一材料完整内能、光子压力和均匀Einstein约束共用变分与压力功|期望值／均匀／单极层；未认证真实非均匀源、散射及量子几何反馈|",
        "C23":"|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964|964膨胀持续时实际材料记录关联可回落，光子和材料整体熵不增|膨胀不能替代非平衡准备、可用容量与记录稳定机制；宏观箭头未完成|",
        "C24":"|C24 宇宙初态、Λ与暗部门|1、3、6；303—304、325、331—340、961、964|964给原材料、自由红移、局部内谱及平均能源的同一膨胀解|Einstein／3+1／稳定物体／初态为输入；真实传播交互、暗部门及生成机制仍未完成|",
        "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—964|964同一材料菜单区分实际频比变化与全速率重标，保持原读数并检查熵账|不是宇宙观测或完整E_nat；965优先核真实物质—传播场共同交互|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living entries through 964; scientific history preserved.")

