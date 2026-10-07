"""Update living navigation after the worldline-interface round."""
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
    heading="## 953：有效物体的协变作用与共同能源交换"
    assert heading not in s
    block=(heading+"\n\n"+f"[953报告]({pre}research_note_953.md)将原h、B接入中性世界线作用，"
        "同一项给内部酉、标量源及四维能源交换；有限尺寸阻止直接继承947误差。"
        +f"[结果]({pre}953/worldline_material_bridge_results.json) · [核验]({pre}953/research_round_953_checks.json)。"
        +"正式953／累计3738。协变单极接口已核，完整共同模型未完成；停止世界线修补，"
        +f"接[954]({pre}954/drafts/STATUS.md)审门户是否恢复所声称的物理极点与响应。\n\n")
    s,n=re.subn(r"## 953工作期执行更正：先判共同接合，暂缓动态细化\n\n.*?\n\n",
        lambda m:block,s,count=1,flags=re.S)
    assert n==1,p
    # These are living next-action clauses, not frozen past reports.
    s=re.sub(r"按\[953(?:取舍)?复核\]\([^)]+priority_reaudit\.md\)先判共同接合，暂缓(?:原定)?动态细化",
        lambda m:f"953共同作用接口已核，接[954]({pre}954/drafts/STATUS.md)审门户物理恢复",s)
    s=s.replace("001—952轮共952份","001—953轮共953份")
    if p in paths[:4] or p==paths[-1]:
        # Current narrative retains previous 952 findings; update its status only.
        s=s.replace("正式952／累计3737，整体目标未完成","正式953／累计3738，整体目标未完成")
    if p==paths[4]:
        s=s.replace("当前正式952／累计3737，952已结项",
                    "当前正式953／累计3738，953已结项")
    if p==paths[5]:
        s=s.replace("# 231—952轮阶段成果总览","# 231—953轮阶段成果总览",1)
        s=s.replace("231—952的722份","231—953的723份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至952）",
                    "六条共同协议：全局缺口对应与检验优先级（截至953）",1)
        rows={
            "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—953|953中性单极材料与原规范不变门户共用作用；没有新增直接规范荷|非零门户改变物理传播，规范不变不等于完整SM预测恢复；群和场表仍输入|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、854、898、944—953|953明确协变单极与原Gaussian材料的有限尺寸字典；不能直接继承947界|保形状或长波匹配各有范围；完整共同误差未认证，不自动展开点源UV|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—953|953同一h、B交互给世界线标量源及四维应力交换，光滑外场单极阶Ward身份闭合|尚非耦合量子物体—场—动态几何过程；Einstein和3+1是输入|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M)
            assert n==1,key
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 953.")
