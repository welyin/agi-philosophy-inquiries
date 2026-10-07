"""Register one material/Maxwell interface and stop that technical branch."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
    STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
    "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 960：原生材料的电磁色散接口与有限有效域"
    old="## 959：Gauss父表示与共同恢复表"
    assert heading not in s and s.count(old)==1
    block=(heading+"\n\n"+f"[960报告]({pre}research_note_960.md)把实际电子材料的领先记录相位"
        +"接到已知Maxwell色散响应及同一距离来源；给非延迟能量／力界与高频加权尾界。"
        +"旧四点几何直接换点偶极会有约17.20%领先能量差，原实时记录界不自动继承。"
        +f"[结果]({pre}960/material_dispersion_bridge_results.json) · [核验]({pre}960/research_round_960_checks.json)。"
        +"正式960／累计3745，整体目标未完成。已有采用价值，停止色散力和器件优化；"
        +f"接[961]({pre}961/drafts/STATUS.md)按完整共同恢复选择同一有效域。\n\n")
    s=s.replace(old,block+old,1)
    s=s.replace("001—959轮共959份","001—960轮共960份")
    s=s.replace("959静电父接口已核，接[960]("+pre+"960/drafts/STATUS.md)核完整共同父描述的有效域",
        "960领先材料／Maxwell接口已核，接[961]("+pre+"961/drafts/STATUS.md)选择完整共同父描述的有效域")
    if p==paths[4]:
        s=s.replace("当前正式959／累计3744，959已结项","当前正式960／累计3745，960已结项")
    if p==paths[5]:
        s=s.replace("# 231—959轮阶段成果总览","# 231—960轮阶段成果总览",1)
        s=s.replace("231—959的729份","231—960的730份",1)
    if p==paths[-1]:
        s=s.replace("六条共同协议：全局缺口对应与检验优先级（截至959）",
                    "六条共同协议：全局缺口对应与检验优先级（截至960）",1)
        rows={
            "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—960|959静电Gauss精确字典；960实际材料谱的领先静态Maxwell色散匹配|3+1 Maxwell、电偶极及群为物理输入；不是完整规范物质或实时场论认证|",
            "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—960|959静电父字典精确保源；960有能量、距离力的非延迟误差和加权高频尾界|尾界需响应包络和低频匹配；不自动运输958实时界，不要求全UV|",
            "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；434—459、929、939、951、956—960|960同一电荷跃迁谱给958领先条件相位及成熟电磁响应，支持原生材料接口|制备、寿命及六协议仍按共同合同；已有方向价值，停止色散及器件优化|",
            "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—960|959不双计静电能；960从同一延迟能量微分保留传播核距离源|只核领先定态力及静态平衡，不等于完整应力或动态Einstein反馈|",
            "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—960|960旧R/a=3几何直接换点偶极令领先能量差约17.20%，限定实际可继承范围|内部匹配对照不是新宇宙预测，不因此要求修完所有几何或高阶问题|"}
        for key,row in rows.items():
            s,n=re.subn(r"^\|"+key+r" [^\n]*$",lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r"^4\. 957[^\n]*$",lambda m:
            "4. 957假说v0.2保持；959静电父表示和960领先Maxwell材料匹配各按原范围采用。"
            +"[共同恢复表v0.6](../../959/drafts/common_recovery_v0_6.md)保全部部门；"
            +"接[961](../../961/drafts/STATUS.md)按完整有效域判断共同采用，停止色散和器件优化。",
            s,count=1,flags=re.M);assert n==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Updated eight living navigation files for round 960.")
