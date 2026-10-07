"""Publish 967 and verify preservation, exact energy intervals and shared inputs."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"966/verify_round966.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"comparison_loop_results.json")',1)[0]
prefix=prefix.replace("for n in range(776,966)","for n in range(776,967)")
prefix=prefix.replace("Delivery checks for 966","Delivery checks for 967")
prefix=prefix.replace("research_round_966_checks.json","research_round_967_checks.json")
checks=r'''    result=read(HERE/"collective_material_results.json")
    assert result["round"]==967 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["parameters"]==dict(U=1.,v=.1,absolute_kappa=.2,
        active_dimension=216,full_with_old_spectator_spins=13824)
    from fractions import Fraction as F
    assert F(result["loop_sign_energy_difference_interval"]["lo"])>0
    bound=2*F(9,250)*F(347,1000)/F(219,500)
    assert result["fixed_menu_rows"][0]["connected_energy"]>0
    assert result["fixed_menu_rows"][1]["connected_energy"]<0
    for row in result["fixed_menu_rows"]:
        c=row["ground_interval"];z=row["connected_energy_interval"]
        assert c["lower_inertia"]==[0,27] and c["upper_inertia"]==[1,26]
        assert F(c["hi"])-F(c["lo"])==F(5,10**12)
        assert F(z["lo"])*F(z["hi"])>0
        assert F(row["all_time_uniform_input_error_bound"])==bound
        assert row["direct_uniform_error"]<float(bound)
        assert max(row["intertwining_error"],row["isometry_error"])<2e-13
        assert row["pair_only_prediction_probability"]==1.
        assert row["actual_receiver_probability"]<2e-5
        assert row["certified_prediction_discrepancy_lower"]>.885
        assert row["full_trace"]==216.
        assert row["energy_current_balance"]<1e-12
    assert not result["scope"]["independent_three_body_coupling_introduced"]
    assert not result["scope"]["autonomous_geometry_or_binding_proved"]
    # Recompute actual CAR-to-singlet mapping, rational LDL signs, all 64
    # dressed columns, the physical receiving effect and shared sources.
    import importlib.util
    spec=importlib.util.spec_from_file_location("core967",HERE/"collective_material.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_967.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","第三份材料","精确双体项","接触系数仍为外部固定输入",
                 "同一H的不同记录部门","不反证量子局部层析"):
        assert term in prose,term
    newfiles=[note,HERE/"collective_material.py",HERE/"collective_material_results.json",
        Path(__file__),HERE/"drafts/collective_material_decision.md",
        HERE/"drafts/collective_mechanism_increment.md",HERE/"drafts/publish967.py",
        STAGE/"968/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
tail=tail.replace("截至966","截至967").replace("n<=966","n<=967")
tail=tail.replace("list(range(1,967))","list(range(1,968))")
tail=tail.replace("001—966轮共966份","001—967轮共967份")
tail=tail.replace("231—966的736份","231—967的737份")
tail=tail.replace("966：局部比较的闭路生成与表示筛选","967：原生材料的集体作用与共同来源")
tail=tail.replace('STAGE/"965/research_round_965_checks.json"','STAGE/"966/research_round_966_checks.json"')
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=967,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=967,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_966=prev["cumulative_numbered_test_groups_from_965"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_rational_ground_intervals=True,
        same_native_material_and_density_rule=True,
        no_independent_three_body_parameter=True,
        all_low_inputs_and_passive_references_bound=True,
        autonomous_geometry_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"967/verify_round967.py";assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 967：原生材料的集体作用与共同来源"
    oldheading="## 966：局部比较的闭路生成与表示筛选"
    assert heading not in s and s.count(oldheading)==1
    block=(heading+"\n\n"
       +f"[967报告]({pre}research_note_967.md)把435的闭合虚过程接到958原材料："
       +"同一密度规则生成非加性集体能量，不增独立三体系数；"
       +"孤立双体谱相同的两个符号菜单有严格不同的共同谱。"
       +"只加精确双体低有效项会导致>.885的有限接收概率偏差，同一电荷来源也需集体修正。"
       +f"[结果]({pre}967/collective_material_results.json) · "
       +f"[核验]({pre}967/research_round_967_checks.json)。正式967／累计3752，整体目标未完成。\n\n"
       +f"[机制增量]({pre}967/drafts/collective_mechanism_increment.md)保留组织、虚电荷、集体记录及来源的共同链；"
       +"三块图与接触仍输入，未完成自主几何或965场的共同接合。"
       +f"停止三体和装置优化，接[968]({pre}968/drafts/STATUS.md)先判定内部关系更新变量与实际传播的机制。"
       +"应用目标及任务设置保持。\n\n")
    s=s.replace(oldheading,block+oldheading,1)
    s=s.replace("001—966轮共966份","001—967轮共967份")
    if p==paths[4]:
        s=s.replace("当前正式966／累计3751，966已结项","当前正式967／累计3752，967已结项")
    if p==paths[5]:
        s=s.replace("# 231—966轮阶段成果总览","# 231—967轮阶段成果总览",1)
        s=s.replace("231—966的736份","231—967的737份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至966）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至967）",1)
        rows={
         "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—967|967沿同一旧关系码给三方接收预测；机械相加双体项的偏差>.885|初态与效果可访问性输入；原双体时间表不自动保持|",
         "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、568、930—934、959—967|967电荷取向重命名保整环乘积，闭路虚过程改变实际共同谱|不是SM群或磁plaquette推导；二轨道与固定接触输入保留|",
         "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—967|967精确衣着给全64维低输入的全时间界；有理LDL认证非加性能量|同材料家族三块扩展；没有把965场或964背景自动并入|",
         "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—967|967同一原生电荷规则固定集体记录项，复用435机制而不另设三体门|接触、制备、旁观材料及寿命仍输入；停止三体工程|",
         "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—967|967集体能量与Q_iQ_j来源来自同一H，不能沿用孤立双边源相加|κ仍为固定菜单；参数来源尚非关系移动、共同度规或Einstein动力学|",
         "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—967|967同双体谱的环符号菜单有严格共同谱差；原记录检验排除机械相加|有限模型而非经验确认；完整统一与自主关系更新仍未完成|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        before="4. 957假说v0.2保持；[961整体机制图](../../research_note_961.md)及[966增量v0.2](../../966/drafts/mechanism_map_v0_2.md)为解释入口，与[965物理恢复表](../../965/drafts/common_recovery_v0_8.md)分开使用。966候选经表示筛选降为说明性例；接[967](../../967/drafts/STATUS.md)优先回查既有内容幅度机制与共同材料，停止局部扩建。"
        after="4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[966表示筛选](../../966/drafts/mechanism_map_v0_2.md)及[967共同机制增量](../../967/drafts/collective_mechanism_increment.md)与[965物理恢复表](../../965/drafts/common_recovery_v0_8.md)分开使用。967已接通同材料的集体比较与来源；接[968](../../968/drafts/STATUS.md)先审实际内部关系变量，不扩建三体或器件。"
        assert before in s;s=s.replace(before,after,1)
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 967 preservation verifier.")

