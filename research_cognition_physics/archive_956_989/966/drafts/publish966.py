"""Publish 966; preserve every prior scientific file and frozen entry."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"965/verify_round965.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"material_field_window_results.json")',1)[0]
prefix=prefix.replace("for n in range(776,965)","for n in range(776,966)")
prefix=prefix.replace("Delivery checks for 965","Delivery checks for 966")
prefix=prefix.replace("research_round_965_checks.json","research_round_966_checks.json")
checks=r'''    result=read(HERE/"comparison_loop_results.json")
    assert result["round"]==966 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert result["physical_sector_dimension"]==64
    assert result["low_organization_dimension"]==16
    assert sum(r["survives"] for r in result["exact_third_order_surviving_words"])==6
    assert result["tree_third_order_zero"]
    assert result["fixed_joint_dictionary_factorizes_H_exactly"]
    assert result["same_dictionary_makes_original_loop_effect_record_dependent"]
    assert not result["content_dependent_spectrum_generated"]
    from fractions import Fraction
    for row in result["fixed_menu_rows"]:
        e=Fraction(str(row["epsilon"]))
        eta=Fraction(63,8)*e+Fraction(117,16)*e*e+Fraction(3,2)*e**3
        assert eta==Fraction(row["all_input_uniform_error_bound"])
        assert Fraction(row["certified_loop_contrast_lower"])==Fraction(4,5)-4*eta
        assert Fraction(row["certified_mean_replacement_discrepancy_lower"])==Fraction(4,5)-2*eta
        assert row["actual_loop_contrast"]>=float(Fraction(4,5)-4*eta)>0
        assert row["direct_uniform_operator_error"]<float(eta)
        assert row["mixed_record_negative_flux_probability"]>float(Fraction(4,5)-2*eta)
        for entry in row["energy_parts"]:
            assert abs(sum(entry["changes"]))<1e-12
        assert row["full_energy_current_residual"]<1e-12
    assert result["fixed_menu_rows"][1]["certified_loop_contrast_lower_float"]>=.641768
    assert not result["scope"]["SM_or_GR_derived"]
    assert not result["scope"]["actual_spacetime_propagation_proved"]
    # Re-execute the complete physical-sector reconstruction and algebra.
    import importlib.util
    spec=importlib.util.spec_from_file_location("core966",HERE/"comparison_loop.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.compute(),result)
    note=STAGE/"research_note_966.md"
    prose=note.read_text("utf-8")
    assert prose.count("$$")==18
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,10)]
    for term in ("整体目标未完成","降为说明性例","固定联合变换","额外建模输入",
                 "不是对FUCP","不能保持原准备和全部效果不变而删除记录因子"):
        assert term in prose,term
    newfiles=[note,HERE/"comparison_loop.py",HERE/"comparison_loop_results.json",
        Path(__file__),HERE/"drafts/mechanism_priority_entry.md",
        HERE/"drafts/mechanism_map_v0_2.md",HERE/"drafts/publish966.py",
        STAGE/"967/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
tail=tail.replace("截至965","截至966").replace("n<=965","n<=966")
tail=tail.replace("list(range(1,966))","list(range(1,967))")
tail=tail.replace("001—965轮共965份","001—966轮共966份")
tail=tail.replace("231—965的735份","231—966的736份")
tail=tail.replace("965：原生材料的共同电磁读口","966：局部比较的闭路生成与表示筛选")
tail=tail.replace('STAGE/"964/research_round_964_checks.json"','STAGE/"965/research_round_965_checks.json"')
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=966,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=966,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_965=prev["cumulative_numbered_test_groups_from_964"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        exact_comparison_induced_record_loop=True,
        finite_time_all_input_bound=True,
        exact_joint_dictionary_factorization=True,
        candidate_deprioritized_after_representation_screen=True,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"966/verify_round966.py"
assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 966：局部比较的闭路生成与表示筛选"
    oldheading="## 965：原生材料的共同电磁读口"
    assert heading not in s and s.count(oldheading)==1
    block=(heading+"\n\n"
       +f"[966报告]({pre}research_note_966.md)提出并筛选一份内部比较规则："
       +"局部比较生成记录—闭路项，固定窗口概率差解析下界.641768；"
       +"但完整联合换表示能从H消去记录符号，原实际效果仍含记录。"
       +"候选保留为说明性例，停止扩建，不把可读响应当作内容选择内在动力学的证明。"
       +f"[结果]({pre}966/comparison_loop_results.json) · "
       +f"[核验]({pre}966/research_round_966_checks.json)。正式966／累计3751，整体目标未完成。\n\n"
       +f"[机制图v0.2]({pre}966/drafts/mechanism_map_v0_2.md)区分有限响应、表示变换与内容依赖谱源；"
       +"三体比较、群、能隙及初态仍为输入，未接成958／965同一材料。"
       +f"接[967]({pre}967/drafts/STATUS.md)回查435等已有内容改变作用幅度的机制，"
       +"先判断共同接合价值。现行目标已包含机制优先，应用目标及任务设置保持。\n\n")
    s=s.replace(oldheading,block+oldheading,1)
    s=s.replace("001—965轮共965份","001—966轮共966份")
    if p==paths[4]:
        s=s.replace("当前正式965／累计3750，965已结项","当前正式966／累计3751，966已结项")
    if p==paths[5]:
        s=s.replace("# 231—965轮阶段成果总览","# 231—966轮阶段成果总览",1)
        s=s.replace("231—965的735份","231—966的736份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至965）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至966）",1)
        rows={
         "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—966|965共同读口保留；966同一记录标签与闭路效果有有限窗口响应|966换表示须共同变换效果；原标签保持不是未知裸态隔离|",
         "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、568、930—934、959—966|966输入Z2下从比较生成闭路项，系数由J与能隙决定|不是568的SU(2)材料诱导；群、Gauss扇区、三体比较仍输入|",
         "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—966|965无限占据有界效果界保留；966精确有理残差给全低部门有限时间界|两个不同H，不拼接为同一实现；均不要求全UV完成|",
         "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—966|原生材料结果保留；966具体比较保三个逻辑记录并给关系响应|候选经表示筛选降为说明性例；未与原生材料合并，停止装置扩建|",
         "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—966|966同一H给组织、比较及关系能源回流和一致参数来源|有限内部守恒不产生共同度规；符号记录不构成不可消去的谱源|",
         "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—966|966闭路差下界.641768；错误均值替换的有限概率误差下界.720884|精确换表示保持预测却削弱生成解释；967优先回查已有内容幅度机制|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
        before="4. 957假说v0.2保持；[961整体机制图](../../research_note_961.md)为现行解释入口，与[959物理恢复表](../../959/drafts/common_recovery_v0_6.md)分开使用。按草图→冲突→最小检验→修订→系统验证推进；接[962](../../962/drafts/STATUS.md)，暂缓局部优化。"
        after="4. 957假说v0.2保持；[961整体机制图](../../research_note_961.md)及[966增量v0.2](../../966/drafts/mechanism_map_v0_2.md)为解释入口，与[965物理恢复表](../../965/drafts/common_recovery_v0_8.md)分开使用。966候选经表示筛选降为说明性例；接[967](../../967/drafts/STATUS.md)优先回查既有内容幅度机制与共同材料，停止局部扩建。"
        assert before in s
        s=s.replace(before,after,1)
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 966 preservation verifier.")

