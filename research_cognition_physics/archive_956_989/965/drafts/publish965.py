"""Publish 965 using the frozen 964 preservation machinery."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/"964/verify_round964.py").read_text("utf-8")
prefix=old.split('    result=read(HERE/"material_cosmology_results.json")',1)[0]
prefix=prefix.replace("for n in range(776,964)","for n in range(776,965)")
prefix=prefix.replace("Delivery checks for 964","Delivery checks for 965")
prefix=prefix.replace("research_round_964_checks.json","research_round_965_checks.json")
checks=r'''    result=read(HERE/"material_field_window_results.json")
    physical=read(HERE/"physical_mode_effect_results.json")
    assert result["round"]==965 and result["all_scientific_checks_passed"]
    assert physical["round"]==965 and physical["all_checks_passed"]
    for data in (result,physical):
        for rel,digest in data["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    previous=read(STAGE/"958/capacitive_material_write_results.json")
    assert result["T"]==previous["exact_interface"]["controlled_phase_time"]
    assert result["old_material_parameters"]==dict(U=1.,v=.1,kappa=.2)
    assert [r["cutoff"] for r in result["fixed_menu_rows"]]==[7,10]
    for row in result["fixed_menu_rows"]:
        eps=row["residual_bound"]["total"]
        assert eps<2e-6 and row["residual_bound"]["defining_matrix_operator_error"]<1e-12
        assert row["probabilities"]["receiver_contrast"]-4e-6>.998723
        assert row["algebra"]["low_charge_projection"]<1e-14
        assert row["algebra"]["projected_square_norm"]>.07
        assert row["algebra"]["charge_internal_commutator"]>.39
        for energy in row["sources"]["energy_rows"]:
            assert abs(sum(energy["after"])-sum(energy["before"]))<1e-11
        assert row["sources"]["current_identity_residual"]<1e-12
        assert not row["scope"]["local_propagation_or_full_Maxwell_recovered"]
        assert not row["scope"]["unbounded_energy_observable_error_transported"]
    for row in physical["rows"]:
        assert row["conservative_reported_lower"]>.049871
        assert row["maximum_compressed_effect_norm"]<=1+1e-12
    # Independently check that the near-unitary normalization defect is covered
    # by the explicit endpoint guard, on the WHOLE declared input isometry.
    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location("core965",HERE/"material_field_window.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    m,Q,Hm,S,W0,local=mod.material()
    n=11;g=.003;omega=.5
    a=np.diag(np.sqrt(np.arange(1,n)),1)
    H=np.kron(Hm,np.eye(n))+np.kron(np.eye(36),omega*np.diag(np.arange(n)))
    H+=g*np.kron(S,a+a.T)+g*g/omega*np.kron(S@S,np.eye(n))
    ev,V=np.linalg.eigh(H);Vin=np.kron(W0,np.eye(n)[:,:2])
    output=V@(mod.phases(ev,result["T"])[:,None]*(V.T@Vin))
    gram_defect=float(np.linalg.norm(output.conj().T@output-np.eye(32)))
    assert gram_defect<1e-10
    note=STAGE/"research_note_965.md"
    prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("单模读口不等于空间传播","联合算符","整体目标未完成",
                 "没有用有界读数的全态误差自动认证无界"):
        assert term in prose,term
    newfiles=[note,HERE/"material_field_window.py",HERE/"material_field_window_results.json",
        HERE/"physical_mode_effect.py",HERE/"physical_mode_effect_results.json",
        Path(__file__),HERE/"drafts/electromagnetic_decision.md",
        HERE/"drafts/common_recovery_v0_8.md",HERE/"drafts/publish965.py",
        STAGE/"966/drafts/STATUS.md"]
'''
tail=old.split('    nav=[STAGE.parent/n for n in',1)[1]
tail='    nav=[STAGE.parent/n for n in'+tail
tail=tail.replace("截至964","截至965").replace("n<=964","n<=965")
tail=tail.replace("list(range(1,965))","list(range(1,966))")
tail=tail.replace("001—964轮共964份","001—965轮共965份")
tail=tail.replace("231—964的734份","231—965的735份")
tail=tail.replace("964：同一材料的膨胀、记录与热账","965：原生材料的共同电磁读口")
tail=tail.replace('STAGE/"963/research_round_963_checks.json"','STAGE/"964/research_round_964_checks.json"')
start=tail.index("    out=dict(");end=tail.index("    if not writing:",start)
tail=tail[:start]+'''    out=dict(round=965,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=965,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_964=prev["cumulative_numbered_test_groups_from_963"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,uniform_input_gram_defect=gram_defect,
        old_material_and_receiver_reused=True,
        infinite_occupation_residual_bound_used=True,
        field_observable_includes_polarization=True,
        single_mode_not_local_propagation=True,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
verifier=STAGE/"965/verify_round965.py"
assert not verifier.exists()
verifier.write_text(prefix+checks+tail,encoding="utf-8")

paths=[RESEARCH/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
 STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                  "_shared/notes/unified_physics_condition_ledger_current.md")]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode("utf-8-sig").replace("\r\n","\n")
    pre="archive_764_/" if p.parent==RESEARCH else "../../" if p==paths[-1] else ""
    heading="## 965：原生材料的共同电磁读口"
    oldheading="## 964：同一材料的膨胀、记录与热账"
    assert heading not in s and s.count(oldheading)==1
    block=(heading+"\n\n"
       +f"[965报告]({pre}research_note_965.md)保持958原材料和接收效果，"
       +"同一电荷接入量子电磁模式及自极化；原记录差>0.998723，"
       +"极化修正场效果差>0.049871，使用无限占据残差界而非只比较截断。"
       +f"[主结果]({pre}965/material_field_window_results.json) · "
       +f"[场读口]({pre}965/physical_mode_effect_results.json) · "
       +f"[核验]({pre}965/research_round_965_checks.json)。正式965／累计3750，整体目标未完成。\n\n"
       +f"[共同恢复v0.8]({pre}965/drafts/common_recovery_v0_8.md)保留单模、表示与读取权限边界；"
       +"不是空间传播或实体光子探测，不能直接合并964不同场准备。"
       +f"停止模式和器件优化，接[966]({pre}966/drafts/STATUS.md)回同一交互、尺度与物理来源。"
       +"应用目标及任务设置保持。\n\n")
    s=s.replace(oldheading,block+oldheading,1)
    s=s.replace("001—964轮共964份","001—965轮共965份")
    if p==paths[4]:
        s=s.replace("当前正式964／累计3749，964已结项","当前正式965／累计3750，965已结项")
    if p==paths[5]:
        s=s.replace("# 231—964轮阶段成果总览","# 231—965轮阶段成果总览",1)
        s=s.replace("231—964的734份","231—965的735份",1)
    if p==paths[-1]:
        s=s.replace("## 六条共同协议：全局缺口对应与检验优先级（截至964）",
                    "## 六条共同协议：全局缺口对应与检验优先级（截至965）",1)
        rows={
         "C03":"|C03 事件身份、关系记录、访问|1、2、3、4；391、434、923、929、947、956—965|965原接收效果与含极化的固定场效果均有正概率差下界|模式参考和末读权限输入；联合效果不是远处裸光子计数，未认证顺序双读|",
         "C14":"|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—965|965原Q给非对易量子电磁模式接口，保自极化与一致表示字典|有限材料和单模物理输入；不等于完整QED规范匹配或Maxwell光锥|",
         "C20":"|C20 跨尺度、误差与共同来源|4、6、H3；592、625、854、898、955—965|965用真实边界与有限谱残差控制无限光子占据的有界读口|声明初态／单模域；不代替未选模式、真实材料及无界应力的误差|",
         "C21":"|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—965|965同一旧材料保持记录并把内容接到电磁模式，无新基本记录荷|没有实体探测、空间飞行或整套六协议；停止局部模式工程|",
         "C22":"|C22 应力、守恒与反作用|1、3、4、6、H2；935—965|965同一自治H含材料、模式、线性作用、自极化和参数源|内部能源闭合不等于腔体支持和全几何源；964自由压力不自动用于新交互|",
         "C27":"|C27 共同数据与可区分预测|1、3、5、6；945、958—965|965固定参数有受控原记录与极化修正场效果的共同窗口|单模有效存在见证，不是实验／全共同恢复；966回同一尺度来源规则|"}
        for key,row in rows.items():
            s,count=re.subn(r"^\|"+key+r" [^\n]*$",lambda _:row,s,count=1,flags=re.M)
            assert count==1
    if b"\r\n" in raw:s=s.replace("\n","\r\n")
    updates[p]=(b"\xef\xbb\xbf" if raw.startswith(b"\xef\xbb\xbf") else b"")+s.encode("utf-8")
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print("Published eight living navigation files and the 965 preservation verifier.")

