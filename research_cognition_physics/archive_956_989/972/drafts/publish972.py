"""Publish 972 once, preserving all previous scientific files."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'971/verify_round971.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"native_covariant_source_results.json")',1)[0]
for a,b in [('range(776,971)','range(776,972)'),('Delivery checks for 971','Delivery checks for 972'),
            ('research_round_971_checks.json','research_round_972_checks.json')]:prefix=prefix.replace(a,b)
checks=r'''    result=read(HERE/"native_thermal_record_results.json")
    assert result["round"]==972 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    assert [x["commutant_dimension"] for x in result["algebra"]]==[2,1]
    assert result["algebra"][0]["sector_commutator"]==0
    assert result["algebra"][1]["sector_commutator"]>.02
    for row in result["rows"]:
        assert abs(sum(row["energy_changes"]))<1e-11
        assert abs(sum(row["energy_currents"]))<1e-12
        assert row["landauer_identity_error"]<1e-12
        assert row["entropy_production"]>row["entropy_production_error_bound"]
        assert row["certificate"]["material_trace_distance_bound"]<4e-8
        if row["eta"]==0:assert abs(row["record_information_loss"])<1e-12
    assert result["witness"]["information_loss_lower"]>2.15e-5
    assert result["rows"][3]["field_heat"]+result["rows"][3]["heat_error_bound"]<0
    assert abs(result["rows"][-1]["field_heat"])<result["rows"][-1]["heat_error_bound"]
    assert abs(result["isolated_material_control_information"]-math.log(2))<1e-12
    for key in ("complete_reset_proved","thermalization_proved",
                "monotone_irreversible_arrow_proved","cosmology_common_window_proved","full_goal_completed"):
        assert result["scope"][key] is False
    import importlib.util
    spec=importlib.util.spec_from_file_location("core972",HERE/"native_thermal_record.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    note=STAGE/"research_note_972.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==16
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,9)]
    for term in ("整体目标未完成","不领取正热流结论","不是材料熵减少的擦除过程",
                 "不是整个24维材料没有守恒量","没有证明完整复位"):
        assert term in prose,term
    newfiles=[note,HERE/"native_thermal_record.py",HERE/"native_thermal_record_results.json",
        Path(__file__),HERE/"drafts/thermal_record_decision.md",HERE/"drafts/probe_contact.py",
        HERE/"drafts/thermal_record_mechanism.md",HERE/"drafts/publish972.py",
        STAGE/"973/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
'''
tail='    nav=[STAGE.parent/n for n in'+old.split('    nav=[STAGE.parent/n for n in',1)[1]
for a,b in [
    ('截至971','截至972'),('n<=971','n<=972'),('list(range(1,972))','list(range(1,973))'),
    ('001—971轮共971份','001—972轮共972份'),('231—971的741份','231—972的742份'),
    ('971：原生材料的共同父接口与重复计数检验','972：原生记录的热接触、保护扇区与内部交换'),
    ('STAGE/"970/research_round_970_checks.json"','STAGE/"971/research_round_971_checks.json"')]:
    tail=tail.replace(a,b)
start=tail.index('    out=dict(');end=tail.index('    if not writing:',start)
tail=tail[:start]+'''    out=dict(round=972,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=972,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_971=prev["cumulative_numbered_test_groups_from_970"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,
        arbitrary_electric_bath_sector_obstruction=True,
        old_native_exchange_lifts_local_protection=True,
        common_infinite_mode_finite_window_error_bound=True,
        finite_information_energy_entropy_checked=True,
        complete_reset_verified=False,irreversible_arrow_verified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
'''+tail[end:]
target=STAGE/'972/verify_round972.py';assert not target.exists()
target.write_text(prefix+checks+tail,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 972：原生记录的热接触、保护扇区与内部交换'
    previous='## 971：原生材料的共同父接口与重复计数检验'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[972报告]({pre}research_note_972.md)证明原生纯电荷作用对任意环境保留P标签，不能将两扇区共同重置。'
        +'968已有内部交换解除局部保护，与同型电磁模式共同作用后，有限记录信息转移下界>2.15×10⁻⁵，'
        +'完整能源及热参考熵账闭合；没有完全擦除或不可逆箭头。'
        +f'[结果]({pre}972/native_thermal_record_results.json) · '
        +f'[核验]({pre}972/research_round_972_checks.json)。正式972／累计3757，整体目标未完成。\n\n'
        +f'[机制增量]({pre}972/drafts/thermal_record_mechanism.md)把保护、更新与热资源约束回同一材料作用，'
        +'不把局部标签信息减少直接当成放热或资料全局消失。'
        +f'停止热接触效率与器件优化，接[973]({pre}973/drafts/STATUS.md)核对一份共同实现的对象、准备和重叠任务域。'
        +'应用目标和957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—971轮共971份','001—972轮共972份')
    if p==paths[4]:s=s.replace('当前正式971／累计3756，971已结项','当前正式972／累计3757，972已结项')
    if p==paths[5]:
        s=s.replace('# 231—971轮阶段成果总览','# 231—972轮阶段成果总览',1)
        s=s.replace('231—971的741份','231—972的742份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至971）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至972）',1)
        rows={
        'C03':'|C03 事件身份、关系记录、访问|1、2、3、4；391、434、929、958—972|972区分严格P标签保护、局部信息进入关联与全局标签保持|只有指定总零材料部门；没有完全复位或信息全局消失|',
        'C19':'|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—957、964—972|972同一材料与一次有限能Gibbs模式准备有真实有限交互及来源|不是无限可重置热浴；准备成本输入，未进入969相同历史|',
        'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；430—459、929、956—972|972以968既有F和原Q解除局部保护，实现有界的热记录信息转移|未制造复位器件或实际读口；不等于全24维材料完全可控|',
        'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—972|972材料、模式和自极化／交互同H保能源，不能把材料能差直接当热|971协变身份范围保留；本轮未合并动态支撑或宇宙几何|',
        'C23':'|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—972|972真实热接触的有限信息转移与参考熵账闭合；旧F解除纯电荷复位禁阻|没有完整擦除、热化或单调箭头；特殊准备和有限容量仍需明示|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—972|972无F标签信息精确保持；有F与场共同作用后减少>2.15×10⁻⁵|量很小且未作器件效率声明；不是实验确认或全部门共同预测|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[961整体图](../../research_note_961.md)、[969共同恢复v0.9](../../969/drafts/common_recovery_v0_9.md)及[971父接口](../../971/drafts/native_parent_map_v0_1.md)保留。[972机制增量](../../972/drafts/thermal_record_mechanism.md)使保护、更新与热接触共享实际材料作用。停止局部器件优化，接[973](../../973/drafts/STATUS.md)核对一份共同作用、状态与重叠任务域，不继续轮换小模型。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published eight living navigation files and the 972 verifier.')
