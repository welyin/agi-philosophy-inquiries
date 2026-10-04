"""Create 714 verification/publication scripts from the established workflow."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def save(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


def rounds(text):
    return text.replace('713','__LATEST__').replace('714','715').replace('712','713').replace('__LATEST__','714')


text=rounds((HERE/'verify_round713.py').read_text('utf8'))
text=text.replace('joint_holonomy_readout_scale','joint_center_thermal_reference')
text=text.replace('joint_vertex_shared_evolution','joint_holonomy_readout_scale')
text=text.replace('==2971','==2986').replace('==2986\n    assert (HERE/main','==3002\n    assert (HERE/main')
text=text.replace("(4,0,0)","(3,0,0)").replace("display_formulas']==20","display_formulas']==16")
names="""names=('unified_physics_condition_ledger_714.md','round714_drafts/research_note_714_draft.md',
           'round714_drafts/final_review.txt','round714_drafts/literature_scope_audit.json',
           'round714_drafts/scope_and_dedup_review.md','round715_drafts/STATUS.md',
           'round714_drafts/center_reference_entry.py','round714_drafts/center_reference_entry_results.json',
           'round714_drafts/center_reference_entry.md','round714_drafts/center_cross_gradient_check.py',
           'round714_drafts/center_cross_gradient_results.json','round714_drafts/check_and_publish_entry.py',
           'round714_drafts/entry_checks.json')
    new="""
text=re.sub(r'names=\(.*?\)\n    new=',lambda _:names,text,count=1,flags=re.S)
text=text.replace('fresh_tests=dict(run=4','fresh_tests=dict(run=3')
text=text.replace('cumulative_numbered_tests=3343','cumulative_numbered_tests=3346')
text=text.replace('cumulative_numbered_scientific_files=1451','cumulative_numbered_scientific_files=1454')
text=text.replace('unchanged_prior_evidence_files=2971','unchanged_prior_evidence_files=2986')
text=text.replace('cumulative_unique_protected_evidence_files=2986','cumulative_unique_protected_evidence_files=3002')
save('verify_round714.py',text)

text=rounds((HERE/'publish_round713.py').read_text('utf8'))
text=text.replace('joint_holonomy_readout_scale','joint_center_thermal_reference')
summary="**第714轮完成：** [中心信息、带荷关联与共同热参考]({p}research_note_714.md)原中心变换的完整H差量仅含带荷边项；原Gibbs平均代价严格正，实际记录、热关联及几何来源受同一约束，来源导数须含参考响应。三组、十六式通过，最新714／3346，1454份编号科学文件、3002份保护证据。[核验]({p}research_round_714_checks.json)、[全条件账]({p}unified_physics_condition_ledger_714.md)。未计算原热环路均值或完成共同连续，目标未完成。"
order="**当前执行顺序（714后，优先于下方历史安排）：** 接[715原热环路信息与共同尺度]({p}round715_drafts/STATUS.md)，核实际观测、原参考与物理尺度映射；不把固定图严格正性当统一细化界，不重复被动性。旧空间及四分支保持。"
text=re.sub(r"summary='[^\n]*'",lambda _:'summary='+repr(summary),text,count=1)
text=re.sub(r"order='[^\n]*'",lambda _:'order='+repr(order),text,count=1)
text=text.replace('## 359.','## 360.').replace('## 264.','## 265.')
text=text.replace('原环路信息、读取预算与正常态电能','原热参考、带荷差量与几何响应')
text=text.replace('旧空间合同复用，读取面积不反推空间维数','旧空间合同复用，热差量不反推空间维数')
text=text.replace('环路读取与状态电能','中心信息与共同热参考')
text=text.replace("'最新科学轮次与检查数为713／3339'","'最新科学轮次与检查数为713／3343'")
text=text.replace("'最新科学轮次与检查数为714／3343'","'最新科学轮次与检查数为714／3346'")
text=text.replace('cumulative_tests=3343','cumulative_tests=3346')
text=text.replace('numbered_scientific_files=1451','numbered_scientific_files=1454')
text=text.replace('unique_protected_evidence_files=2986','unique_protected_evidence_files=3002')
save('publish_round714.py',text)

text=rounds((HERE/'postcheck_round713.py').read_text('utf8'))
text=text.replace('==2986','==3002')
save('postcheck_round714.py',text)
print('Created714 verification and publication scripts.')
