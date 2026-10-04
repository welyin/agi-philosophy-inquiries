"""Create715 checks/publication while preserving all existing scripts."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def save(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


def rounds(text):
    return text.replace('714','__R__').replace('715','716').replace('713','714').replace('__R__','715')


text=rounds((HERE/'verify_round714.py').read_text('utf8'))
text=text.replace('joint_center_thermal_reference','joint_sharp_record_domain')
text=text.replace('joint_holonomy_readout_scale','joint_center_thermal_reference')
text=text.replace('len(history)==2986','len(history)==3002')
text=text.replace('len(history|new|preserved)==3002','len(history|new|preserved)==3016')
names="""names=('unified_physics_condition_ledger_715.md','round715_drafts/research_note_715_draft.md',
           'round715_drafts/final_review.txt','round715_drafts/literature_scope_audit.json',
           'round715_drafts/scope_and_dedup_review.md','round716_drafts/STATUS.md',
           'round715_drafts/observable_gain_entry.py','round715_drafts/observable_gain_entry_results.json',
           'round715_drafts/observable_gain_entry.md','round715_drafts/check_and_publish_entry.py',
           'round715_drafts/entry_checks.json')
    new="""
text=re.sub(r'names=\(.*?\)\n    new=',lambda _:names,text,count=1,flags=re.S)
text=text.replace('cumulative_numbered_tests=3346','cumulative_numbered_tests=3349')
text=text.replace('cumulative_numbered_scientific_files=1454','cumulative_numbered_scientific_files=1457')
text=text.replace('unchanged_prior_evidence_files=2986','unchanged_prior_evidence_files=3002')
text=text.replace('cumulative_unique_protected_evidence_files=3002','cumulative_unique_protected_evidence_files=3016')
save('verify_round715.py',text)

text=rounds((HERE/'publish_round714.py').read_text('utf8'))
text=text.replace('joint_center_thermal_reference','joint_sharp_record_domain')
summary="**第715轮完成：** [实际热读取与共同能源域]({p}research_note_715.md)指定环路锐化在原完整Gibbs上的平均注能为线性阶；记录／后态虽收敛，极限却有无限能源，共形来源同步发散。三组、十六式通过，最新715／3349，1457份编号科学文件、3016份保护证据。[核验]({p}research_round_715_checks.json)、[全条件账]({p}unified_physics_condition_ledger_715.md)。仅关闭固定图锐化接法，统一目标未完成。"
order="**当前执行顺序（715后，优先于下方历史安排）：** 接[716光滑物质模与同一原过程]({p}round716_drafts/STATUS.md)，复用633及666—667，核原有限图sterile读取、共同参考、质量与来源；停止锐化函数扫描，保留自由分支与因果实现限制。旧空间及699范围保持。"
text=re.sub(r"summary='[^\n]*'",lambda _:'summary='+repr(summary),text,count=1)
text=re.sub(r"order='[^\n]*'",lambda _:'order='+repr(order),text,count=1)
text=text.replace('## 360.','## 361.').replace('## 265.','## 266.')
text=text.replace('原热参考、带荷差量与几何响应','原真实热态中的记录收敛与能源域')
text=text.replace('旧空间合同复用，热差量不反推空间维数','旧空间合同复用，锐化极限不代替空间连续')
text=text.replace('中心信息与共同热参考','实际热读取与共同能源域')
text=text.replace("'最新科学轮次与检查数为714／3343'","'最新科学轮次与检查数为714／3346'")
text=text.replace("'最新科学轮次与检查数为715／3346'","'最新科学轮次与检查数为715／3349'")
text=text.replace('cumulative_tests=3346','cumulative_tests=3349')
text=text.replace('numbered_scientific_files=1454','numbered_scientific_files=1457')
text=text.replace('unique_protected_evidence_files=3002','unique_protected_evidence_files=3016')
save('publish_round715.py',text)

text=rounds((HERE/'postcheck_round714.py').read_text('utf8'))
text=text.replace('==3002','==3016')
save('postcheck_round715.py',text)
print('Created715 verification and publication scripts.')
