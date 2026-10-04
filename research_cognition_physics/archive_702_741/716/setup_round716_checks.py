"""Generate716 checks from the previous exclusive publication workflow."""
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,content):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(content)


v=(HERE/'verify_round715.py').read_text('utf8')
# Simultaneous substitutions protect the preceding-round references.
mapping={'715':'716','714':'715','716':'717','3002':'3016','3016':'3031',
         '3349':'3352','1457':'1460',
         'joint_sharp_record_domain':'joint_smooth_mode_contract',
         'joint_center_thermal_reference':'joint_sharp_record_domain',
         "text['display_formulas']==16":"text['display_formulas']==20",
         'observable_gain_entry':'sterile_mode_entry'}
import re
v=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
         lambda m:mapping[m.group()],v)
v=v.replace("           'round716_drafts/entry_checks.json')",
            "           'round716_drafts/entry_checks.json','round716_drafts/first_action_diagnostic.json')")
write('verify_round716.py',v)

p=(HERE/'publish_round715.py').read_text('utf8')
mapping={'715':'716','714':'715','716':'717','3349':'3352','3346':'3349',
         '1457':'1460','3016':'3031','361.':'362.','266.':'267.',
         'joint_sharp_record_domain':'joint_smooth_mode_contract'}
p=re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
         lambda m:mapping[m.group()],p)
old=next(line for line in p.splitlines() if line.startswith('summary='))
new="summary='**第716轮完成：** [光滑物质模的共同资源与来源]({p}research_note_716.md)原模式的局部势矩与跳跃共同控制读取、时间变化及指定几何导数；非均匀体积改变仪器，来源需同步。三组、二十式通过，最新716／3352，1460份编号科学文件、3031份保护证据。[核验]({p}research_round_716_checks.json)、[全条件账]({p}unified_physics_condition_ledger_716.md)。实际共同连续与因果实现仍开放。'"
p=p.replace(old,new)
old=next(line for line in p.splitlines() if line.startswith('order='))
new="order='**当前执行顺序（716后，优先于下方历史安排）：** 接[717原记录后参考与共同传播]({p}round717_drafts/STATUS.md)，先核实际非平稳态和有限时间过程，复用623—625、632—637及667。停止波包／范数优化，旧空间和699范围保持。'"
p=p.replace(old,new)
p=p.replace('原真实热态中的记录收敛与能源域','原物质模的体积、资源与几何共同条件')
p=p.replace('旧空间合同复用，锐化极限不代替空间连续','旧空间合同复用，波包归一不选择物理维数')
p=p.replace('[实际热读取与共同能源域](research_note_716.md)',
            '[光滑物质模的共同资源与来源](research_note_716.md)')
write('publish_round716.py',p)

p=(HERE/'postcheck_round715.py').read_text('utf8')
mapping={'715':'716','716':'717','3016':'3031'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
write('postcheck_round716.py',p)
print('Generated716 checks and exclusive publication scripts.')
