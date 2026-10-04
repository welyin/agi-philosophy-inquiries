"""Preserve726 joint classical-wall result and full-matter limitation."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_725.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：726同源记录后的经典壁与完整物质涨落边界\n'+rest
ledger=ledger.replace('接[724全账](unified_physics_condition_ledger_724.md)，回填[725报告](research_note_725.md)。[结果](joint_relational_normal_form_results.json)、[核验](research_round_725_checks.json)。',
                      '接[725全账](unified_physics_condition_ledger_725.md)，回填[726报告](research_note_726.md)。[结果](joint_recorded_classical_wall_results.json)、[核验](research_round_726_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**726当前增量：** 同源574 Gauss波包和724真实位置读取在固定图窗口σ→0、ℏ/σ→0下共同保位置、读后法向谱及玻色来源。原差分误差明示；598空CAR参考虽保原能源／来源均值，却留下正Majorana能量涨落，不满足完整来源的无涨落经典接法。下一项转物质参考与传播相容性。\n\n'+marker)
updates={
    'C03':'726同一真实K_y后态的有限玻色菜单残差共同消失，大概率记录条件态集中，非联合锐测量',
    'C09':'726原闭Q读后谱重量有错误符号上界；原651壁保持且不替换前向差分',
    'C15':'726原Majorana在空CAR参考留下严格正完整H涨落；未删耦合，其他参考未排除',
    'C19':'726原Gauss半经典来源可复用，但574玻色包乘空CAR不自动给完整物质经典参考',
    'C20':'726只证明固定图半经典及后续经典采样顺序，不签收固定ℏ共同空间极限',
    'C22':'726原完整瞬时来源均值可匹配，完整H涨落和变化准备总导数需另核'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 726真实记录、条件性经典壁与原物质涨落

- 574已有同源严格Gauss管状波包，不重做一般半经典存在；724 K_y保真实记录，725 Q在原h>0管内与局部表达一致。
- Gaussian一二阶核导数给读后玻色有限菜单的共同平方残差；σ→0且ℏ/σ→0同时压低位置误差与法向谱误判上界。
- 非交换菜单以同一条件态的多个谱尾描述，不宣布联合锐测量；壁上节点不强行确定到某一侧。
- 原前向差分保留，651原点给精确正符号与明确采样误差。先固定图再经典采样，不交换两极限。
- 598原空CAR参考对完整质量／跳跃均值为零，但Majorana配对给正能源方差。旧712／719机制用于当前共同极限，不称为新配对发现。
- 原完整H及瞬时来源均值匹配不足推出完整物质无涨落极限、有限等待或Einstein自洽；换参考和联合态仍可检验。

## 本轮合并与下一项

共同记录后的经典几何参考取得固定图条件连接，同时暴露C15／C19／C22费米参考的真实兼容要求。停止波包宽度、阈值、法向读口或差分方案优化。

接[727](round727_drafts/STATUS.md)：先回查602—604、614、630—633、666—667及719—720，检验同一完整物质参考能否保原来源、实际记录与空间传播。换态须交代准备／来源，不以删Majorana或重选空真空解决。649、699及旧空间范围保留，统一目标开放。
'''
write('unified_physics_condition_ledger_726.md',ledger)
write('round726_drafts/research_note_726_draft.md',(HERE/'research_note_726.md').read_text('utf8'))
write('round727_drafts/STATUS.md','''# 第727轮入口：原物质参考、几何来源与空间传播

接[726](../research_note_726.md)、[全账](../unified_physics_condition_ledger_726.md)。原真实位置记录与玻色法向／来源可共同接回声明的经典壁；空CAR参考在保原均值时留下非消失Majorana涨落。因此优先核原物质参考，不再优化波包或阈值。

1. 回查602—604的原物质参考、614字典、630—633连续自由分支、666—667变换以及719—720共享／自旋参考。已有真空非定态、粒子空穴区别和保护自旋限制直接复用。
2. 检验是否存在合格原Gauss费米参考或联合态，在同一玻色来源上保持目标几何来源，并控制完整能源／应力涨落；区分共同核态、基态、热态与半经典态。
3. 任何保护占据／自旋构造需同步核原空间跳跃、Weyl传播及原Dirac／Majorana。不能将某固定图spin标量分支直接升级为完整手征物理。
4. 换参考可能改变能源、几何及准备代价，必须同账；不可把换字典后的空真空当原态，不删质量项也不默认所有参考失败。
5. 保382—386、425、522—524旧空间、649引力匹配、699限定路线边界；不设置新目标，不先做装置工程。仅真正关闭共同接口或取得限定反例才编号。
''')
write('round726_drafts/literature_scope_audit.json',json.dumps(dict(
    external_sources_newly_imported=[],
    inherited='574 original Gauss semiclassical tube proof;651 original causal wall;724 actual Gaussian instrument;725 original normal form;712/719 original Majorana vacuum action.',
    own_mapping='Actual postmeasurement joint bosonic residual estimate and original wall difference bound, with full original CAR variance floor retained.',
    excluded='General coherent-state existence counted again; noncommuting joint sharp measurement; full-matter fluctuation-free limit inferred from bosonic concentration; finite-hbar spatial continuum.'
),ensure_ascii=False,indent=2)+'\n')
write('round726_drafts/scope_and_dedup_review.md','''# 726范围与去重

主代理审查，无新增独立代理。

- 574态存在及读前O(hbar)平方残差直接继承；新增为实际K_y完整后态及同一原壁任务，不以隐藏随机相位当记录。
- 二阶算符交换子只在原紧h>0/F>0管上使用；K_y保支撑。全图常数不声称一致。
- 菜单为标量玻色二阶符号；598原CAR不能不经检验地写成相同标量极限。
- 错符号概率是终端Q谱重量上界，不是已造锐Q仪器。罕见记录、壁上左右分类和非交换菜单均留范围。
- 既定前向差分误差保留，数值固定同一点正几何不是重新完成全图量子细化。
- 712／719已知Majorana空参考非定态，本轮把正方差准确接回原完整H与半经典族；不称所有参考失败。
- 本轮固定原CAR与Yukawa的hbar族为明确输入；几何来源只是固定canonical瞬时均值，不借704替代波包准备总导数。
- 数值为原局部径向条件诊断、解析原场差分和原CAR校准；完整Gauss证明解析，非全图量子数值。
- 旧空间及649／699不变，统一目标未完成。
''')
names=('research_note_726.md','joint_recorded_classical_wall.py','joint_recorded_classical_wall_results.json',
       'unified_physics_condition_ledger_726.md')
write('round726_drafts/final_review.txt',
      'Main-agent review only. Actual K_y second-derivative estimate and postrecord spectral bound; fixed-graph constants; canonical source and forward-difference error; full physical-right CAR vacuum pair norm; retained Majorana variance; limited instantaneous source means reviewed. No autonomous preparation or full quantum Einstein limit.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
v=replace((HERE/'verify_round725.py').read_text('utf8'),{
    '725':'726','724':'725','726':'727','3144':'3158','3158':'3172','3380':'3383','1487':'1490',
    'joint_relational_normal_form':'joint_recorded_classical_wall',
    'joint_material_region_records':'joint_relational_normal_form','wall_domain_entry':'common_wall_limit_entry'})
write('verify_round726.py',v)
p=replace((HERE/'publish_round725.py').read_text('utf8'),{
    '725':'726','724':'725','726':'727','3380':'3383','3377':'3380','1487':'1490','3158':'3172',
    '371.':'372.','276.':'277.','joint_relational_normal_form':'joint_recorded_classical_wall'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第726轮完成：** [同源记录后的经典壁与完整物质涨落边界]({p}research_note_726.md)原真实记录在固定图半经典窗口共同恢复位置、法向谱及玻色来源；空CAR参考保均值却留正Majorana能源涨落，完整物质不能直接签收。三组、十八式通过，最新726／3383，1490份编号科学文件、3172份保护证据。[核验]({p}research_round_726_checks.json)、[全条件账]({p}unified_physics_condition_ledger_726.md)。原差分误差、量子连续及引力边界保留。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（726后，优先于下方历史安排）：** 接[727原物质参考与空间传播]({p}round727_drafts/STATUS.md)，回查602—604、614、630—633及719—720，核费米参考、同一来源、质量与传播能否共同相容。停止波包与法向优化；旧空间、649／699及统一目标保持。'")
p=p.replace('原壁闭法向、实际记录与共同来源','读后经典壁与完整物质参考')
p=p.replace('旧空间合同复用，量子法向仍须接原经典壁','旧空间合同复用，玻色经典壁不等于完整来源')
p=p.replace('[原关系壁的闭法向平方与同一记录来源](research_note_726.md)',
            '[同源记录后的经典壁与完整物质涨落边界](research_note_726.md)')
write('publish_round726.py',p)
write('postcheck_round726.py',replace((HERE/'postcheck_round725.py').read_text('utf8'),{
    '725':'726','726':'727','3158':'3172'}))
print('Prepared726 publication and727 full-matter task.')
