"""Preserve724 evidence and prepare conflict-safe navigation publication."""
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


ledger=(HERE/'unified_physics_condition_ledger_723.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：724原物质选区、规范拼接与移动记录来源\n'+rest
ledger=ledger.replace('接[722全账](unified_physics_condition_ledger_722.md)，回填[723报告](research_note_723.md)。[结果](joint_reference_process_transport_results.json)、[核验](research_round_723_checks.json)。',
                      '接[723全账](unified_physics_condition_ledger_723.md)，回填[724报告](research_note_724.md)。[结果](joint_material_region_records_results.json)、[核验](research_round_724_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**724当前增量：** 保651原精确h壁，用同一次有限精度字段记录选择节点区域；全Gauss匹配、原能源及来源随记录运输。合并同次记录与拼回区域交换；移动记录边界必须计入分支来源。未得到时空法向、量子引力拼接、自主装置或空间细化。\n\n'+marker)
updates={
    'C03':'724同次字段记录给真实有限区域instrument；区间效果的Luders替换一般改变后态',
    'C09':'724保651精确h壁实现节点选区，法向／因果类型尚须共同导数域',
    'C19':'724实际记录保原完整平均能源；分辨率为任务输入，记录装置代价未被消去',
    'C20':'724同次记录合并与Gauss区域拼接交换，非物理图细化',
    'C22':'724移动读数区间的来源通量逐分支保留，拉回同一原空间后非选择相消'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 724同一物质记录、区域组织与来源

- 563弱径向、617全规范匹配、647—651经典边界、652参考、704／723热准备与有限实际过程直接继承。
- 原f=s−λh保持实际壁；h的弱梯度给完整原能源注能恒等式，没有宣布速度自伴或h全局光滑。
- 同次Gaussian字段读数经有限区间分类，保实际CP后态；全CAR、全部边界匹配和补集保留。
- 同次细记录合并与原Gauss区域拼回严格交换，并保原能源／来源；不是重复测量或图细化。
- 阈值变化贡献真实分支来源；非选择相消只在拉回同一原空间后成立，不等于GHY或完整引力边界源。
- 649重复几何动能降域障碍、699限定路线结果以及旧空间合同保持。

## 本轮合并与下一项

C03／C09／C19／C20／C21／C22的原物质字段、真实区域、规范边界资料、能源和一阶来源接入同一对象。停止区间精度及阈值优化。

接[725](round725_drafts/STATUS.md)：核原壁位置读取与法向／因果类型能否在同一量子域和同一实际任务实现；先复用563、651—652及723—724，不从h可测性直接断言h速度自伴。空间接口382—386、425、522—524保持；统一目标开放。
'''
write('unified_physics_condition_ledger_724.md',ledger)
write('round724_drafts/research_note_724_draft.md',(HERE/'research_note_724.md').read_text('utf8'))
write('round725_drafts/STATUS.md','''# 第725轮入口：实际关系区域与原壁法向

接[724](../research_note_724.md)、[全账](../unified_physics_condition_ledger_724.md)。同次物质记录已能划分节点区域并一致合并，保原Gauss、能源及一阶记录来源。这不自动给时空法向。

1. 回查563弱径向导数、651原h壁、652光滑T／s速度及平方域；保原f=s−λh，不改成相切的线性T壁。
2. 核“位置函数可测／保能源”与“其法向速度或平方可作为同一自伴参考”之间是否已有证明。h轴的弱导数足够724，但可能不足完成流／自伴性；先查旧轮再检验。
3. 若存在域障碍，明确是否仅阻断直接h速度读法，能否由旧T／s与混合资料在h>0任务域重建；不混成原认知原则反证。
4. 与724实际区域通道及原完整H共同核后态、法向资料和来源。原几何仍外部，649与699保持；不将标签端点通量称为引力边界方程。
5. 停止Gaussian函数、阈值、尾部或常数扫描。复用382—386、425、522—524，不恢复384已消去的Lipschitz；386／425替代桥不合并为同时必要。
''')
write('round724_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://arxiv.org/pdf/1109.0036',author='William Donnelly',
                  read='Sec. II, original Hilbert embedding and boundary representation data.',
                  scope='Gauge-region embedding only, already inherited through617. No gravity or continuum claim imported.')],
    inherited='563 weak radial method;617 same full gauge energy gluing;647-651 classical relational walls;652 references;704/723 actual own-Gibbs source.',
    own_mapping='Same original wall selects regions through one instrument; record coarsening commutes with original gluing; moving record endpoints supply branch source flux.',
    not_imported='Autonomous apparatus, independent physical-region factorization, quantum spacetime normal, physical graph refinement or quantum Einstein boundary dynamics.'
),ensure_ascii=False,indent=2)+'\n')
write('round724_drafts/scope_and_dedup_review.md','''# 724范围与去重

主代理审查；无新增独立代理。

- 弱径向、Gaussian单次反作用和规范切口分别已有563／617等结果，不称为新基础定理。
- 真正新增对象：原精确关系壁选择的实际区域记录，及同次记录组织合并与全Gauss拼接的交换、能源与移动标签来源。
- 有限图所有区域和原CAR保留。标签分区不是图细化，切边中间未生成连续位置。
- 32维系数测试含被动参考与非阿贝尔切口，不称为全端点Gauss／CAR仿真；能源探针身份与完整H解析身份区分。
- 64维旧诊断测试来源链式关系，不把其有限差分动能等同连续Gaussian注能恒等式。
- 位置h允许弱梯度，不据此声称其速度自伴；法向、因果类型、649引力拼接与699仍保持。
- 旧空间382—386、425、522—524复用，384Lipschitz不恢复。
''')
main=('research_note_724.md','joint_material_region_records.py',
      'joint_material_region_records_results.json','unified_physics_condition_ledger_724.md')
write('round724_drafts/final_review.txt',
      'Main-agent review only. Same exact wall; weak-gradient original-energy identity; actual Gaussian-bin CP map; full matched gauge ranges; same-record coarsening; energy-weighted moving-endpoint derivative reviewed. No quantum spacetime normal or Einstein gluing claim.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')

v=replace((HERE/'verify_round723.py').read_text('utf8'),{
    '723':'724','722':'723','724':'725','3116':'3130','3130':'3144',
    '3374':'3377','1481':'1484',
    'joint_reference_process_transport':'joint_material_region_records',
    'joint_squared_reference_readout':'joint_reference_process_transport',
    'joint_reference_entry':'relational_wall_entry'})
write('verify_round724.py',v)
p=replace((HERE/'publish_round723.py').read_text('utf8'),{
    '723':'724','722':'723','724':'725','3374':'3377','3371':'3374',
    '1481':'1484','3130':'3144','369.':'370.','274.':'275.',
    'joint_reference_process_transport':'joint_material_region_records'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第724轮完成：** [原物质选区、规范拼接与移动记录来源]({p}research_note_724.md)原精确物质壁经同次有限读数选择区域；记录合并与全Gauss拼接交换，保原能源及来源；移动记录边界须计分支通量。三组、十八式通过，最新724／3377，1484份编号科学文件、3144份保护证据。[核验]({p}research_round_724_checks.json)、[全条件账]({p}unified_physics_condition_ledger_724.md)。未完成法向、空间细化或量子引力拼接。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（724后，优先于下方历史安排）：** 接[725实际区域与原壁法向]({p}round725_drafts/STATUS.md)，回查563及651—652的导数域，核原h壁位置与法向资料能否共同实现。停止区间、阈值和读口优化；旧空间、649／699范围及统一目标保持。'")
p=p.replace('完整物质参考、实际历史与来源','原物质区域、规范组织与移动来源')
p=p.replace('旧空间合同复用，完整内部过程仍须接实际区域','旧空间合同复用，实际区域仍须接法向与因果类型')
p=p.replace('[完整物质参考记录与同一有限量子过程](research_note_724.md)',
            '[原物质选区、规范拼接与移动记录来源](research_note_724.md)')
write('publish_round724.py',p)
write('postcheck_round724.py',replace((HERE/'postcheck_round723.py').read_text('utf8'),{
    '723':'724','724':'725','3130':'3144'}))
print('Prepared724 review, ledger, verification and725 domain audit.')
