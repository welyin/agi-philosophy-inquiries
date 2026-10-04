"""Prepare725 evidence; preserve historical science and current goal."""
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


ledger=(HERE/'unified_physics_condition_ledger_724.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：725原关系壁的闭法向平方与同一记录来源\n'+rest
ledger=ledger.replace('接[723全账](unified_physics_condition_ledger_723.md)，回填[724报告](research_note_724.md)。[结果](joint_material_region_records_results.json)、[核验](research_round_724_checks.json)。',
                      '接[724全账](unified_physics_condition_ledger_724.md)，回填[725报告](research_note_725.md)。[结果](joint_relational_normal_form_results.json)、[核验](research_round_725_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**725当前增量：** 原曲目标全域轴点控制给原壁闭导数平方及完整能源界，不先要求带符号速度自伴。同一实际选区的法向矩、规范输送、有限过程及一阶来源接通；位置读取一般降低同一Q的非选择均值。直接Q仪器、经典因果类型及引力拼接仍分范围。\n\n'+marker)
updates={
    'C03':'725实际原位置记录保有限法向矩，但一般改变原Q；不能沿用读前全部几何资料',
    'C09':'725原h壁由原核闭包的C-dagger-C定义法向平方，保s/h混合项，带符号速度自伴非前置',
    'C19':'725原全域Hardy映射使同一原能源控制轴点及法向绝对矩，原Gibbs无需换态',
    'C20':'725闭法向形式沿724实际区域、617匹配及723能源夹权有限过程共同运输',
    'C22':'725法向形式本身变化、准备变化和实际标签端点均入一阶来源'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 725原壁法向平方、记录反作用与来源

- 原曲目标的全域Hardy余项控制F/h²，轴外原核在完整能源形式域稠密；722局部工具不充当全域常数。
- 原一阶导数闭包C受原完整H能源控制，C-dagger-C给正自伴平方；不假定C自身自伴，不宣称物理量子化唯一。
- 同一精确h壁的空间差分及s/h混合项全部保留，Q为自伴候选；原正常有限能源态有有限绝对矩。
- 724位置instrument对法向平方有精确有限反作用。它否定指定读取不扰动法向的接法，不否定联合有限精度或半经典任务。
- 原区域匹配、同次标签合并和能源夹权有限过程保实际后态上的终端法向形式；一阶来源保Q、准备及端点三项。
- 直接Q谱读取的能源保持未证明，722的自伴速度群论证不能直接挪用。未得经典因果类型、空间细化或量子Einstein拼接。

## 本轮合并与下一项

C03／C09／C19／C20／C22的原精确壁、闭法向、实际区域、能源及来源共同实现，减去“先证明带符号速度自伴”这一不必要前置；闭包及差分仍明示为候选定义。

接[726](round726_drafts/STATUS.md)：回查574／575、588及647—651原同一经典来源，检验实际选区与法向在同一受控映射能否共同恢复原严格壁类型。停止更多谱函数、阈值和反作用优化；任何半经典参数或集合尺度均标输入，649、699及旧空间范围保持，统一目标开放。
'''
write('unified_physics_condition_ledger_725.md',ledger)
write('round725_drafts/research_note_725_draft.md',(HERE/'research_note_725.md').read_text('utf8'))
write('round726_drafts/STATUS.md','''# 第726轮入口：同一原来源、实际选区与经典壁

接[725](../research_note_725.md)、[全账](../unified_physics_condition_ledger_725.md)。原精确h壁的闭法向平方和实际区域已有共同能源、后态及一阶来源；位置读取一般改变法向资料。停止单一读口与域变种优化，转回原真实几何任务。

1. 回查574／575原半经典Gauss来源、588局部流及647—651的同一经典解、类时壁和完整边界资料。先辨认哪些局部／有限图量已经连接，勿重做一般相干态极限。
2. 核同一个状态族与同一次允许仪器中，字段位置误差、法向平方、混合项、测量反作用及原来源能否共同控制。
3. 若采用ℏ→0或指定尺度转换，明示这是原模型已有或新增的受控极限，不能宣称已得固定ℏ宇宙的物理连续极限。不得更换原壁为相切线性T壁。
4. 位置记录与Q的终端形式可复用725；实际直接Q读法仍不能自动沿用722能源群界。先满足任务所需范围，不为函数变种独立编号。
5. 382—386、425、522—524旧空间、649引力匹配及699限定结论保留；不存在从字段标签或均值符号直接推出完整时空的步骤。统一目标保持。
''')
write('round725_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://arxiv.org/pdf/0803.0503',authors='Frank and Seiringer',
                  read='Introduction equation(1.1), local Hardy inequality.',
                  scope='Mature square-completion method only; no fractional theorem or curved-target theorem imported.')],
    inherited='556/563 radial domains;722 local Hardy;652 original smooth references and W prescription;617 gauge gluing;704/723 energy-weighted process;724 actual material regions.',
    own_mapping='Global curved-target F/h^2 remainder; original-core closed square controlled by full energy; same region measurement normal backaction and terminal source transport.',
    not_claimed='Self-adjoint signed wall velocity; uniqueness of physical ordering; Q spectral instrument energy preservation; classical causal signature; quantum Einstein gluing.'
),ensure_ascii=False,indent=2)+'\n')
write('round725_drafts/scope_and_dedup_review.md','''# 725范围与去重

主代理审查；无新增独立代理审查。

- 725入口只审计旧556／563域警告的原壁映射，不计正式科学增量。
- 722已有局部Hardy；本轮具体原K给全域F/h²恒等式，再接全部原能源及实际区域法向。
- 初始核轴外稠密需四维横向容量估计，非零测度论。C可闭由核对称得到；C-dagger-C自伴不要求C自伴。
- 所选闭包是明确候选定义，未排除其他量子化。空间差分仍是652处方，原H不更换。
- Gaussian后法向变化是同一原壁同一读口的具体共同任务结果，不称通用测量下界。
- 解析保全CAR／Gauss，连续积分只校准径向族，64维用旧H和另明微分形式，非全图谱。
- 终端形式表示和实际Q测量不同；722仪器能源结论不自动继承。
- 原几何外部，未得经典严格符号、时空连续、649引力商空间或699修复。旧空间假设不恢复。
''')
names=('research_note_725.md','joint_relational_normal_form.py','joint_relational_normal_form_results.json',
       'unified_physics_condition_ledger_725.md')
write('round725_drafts/final_review.txt',
      'Main-agent review only. Global original-measure Hardy identity; axis-core energy density; closable symmetric derivative; C-dagger-C closed form; original full-energy bound; actual Gaussian normal bias; matched-range and terminal-source transport checked. Direct Q instrument energy preservation and causal signature are not claimed.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
v=replace((HERE/'verify_round724.py').read_text('utf8'),{
    '724':'725','723':'724','725':'726','3130':'3144','3144':'3158','3377':'3380','1484':'1487',
    'joint_material_region_records':'joint_relational_normal_form',
    'joint_reference_process_transport':'joint_material_region_records',
    'relational_wall_entry':'wall_domain_entry'})
write('verify_round725.py',v)
p=replace((HERE/'publish_round724.py').read_text('utf8'),{
    '724':'725','723':'724','725':'726','3377':'3380','3374':'3377',
    '1484':'1487','3144':'3158','370.':'371.','275.':'276.',
    'joint_material_region_records':'joint_relational_normal_form'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第725轮完成：** [原关系壁的闭法向平方与同一记录来源]({p}research_note_725.md)原曲目标全域轴点控制给原壁闭平方和完整能源界；同一位置记录保有限法向资料但改变其均值，区域、有限过程及一阶来源共同保持。三组、十八式通过，最新725／3380，1487份编号科学文件、3158份保护证据。[核验]({p}research_round_725_checks.json)、[全条件账]({p}unified_physics_condition_ledger_725.md)。直接Q仪器与经典壁类型尚未完成。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（725后，优先于下方历史安排）：** 接[726同一来源、选区与经典壁]({p}round726_drafts/STATUS.md)，回查574／575及647—651，核同一受控极限中的位置、法向、测量反作用与来源。停止谱函数、阈值和域变种优化；旧空间、649／699及统一目标保持。'")
p=p.replace('原物质区域、规范组织与移动来源','原壁闭法向、实际记录与共同来源')
p=p.replace('旧空间合同复用，实际区域仍须接法向与因果类型','旧空间合同复用，量子法向仍须接原经典壁')
p=p.replace('[原物质选区、规范拼接与移动记录来源](research_note_725.md)',
            '[原关系壁的闭法向平方与同一记录来源](research_note_725.md)')
write('publish_round725.py',p)
write('postcheck_round725.py',replace((HERE/'postcheck_round724.py').read_text('utf8'),{
    '724':'725','725':'726','3144':'3158'}))
print('Prepared725 evidence and726 common-source task.')
