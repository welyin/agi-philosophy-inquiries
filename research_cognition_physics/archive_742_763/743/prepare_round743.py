"""Prepare743 reports and verified publication without rewriting old evidence."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,values):
    return re.sub('|'.join(re.escape(k) for k in sorted(values,key=len,reverse=True)),lambda m:values[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_742.md').read_text('utf8')
ledger='# 联合条件总账：743原量子标量与严格Gauss关联\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[741全账](unified_physics_condition_ledger_741.md)，回填[742报告](research_note_742.md)。[结果](joint_record_gaussian_boundary_results.json)、[核验](research_round_742_checks.json)。',
    '接[742全账](unified_physics_condition_ledger_742.md)，回填[743报告](research_note_743.md)。[结果](joint_native_quantum_record_results.json)、[核验](research_round_743_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**743当前增量：** 在原完整固定图H、正常Gauss玻色包乘CAR真空上，原Majorana耦合产生严格非高斯四点缺陷，其短时系数为|Ys|²Var(x5)/ℏ²。同一方差固定标量—配对关联及纠缠。717读口／来源桥直接复用；不同准备未冒充730连续参考，理想仪器及自主末读仍开放。\n\n'+marker)
updates={'C03':'743原作用内生成非高斯记录关联；实际未知态仪器和末读仍需核',
         'C19':'743正常严格Gauss波包有正标量方差，同一准备决定关联和噪声；不同于730连续参考',
         'C22':'743同一原H保总能源，来源使用同一后态；717读口桥复用，动态引力尚未接通该新准备'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 743原完整Gauss过程的关联资源

- 入口条件准备的标量协方差jet保持；正式补原完整H中的正常Gauss初态，未关掉Hb、Dirac及跳跃。
- CAR真空上的第一步只有原Majorana配对；同一标量方差决定四点缺陷、标量—费米关联及初始纠缠。
- 原正常R=0波包直接复用717；平均配对为零不妨碍四点关联非零。
- 717的sin s质量力和719的一阶配对生成是旧结果，不重新计轮次；新的是严格物理准备中的方差型四点与共同关联。
- 未证明633原n_f的Lüders实现；未把新的有限图准备换称730的连续Hadamard参考。
- 同一自治H与623来源域保持，但没有求解此准备的动态量子几何。

## 本轮合并与下一项

C03／C15／C19／C22共用原质量和正常Gauss准备。742限制不要求新增物种；原相互作用已提供超出二次背景的资源。

接[744](round744_drafts/STATUS.md)：核原sin s末读诱导的实际未知态仪器、联合后态及来源；不只比较均值，不添加理想新装置。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_743.md',ledger)
write('round743_drafts/research_note_743_draft.md',(HERE/'research_note_743.md').read_text('utf8'))
write('round744_drafts/STATUS.md','''# 第744轮入口：原相互作用诱导的实际仪器

接[743](../research_note_743.md)、[条件账](../unified_physics_condition_ledger_743.md)。原完整Gauss过程已有非高斯关联；旧717给Majorana向sin s末读的响应，不能把两个可分均值当成原占据仪器已实现。

1. 复用743正常Gauss波包及同节点sterile偶编码span{vac,pair}，保持原全H、原sin s效果及末读操作输入。
2. 在未知编码态可带外部被动参考的范围内，求实际诱导效果／仪器及共同后态；区分读取配对相干与读取n_f。
3. 优先连接有限时间、完整CP输出、实际来源与反作用；复用577、717—718、623—625的域和历史，避免重做一般仪器定理。
4. 先保全部原物种与Gauss，不假定量子标量仍是独立、无成本、不变的参考。准备与最后读口仍须注明操作输入。
5. 不是优化本波包方差或寿命，不重新研究空间维数。原连续、730参考映射、动态几何及全统一仍开放。
''')
write('round743_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/quant-ph/0404180',locations='SectionIV equation17',
    use='Gaussian Wick identity only; nonGaussian generation in original full-Gauss preparation is derived here.')],
    inherited='598 full self-adjoint model;623 common domain/mass coordinates;717 bounded scalar feedback and normal packet;719 pair first derivative;742 restricted Gaussian-dilation theorem.',
    new='Full original finite-graph quantum H in strict Gauss preparations gives variance-driven fermion fourth cumulant and the same scalar-pair correlation/entanglement.',
    excluded='633 occupation instrument realization, autonomous terminal record, the730 continuum reference, all-scale quantum SM, dynamical quantum gravity.'),ensure_ascii=False,indent=2)+'\n')
write('round743_drafts/scope_and_dedup_review.md','''# 743范围与证明审查

上一轮完成742并执行743入口，属于进展。当前未发现活跃Python进程；无新代理或图像检验。入口冻结证据保持。

- 回查发现717已经给原逆度量的Majorana选择、sin s有界力及正常波包；719已给真空配对首项。本轮直接复用，不重新宣布这些结论。
- 条件128维参考不等于严格Gauss态。本轮另取598完整有限图的正规Gauss波包乘CAR真空，并明确这是不同准备，不是730连续真空。
- 原H在该初态上的一阶像：Hb保费米真空、Dirac和跳跃杀真空、Majorana创建原sterile对。后续全部作用保留，未假定有限配对子空间动力学闭合。
- 占据和双占据从t²开始，异常两点从t开始；代Wick给|Ys|²Var(x5)t²/ℏ²。正常波包方差严格正，不靠数值小数证明。
- 标量—配对关联及总费米线性熵共享同一方差；不同节点配对正交，熵系数为方差之和。
- 核态有H的任意有限矩，弱Taylor用于所列期望。不是无界算符全空间范数展开，也不给跨图统一余项。
- 旧η±相干可经同一原读口区分；这一不同输入家族不冒充真空已写出一个完整记录，不以相同均值替代完整instrument。
- 原H自动保Gauss和总能源；原质量／来源一致，但物理端点、连续和动态几何仍需核。
- 两组数值只核条件jet、实际CAR和正常配置积分，无完整H时间模拟。空间382—386、425、522—523和604、649／699保持。
''')
main=('research_note_743.md','joint_native_quantum_record.py','joint_native_quantum_record_results.json','unified_physics_condition_ledger_743.md')
write('round743_drafts/final_review.txt','Primary-agent review only. Native full-H Gauss preparation produces nonGaussian correlations; exact occupation instrument and autonomous record remain open.\n'+
    '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'743':'744','742':'743','741':'742','3422':'3424','3421':'3422',
         '1538':'1541','3387':'3401','3373':'3387','388':'389','293':'294',
         'joint_record_gaussian_boundary':'joint_native_quantum_record','gaussian_record_entry':'quantum_scalar_record_entry'}
publication=remap((HERE/'publish_round742.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第743轮完成：** [原量子标量与严格Gauss关联]({p}research_note_743.md)原完整H的正常Gauss准备产生非高斯四点缺陷；同一标量方差与Ys固定内部关联及纠缠，717读口／来源桥直接复用。无需为解除742限制增加物种，但理想仪器、末读自治及连续参考未签收。两组、十八式通过，最新743／3424，1541份编号科学文件、3401份保护证据。[核验]({p}research_round_743_checks.json)、[条件账]({p}unified_physics_condition_ledger_743.md)。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（743后，优先于下方历史安排）：** 接[744原作用诱导的实际仪器]({p}round744_drafts/STATUS.md)，核未知偶编码输入、原sin s末读、完整后态与共同来源；不只比较均值。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('原实际记录与二次实现的严格边界','原量子标量与严格Gauss记录关联')
publication=publication.replace('旧空间合同保持，平均来源不替代完整记录','旧空间合同保持，原相互作用与实际记录共同核验')
write('publish_round743.py',publication)
write('postcheck_round743.py',remap((HERE/'postcheck_round742.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round742.py').read_text('utf8'),mapping)
verification=verification.replace('==(1,0,0)','==(2,0,0)').replace('run=1,failures=0','run=2,failures=0')
verification=verification.replace("checks['display_formulas']==12","checks['display_formulas']==18")
verification=verification.replace('original_reference_record_partner_and_positive_Gaussian_gap_checked=True','original_full_H_Gauss_preparation_four_point_and_scalar_correlation_checked=True')
write('verify_round743.py',verification)
print('Prepared743 native Gauss process and744 actual instrument interface.')
