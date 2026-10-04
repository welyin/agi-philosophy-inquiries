"""Prepare720 evidence; inherited protected drafts are never overwritten."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_719.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：720传播中的内部参考、原质量与实际记录\n'+rest
ledger=ledger.replace('2026-10-03。接[718全账](unified_physics_condition_ledger_718.md)，回填[719报告](research_note_719.md)。[结果](joint_local_relational_record_results.json)、[核验](research_round_719_checks.json)。',
                      '2026-10-04。接[719全账](unified_physics_condition_ledger_719.md)，回填[720报告](research_note_720.md)。[结果](joint_momentum_spin_reference_results.json)、[核验](research_round_720_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**720当前增量：** 原两质量允许均匀自由分支的FW守恒自旋，但原位置、记录与来源必须共同输送；对所有原径向非均匀质量要求同一平移不变参考时，共同交换子仅为两质量标量，不能承载SU2。这是明确量词下的参考接口分类，不是否定所有联合参考或统一目标。\n\n'+marker)
updates={
'C03 事件记录':'720原CAR记录须沿动量相关FW同构完整输送；冻结旧读口会改变后态',
'C06 Lorentz':'720原均匀有质量Weyl分支有守恒平均spin，原空间局部性不自动继承',
'C15 物种手征':'720保原Dirac和Majorana，径向原场同时缩放全部质量而破坏固定FW保护',
'C19 参考态':'720固定背景自由参考可输送；同一背景无关平移不变SU2不能适应全部原径向质量',
'C22 来源反作用':'720FW几何连接和真实仪器导数同时必需，不当新引力来源'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 720参考、传播、原质量与几何的共同字典

- 原633的两个正且不等质量各给四维Clifford块；同一准确FW变换保原Nambu实结构。
- 自由守恒平均spin存在，不能把720入口的固定字典下界扩大为一般无守恒spin。
- 平均spin的动量乘子非恒定，因此不保持原位置点局部动作；区域代数若变换须明示。
- 原占据数完整instrument及后态须同步输送。实际两动量记录另用256维Fock核对能源。
- 几何来源包含表示连接，真实仪器参数导数不可略去；这些是已有导数工具的实际原质量映射。
- 原径向标量c(x)同时缩放Dirac／Majorana。若同一有界平移不变乘子保全部正背景演化，它只能是两质量投影的标量组合。
- 分类仅限制背景无关的内部乘子，不排除联合场、轨道、随背景变化及有限精度任务。常背景自由结论不是完整相互作用连续证明。

## 本轮合并与下一项

C03／C06／C15／C19／C22须沿同一观测和参考字典匹配。停止继续保护spin或质量分支扫描，转回原物质参考与实际方向记录。

接[721](round721_drafts/STATUS.md)：复用548—554、573—574、647—651及523／524，核原参考图册如何接同一量子记录；不得重做已完成的秩、角点、坐标变分或尾部反例。工程设计后置，旧空间及699范围保持。
'''
write('unified_physics_condition_ledger_720.md',ledger)
write('round720_drafts/research_note_720_draft.md',(HERE/'research_note_720.md').read_text('utf8'))
write('round721_drafts/STATUS.md','''# 第721轮入口：原物质坐标与实际方向记录的共同来源

接[720](../research_note_720.md)、[全账](../unified_physics_condition_ledger_720.md)。均匀自由FW参考可用，但不是所有非均匀物质变化下的固定内部spin；停止保护参考分类和精度扫描。

1. 复用548—554、573—574、647—651已给的h、s、导数参考、共同受约束经典态、关系定位及边界来源。不能再把零阶秩、实际满秩例、定位变化或平方读口高矩反例计作新轮次。
2. 复用382—386、425、522—524的条件性空间接口和实际仪器。384已删Lipschitz不恢复，386／425替代，523共同实现及524给定背景局部标量探针保持。
3. 共同核查同一原物质态能否给空间方向菜单及原CAR记录；首先找到尚缺对象映射，不另设免费独立参考或仅以经典图册代替量子instrument。
4. 可使用707—708已有共同来源／保留资料、716—718已有实际记录资源和704—706来源输送。工程设计后置，不再靠设计静止子系统回避共同传播。
5. 先做去重入口；若只是旧条件的重述，记入口审计而不增加正式轮次。原699、固定图正H、辅助与自由连续各分支范围不变。
''')
write('round720_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://journals.aps.org/pr/abstract/10.1103/PhysRev.78.29',
                  authors='Foldy and Wouthuysen',
                  read='Original publisher abstract: free mean-spin conservation and distinct position/spin operators.',
                  scope='Mature FW interpretation. Matrix formula for the actual two original masses proved directly here.')],
    failed_fetch='MDPI exact-FW page returned HTTP429; no claims rely on that page.',
    inherited='374/381 topology, 375/663 frame transport, 606 projection, 718 parameter connection.',
    own_mapping='Original two masses, Nambu reality, true CAR record, geometric response and uniform-in-profile commutant classification.',
    no_image_checks=True),ensure_ascii=False,indent=2)+'\n')
write('round720_drafts/scope_and_dedup_review.md','''# 720范围与推导审查

准确FW公式为成熟工具，本轮不申报新的自由Dirac理论。两个质量均来自633，K与质量反对易，正质量是本构造的必要范围。全CAR实结构已核，未默认无限体积Fock酉可实施性。

真实占据数记录的非选择两点函数在同一U下输送；混合后态不假定Gaussian。用八模式256维Fock独立验证Nambu半因子及加能，保真实两投影。

原几何形变时的来源与instrument连接直接继承718工具，数值给原质量／记录的非零遗漏量；不是额外引力场。均匀双动量权重时能源仪器项抵消，正文保此例外。

局部性结论只针对原位置净；变换整个净不产生物理超光速操作。不能从非局部乘子直接给最短实现时限。

全正径向背景共同分类先用常数差，再用cos和sin平移；B0可逆使乘子恒定，两个四维Clifford全矩阵块给标量交换子。数值零空间仅校准，非解析替代。量词是同一个背景无关平移不变乘子，不包括全相互作用动态场或场依赖／轨道参考。

初次代码旋转共轭漏转置已修，最终三组及结果逐项重算。旧冻结文件不变。无独立代理审查、无图像检查。
''')
names=('research_note_720.md','joint_momentum_spin_reference.py',
       'joint_momentum_spin_reference_results.json','unified_physics_condition_ledger_720.md')
write('round720_drafts/final_review.txt','Main-agent review only. Positive mass, FW algebra, CAR reality, exact Fock record, locality and profile quantifiers checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared720 evidence and721 interface.')

