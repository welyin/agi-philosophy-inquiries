"""Prepare721 report evidence and next shared-reference audit."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_720.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：721原物质速度读口的完整能源与互补质量反馈\n'+rest
ledger=ledger.replace('接[719全账](unified_physics_condition_ledger_719.md)，回填[720报告](research_note_720.md)。[结果](joint_momentum_spin_reference_results.json)、[核验](research_round_720_checks.json)。',
                      '接[720全账](unified_physics_condition_ledger_720.md)，回填[721报告](research_note_721.md)。[结果](joint_velocity_record_energy_results.json)、[核验](research_round_721_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**721当前增量：** 原T／s速度流的全域微分界及原束缚势使两结果instrument保完整固定图能源形式域，保所有记录和原等待；原流对Dirac／Majorana有互补反馈。固定instrument来源有定义，变化instrument的总导数仍分域说明。不含量子θ、跨图一致或自主装置。\n\n'+marker)
updates={
'C03 事件记录':'721原速度instrument保持全部结果、形式域及原等待的有限历史；不是四坐标联合读取',
'C15 物种手征':'721同一原曲目标使s速度流保Dirac、T速度流保Majorana；其它反馈由原流固定',
'C19 参考态':'721任意正常有限能源态读后仍有限能源，无须重置Gibbs；不推出参考不受扰动',
'C22 来源反作用':'721固定instrument的来源形式受完整能源控制；固定t时须再保仪器导数，普通有限能源不足以自动保证全部能源导数'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 721原物质速度读口的完整资源

- 652的参考自伴性及共同表示直接复用；此次补原未截断流对旧能源域的保持。
- 全域Hessian、散度梯度及logF变化控制半密度动能与原势；598／623形式范数等价带回全部规范、质量和边项。
- 两结果Kraus保Gauss、宇称和有限能源。完整结果历史共用原H等待，不重置态、不丢罕见记录。
- 原嵌入坐标中Xs只动z、XT只动xH，因此原Dirac和Majorana反馈互补，不另设反馈系数。
- 原几何固定instrument来源有有限期望；固定t的参数变化另有Kraus导数，不能把能源期望存在当全部能源总导数存在。
- 结论限固定图及外部正几何，不覆盖额外量子θ或跨图统一常数。合法仪器不等于自主或连续局域装置。

## 本轮合并与下一项

C03／C15／C19／C22的参考、读取、能源和质量共用原过程。下一项转回原导数平方菜单及完整关系任务，停止更多速度谱函数或保守常数优化。

接[722](round722_drafts/STATUS.md)：回查554、652、704—708；核Q_f与真实记录、原来源及方向接口的同一性。只复述旧结果时不计新轮次；保留旧空间和699准确范围。
'''
write('unified_physics_condition_ledger_721.md',ledger)
write('round721_drafts/research_note_721_draft.md',(HERE/'research_note_721.md').read_text('utf8'))
write('round722_drafts/STATUS.md','''# 第722轮入口：导数平方参考与实际关系任务

接[721](../research_note_721.md)、[全账](../unified_physics_condition_ledger_721.md)。原速度读口在固定图保完整能源，但不是四个量子坐标的联合读取。

1. 先回查554、652、704—708。652已有Q_f=W_f−V_f²的自伴／闭形式、同一有限表示及平方尾项，禁止重复“先平方再压缩”的旧反例。
2. 核真正需要的任务：完整参考态、有限精度关系记录或原方向效果。参考可观测量、自伴性、实际instrument与坐标重建分开；不能把另外设计的效果当原任务。
3. 721已保速度读口资源，停止函数、参数和保守界优化。若其与Q_f的连接可独立推进，必须同时保原空间差分项、完整质量、状态和来源；否则先记录可检验的缺项。
4. 旧空间382—386、425、522—524及原647—651直接复用，不恢复已消去假设；固定图、辅助和自由连续的分支保持分开。
5. 工程设计后置，无新应用任务、图像检验或目标修改。入口整理不计新科学轮次。
''')
write('round721_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://arxiv.org/pdf/2410.02383v2',authors='Beauchard and Pozzoli',
                  read='Section5 Definition24 and Lemma25, half-density transport and Jacobian.',
                  scope='Mature representation only; no controllability theorem imported as autonomous implementation.')],
    inherited='588 compactly modified flows;652 original self-adjoint flows;598/623 form equivalence;704 preparation response;712 geometric cancellation.',
    own_mapping='Global old-field form bounds, full original energy preservation, actual recorded histories and complementary mass response.',
    no_image_checks=True),ensure_ascii=False,indent=2)+'\n')
write('round721_drafts/scope_and_dedup_review.md','''# 721范围及证明复核

本轮新增是未加配置截断的原流保原完整形式域，不重证652的本质自伴和588酉表示。沿两个时间方向的界给域等号。所有域常数固定图，θ若量子化则另有导数，当前未覆盖。

原H与T+W移位形式等价包含全部规范、边和CAR项，因此不是用纯玻色H替换原动力学。全Fock矩阵质量无界增长由F的比值及旧束缚势共同控制。数值只核几何、32模式系数和旧径向64维诊断。

实际Kraus保结果及后态；非选择随机酉表达仅用于能源求和，不能改称真实控制器。原H等待保移位能量，有限历史的乘积预算不保证无限序列或空间连续。

Xs xH=0、XT z=0是原曲目标的准确身份，对应两种质量的互补反馈。开发中最初要求全部质量反馈都非零被此身份纠正；最终断言反而核精确保留块，没有删原项。

固定instrument的来源是能源有界形式；固定t的仪器导数只在说明的域使用，未用有限平均H推出全部无界能源导数。固定t/w是另一个明确合同。无自主装置、普遍三维或引力生成结论。
''')
names=('research_note_721.md','joint_velocity_record_energy.py',
       'joint_velocity_record_energy_results.json','unified_physics_condition_ledger_721.md')
write('round721_drafts/final_review.txt','Main-agent review only. Original noncompact flow, half-density derivative, full energy form, recorded history, complementary mass and source domains checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared721 evidence and722 interface.')

