"""Prepare722 evidence, retaining every preceding frozen file."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_721.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：722原平方参考的能源障碍与可容许读取\n'+rest
ledger=ledger.replace('接[720全账](unified_physics_condition_ledger_720.md)，回填[721报告](research_note_721.md)。[结果](joint_velocity_record_energy_results.json)、[核验](research_round_721_checks.json)。',
                      '接[721全账](unified_physics_condition_ledger_721.md)，回填[722报告](research_note_722.md)。[结果](joint_squared_reference_readout_results.json)、[核验](research_round_722_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**722当前增量：** 平方相位instrument可把原正常Gauss有限能源态送出能源域；该失败限指定读法。完整Q_f含原W_f的有理instrument却保原完整固定图能源，且几何一阶Kraus导数在能源域有界。保单参考谱分布与全部记录，不保证四参考联合坐标、自主实现或空间细化。\n\n'+marker)
updates={
    'C03 事件记录':'722完整Q_f的实际有理instrument保全部结果和原等待的有限能源历史；指定平方相位读法有反例',
    'C09 参考':'722单变量Q_f的谱信息可用保资源读口采集；不是不对易四参考的联合锐读或坐标逆',
    'C19 参考态':'722存在正常Gauss有限能源态证明一种有界读口仍会产生无限原能源；替代读口覆盖任意正常有限能源准备',
    'C22 来源反作用':'722有理读口在同一能源域中可微，立即读口的来源、仪器、热准备一阶总响应共同有定义'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
# C09 label follows the historical ledger's exact wording.
if 'C09 参考' not in seen:
    for i,row in enumerate(rows):
        if row.startswith('|C09 '):
            parts=row.split('|');parts[2]+='；'+updates['C09 参考'];rows[i]='|'.join(parts);seen.append('C09 参考')
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 722原平方参考、真实读取与共同能源

- 652自伴性、平方压缩和721速度域结论直接复用。新增反例在原五维测度、原Gauss态和完整H形式上成立，不是截断矩阵反例。
- 局部四维Hardy控制与原T流坐标使平方自由尾部导致无限原能源；只排除指定相位读法。
- 对完整W_f−V_f²，原速度能源群界与有界乘法给预解式能源界；指定有理Kraus保Gauss、宇称、原谱分布及有限能源。
- 空间差分不删去，不要求与速度对易。有理相位属于实际后态，不能换成Luders而沿用证明。
- 几何紧集内，Kraus及一阶导数共同作用于原能源域；与704热准备导数接成立即读取的完整三来源响应。
- 全部结论固定有限图、外部正几何，不覆盖量子θ、自主装置、空间一致或四坐标联合重建。

## 本轮合并与下一项

C03／C09／C19／C22共同使用原Q_f、实际Kraus、完整能源和来源。停止谱函数选择及常数优化，下一项核整个参考任务，而非继续证明更多单变量读口。

接[723](round723_drafts/STATUS.md)：回查554、647—652、704、707—708与旧空间合同。区分同一参考的完整态、有限任务、联合记录与坐标重建；已有定理直接推论不计新轮。旧699与全目标保持。
'''
write('unified_physics_condition_ledger_722.md',ledger)
write('round722_drafts/research_note_722_draft.md',(HERE/'research_note_722.md').read_text('utf8'))
write('round723_drafts/STATUS.md','''# 第723轮入口：同一完整量子参考与关系任务

接[722](../research_note_722.md)、[全账](../unified_physics_condition_ledger_722.md)。原导数平方参考已有保资源读取，不等于整个关系坐标任务完成。

1. 回查554、647—652、704、707—708和382—386／425／522—524。明确旧条件性坐标、热参考、实际方向仪器分别需要哪些共同对象，避免重新加入384已消去的Lipschitz条件，386与425保留替代关系。
2. 保原T、s、Q_T、Q_s、共同状态、Gauss区域、实际记录及来源。有限顺序仪器存在是成熟推论；不把它当无扰动联合锐读或坐标逆。单变量分布信息不等于完整量子参考态。
3. 核同一任务沿625／704共同近似及707／708尺度映射的条件；完整输出态、记录、原等待和资源必须共同保留。能直接用旧结论完成的只记入口，不虚增轮次。
4. 停止更多谱函数、参数、尾部及常数优化。自主工程后置，先检验实际对象映射或真正不相容条件。
5. 固定图、辅助手征、指定连续物质及经典几何分支保持区分，699失败不扩大。无新目标、任务或图像检验。
''')
write('round722_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(url='https://arxiv.org/pdf/0803.0503',authors='Frank and Seiringer',
                  read='Introduction equation(1.1), page1, classical local Hardy inequality N=4,p=2.',
                  scope='Mature inequality only. Fractional main theorem not imported; local square-completion proof included.')],
    inherited='652 self-adjoint original references;721 energy-space velocity group;623 full energy-form comparison;704 energy-weighted preparation derivative;722 entry exact scalar tail.',
    own_mapping='Five-dimensional original Gauss lift and local energy obstruction; complete squared-reference rational instrument with noncommuting spatial term; energy-space source derivative.',
    Stieltjes_argument='Elementary spectral transform and analytic identity/inversion, not a stable finite-data tomography theorem.',
    no_image_checks=True),ensure_ascii=False,indent=2)+'\n')
write('round722_drafts/scope_and_dedup_review.md','''# 722证明、重复与范围审查

原五维体积不能漏Higgs的R³因子；m=MR⁴/(wD^(3/2))已由独立Jacobian及两种积分检查。R=0是内部零测度集合，不是删除构形边界。空CAR真空为Gauss不变态，但不要求原H使它不变。H¹准备足以检验有限能源域，不冒称H²准备。

入口Fourier余项直接复用；局部四维Hardy与原完整H形式等价共同给反例。两个平方流方向都出能源域，实际两Kraus至少一个非零分支无限，不能凭平均公式声称两分支分别无限。反例只涉及指定读法。

正结果保完整空间W；其乘法及梯度在固定图有界，故在能源空间有界。预解式先在能源空间由群积分和Neumann得到，再用原物理自伴逆唯一性识别。没有交换非对易D与V，也没有凭有限矩阵有界证明全域。

K0与K1有实际相位，非Luders。保原谱分布不保全部相干。所有响应参数的理想概率确定单变量测度，不等于有限样本稳定重建。大参数是明确读口选择，不是新物理常数或配置阈值。

几何微分使用共同原D(V²)、紧集统一能源域和预解式恒等式；R′Z的能源有界性关闭立即读口总一阶导数。含H的项按形式解释。含等待的几何总导数仍沿704原合同，不因有限历史预算自动推广。

64维非零W采用条件化邻居值，仅校准代数与来源，未宣称完整双节点Gauss实现或完整热谱。解析全图保全部CAR、规范和两类质量，原图与正外部几何仍是输入。量子θ、连续尺度、自主装置及整体量子坐标仍开放。
''')
names=('research_note_722.md','joint_squared_reference_readout.py',
       'joint_squared_reference_readout_results.json','unified_physics_condition_ledger_722.md')
write('round722_drafts/final_review.txt','Main-agent review only. Original Gauss lift, local Hardy form, actual phase-instrument failure, complete noncommuting squared-reference resolvent instrument and energy-space first source derivative reviewed.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared722 evidence and723 whole-reference interface.')
