"""Prepare717 audited mass-feedback evidence and finite-time interface."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_716.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：717原记录的质量反馈与有界物理读数\n'+rest
ledger=ledger.replace(
    '接[715全账](unified_physics_condition_ledger_715.md)，回填[716报告](research_note_716.md)。[结果](joint_smooth_mode_contract_results.json)、[核验](research_round_716_checks.json)。',
    '接[716全账](unified_physics_condition_ledger_716.md)，回填[717报告](research_note_717.md)。[结果](joint_record_mass_feedback_results.json)、[核验](research_round_717_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**717当前增量：** 真实CAR读取通过原质量力改变后续玻色二阶响应。正常Gauss态族表明有界能源／读取矩不保全部无界势响应；原sin s效果却有由同一Majorana、曲目标和体积固定的全配置二阶界。有限时间共同历史和物理连续仍开放，停止高矩反例扫描。\n\n'+marker)
updates={
 'C03 事件记录':'717原sin s效果对sterile读取的二阶概率反馈有态无关精确上界，尚无统一时长余项',
 'C15 物种手征':'717已有Majorana决定原singlet读口二阶反馈；非新增物种或耦合',
 'C19 参考态':'717正常非热Gauss态反例区分有界能源与无界势反馈，未否定原Gibbs',
 'C20 尺度映射':'717物理报告权q/w有界可统一旧有界效果的二阶系数，不等于连续历史存在',
 'C22 来源反作用':'717同一原质量力固定读后标量响应和均匀体积导数，无独立反馈强度'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 717真实记录的质量反馈与有界读数

- 原CAR读取保持整个即时玻色边缘及配置量的一阶时间变化，二阶差却由原曲目标梯度收缩质量力决定。
- 600已有质量力、623已有目标多项式、669已有增长界、712已有交换子工具。本轮只计实际仪器、原参考族及原效果之间的新连接。
- 同一固定图的正常Gauss混合态保持有界总能源、P_f和D_f四阶矩，P_f的二阶反馈仍可无界；没有交换状态族与小时间极限，不称有限时爆炸。
- 原sin s效果的Dirac力消去，Majorana力给精确全配置范数。物理报告权q/w控制跨图二阶界，不再逐个追加高矩反例。
- 同一均匀体积变化给来源因子-6；非均匀模式／权重变化沿716另加。动态量子度规及完整有限时因果实现仍未证明。
- 四组校准含原96模式、正常H5波包、物理求积及原效果；完整相互作用Gibbs未数值模拟。旧空间及699范围不变。

## 本轮合并与下一项

C03／C15／C19／C20／C22共用原质量力、真实读后态及体积。无界势响应失败与有界物理读数可用应分开。

接[718](round718_drafts/STATUS.md)：核同一读取的有限时间比较表示H+2D=RHR，将原有界读数与623—625／704实际历史、634来源守恒一同接入。停止高矩反例及波包优化。
'''
write('unified_physics_condition_ledger_717.md',ledger)
write('round717_drafts/research_note_717_draft.md',(HERE/'research_note_717.md').read_text('utf8'))
write('round718_drafts/STATUS.md','''# 第718轮入口：原有界记录的有限时间共同历史

接[717](../research_note_717.md)、[全账](../unified_physics_condition_ledger_717.md)。同一质量力把原sterile读取与后续玻色读数接通；能源界不保全部无界势响应，旧sin s效果却有态无关二阶系数界。

1. 原效果A与R_f对易，真实读后概率差可精确写成同一初态下H与R_f H R_f=H+2D_f的演化比较。这个比较不是另加一条实际自然规律；先核原Gauss、域、质量／跳跃和几何来源均随同输送。
2. 回查623—625、634、704的真实历史和来源导数，区分小来源二阶响应与有限改变。只补现有比较表示真正缺的有限时间／共同状态接口。
3. 不再制造高矩反例、优化波包或重复一般Zeno。已有二阶界不自动给统一有限时余项，更不签收空间连续。
4. 保留原完整H、实际初态和非选择后态，不在读取后重热化；若引入新的控制或因果探针，须逐项交代其内部耦合和总来源，工程设计后置。

旧空间382—386、425、522—524、633因果实现、634来源拼接及699固定候选限制继续保留。目标是四分支在同一真实过程中的联合连接。
''')
write('round717_drafts/literature_scope_audit.json',json.dumps(dict(
    source='https://arxiv.org/pdf/1003.3372',authors='Gero Friesecke and Bernd Schmidt',
    read='Theorem 1 and domain discussion, PDF pages1-3; HTML unavailable, primary PDF read as text.',
    scope='Domain invariance and self-adjoint observables yield first-order expectation differentiability.',
    own_extension='Second-order statements restricted to compact smooth Gauss core using the inherited graph bounds.',
    inherited='600 mass force; 623 domains and exact potential; 669 growth; 712 metric contraction.',
    not_imported='A general all-states second-order Ehrenfest theorem or a classical/hybrid replacement.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round717_drafts/scope_and_dedup_review.md','''# 717范围和推导审查

原质量力不是新发现：600已给；669和712的成熟工具复用。新增为实际CAR读取后非平稳态、原势响应及旧sin s效果的共同反馈系数。

双交换子符号为delta A''=-Tr rho delta_f F_A，正反例用原Majorana相位与完整CAR校准。K^-1全部交叉项保留，否则Dirac消去会错。

反例用H5正常Gauss包和偶sterile相干，混合权R^-4。整个图额外边势O(log R)^2，质量O(R)，并未从点态大系数冒充正常物理态。不能交换R和t极限，不称有限时爆炸。

原有界效果的精确范数含1/4和q/w；φ=0的Gauss偶态逼近给上确界。只给二阶系数及核态Taylor，不声称统一时长、全态Taylor或动态GR。

草稿三组结果及当前四组结果均保留。范围布尔元数据在正式收据前修正。四组代码、数值、报告表格及公式逐项核对，无独立代理复核。

旧空间384删去Lipschitz保持；386与425替代；523共同实现和524局域UV不重证。完整物理端点及四分支同一连续过程仍开放。
''')
names=('research_note_717.md','joint_record_mass_feedback.py',
       'joint_record_mass_feedback_results.json','unified_physics_condition_ledger_717.md')
write('round717_drafts/final_review.txt',
      'Main-agent review only. Actual time-domain scope, full mass force sign, normalized H5 packets, '
      'total energy growth, original bounded-effect norm and geometry weights checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared717 evidence and718 interface.')
