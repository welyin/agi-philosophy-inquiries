"""Prepare718 report evidence; keep all preceding material frozen."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name
    p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_717.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：718保留真实记录的有限时间过程与几何输送\n'+rest
ledger=ledger.replace(
    '接[716全账](unified_physics_condition_ledger_716.md)，回填[717报告](research_note_717.md)。[结果](joint_record_mass_feedback_results.json)、[核验](research_round_717_checks.json)。',
    '接[717全账](unified_physics_condition_ledger_717.md)，回填[718报告](research_note_718.md)。[结果](joint_record_history_transport_results.json)、[核验](research_round_718_checks.json)。')
marker='## 当前共同对象及仍存在的分支'
assert marker in ledger
ledger=ledger.replace(marker,'**718当前增量：** 同一局部资源给完整首次记录／后续记录／末态的有限时间移动界；原固定图有限能源与已知读口注能控制无条件前缀，不重置Gibbs。几何变化必须同步交叉端点、仪器与准备态；旧共同域和二阶来源直接复用。空间一致性及因果实现仍开放。\n\n'+marker)
updates={
    'C03 事件记录':'718首次结果、后续结果及末态共同受有限时间资源界控制，不除罕见分支概率',
    'C19 参考态':'718使用真实无条件历史前缀；原Gibbs只在读取前平稳，来源含准备响应',
    'C20 尺度映射':'718给实际局部势矩一致时的跨图条件界；固定图能源常数不提供共同空间极限',
    'C22 来源反作用':'718移动R的交叉端点及前两阶连接由同一几何确定，不能只变换H来源'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 718完整记录的有限时间与共同几何

- 原H和反射比较共用完整Gauss与623域；Duhamel加716资源控制实际记录位置的全cq误差。
- 多次原读口沿无条件前缀累积资源，首次及后续记录均保留，不除稀有记录概率。
- 598注能和形式下界保证固定图有限历史可用；图间常数、实际参考和共同极限另验。
- 几何来源共同输送R、H、Kraus及准备态，交叉项端点不可丢弃；二阶共同域与704实际来源复用。
- 三组条件原CAR校准通过，不作为全玻色或完整Gibbs数值模拟；717无界二阶反例保持。
- 一般迹距离和表示变换不是新定理，633因果限制不因有限时间连续性而消失。

## 本轮合并与下一项

C03／C19／C20／C22共用原模式、实际前缀态和几何路径，已经获得固定图有限时间控制。

接[719](round719_drafts/STATUS.md)：联合原区域操作与633因果限制，核同一记录能否进入局部过程；先复用旧结果排除仅靠表示输送的接法，工程设计后置。不再重复距离界或来源微分。
'''
write('unified_physics_condition_ledger_718.md',ledger)
write('round718_drafts/research_note_718_draft.md',(HERE/'research_note_718.md').read_text('utf8'))
write('round719_drafts/STATUS.md','''# 第719轮入口：真实记录、区域因果与共同过程

接[718](../research_note_718.md)、[全账](../unified_physics_condition_ledger_718.md)。实际首次记录、后续历史和几何来源已获得固定图有限时间合同；这个时间界不是传播光锥。

1. 先回查524／633因果读取、634总来源、667区域输送以及705—706区域Gauss合同，避免重做已知瞬时传信反例。
2. 核633原分布模式的因果限制是否在允许的共同参考／表示变换下保留；明确空间支持、测量时间窗与准确度，不把小时间连续性当作实现。
3. 如有可检验的有限精度必要条件，用同一原CAR读取及真实记录量化；认知工程后置，不凭空追加独立装置动力学。
4. 继续共同物理过程目标，保留382—386、425、522—524空间接口以及699限定反例；不重开高矩、一般反射、Taylor余项或谱精度扫描。
''')
write('round718_drafts/literature_scope_audit.json',json.dumps(dict(
    source='https://arxiv.org/abs/quant-ph/9712042',
    authors='Christopher A. Fuchs; Jeroen van de Graaf',
    read='Primary arXiv abstract and metadata only; related distinguishability background.',
    proof='Special pure-state distance bound and actual instrument isometry identity proved directly in note.',
    inherited='598 energy; 623-625 domains and source histories; 704 preparation; 716 resource; 718 entry cross term.',
    not_imported='No general all-states Taylor, uniform spatial-limit or causal-instrument theorem.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round718_drafts/scope_and_dedup_review.md','''# 718范围与去重审查

新增为原完整H、首次真实记录、后续历史、同一局部资源与几何来源的连接。Duhamel、纯化距离、反射身份及623—625／704来源框架均不重算一般新定理。

[U,R]=R(Uprime-U)确定式(6)常数，不要求初态与R对易。多历史滑动使用尚未插入首次读取的无条件前缀，保全部结果；不把分支条件态作统一结论。

有限形式能源通过598下界和已知原L注能控制前缀。固定图常数不能当空间一致界。Gibbs只作读取前准备，不在记录后重置。

几何为实验参数，等待期间固定。R及前两阶导数与参考T和W对易；RdotR为连接。二阶式中2[Gamma,RGR]和Gamma-dot均保留。交叉端点Rdot与热准备导数独立不可省。

数值先生成原64模式系数，再取声明配置下精确不变的八中性模；不是完整玻色历史，也不是全Gauss Gibbs。多历史积分为浮点诊断，非区间证书。

旧空间384消去的Lipschitz不恢复；386／425替代桥保持。633有限时因果装置与实际共同空间极限仍缺。无独立代理复核。
''')
names=('research_note_718.md','joint_record_history_transport.py',
       'joint_record_history_transport_results.json','unified_physics_condition_ledger_718.md')
write('round718_drafts/final_review.txt','Main-agent review only. Full-H domains, CQ constants, '
      'unconditional prefixes, actual energy injection, moving-source endpoints and fixture scope checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared718 evidence and719 interface.')
