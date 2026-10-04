"""Preserve716 evidence and prepare the original nonequilibrium interface."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_715.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：716光滑物质模的共同资源与来源\n'+rest
ledger=ledger.replace(
    '接[714全账](unified_physics_condition_ledger_714.md)，回填[715报告](research_note_715.md)。[结果](joint_sharp_record_domain_results.json)、[核验](research_round_715_checks.json)。',
    '接[715全账](unified_physics_condition_ledger_715.md)，回填[716报告](research_note_716.md)。[结果](joint_smooth_mode_contract_results.json)、[核验](research_round_716_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**716当前增量：** 原sterile波包的同一局部势矩与跳跃作用共同控制读取差量、时间变化及指定几何导数；物理体积归一给无显式节点数的条件界。非均匀几何还改变模式，来源不能漏仪器导数。完整参考的跨图矩、实际传播及因果实现仍未证明。\n\n'+marker)
updates={
 'C03 事件记录':'716原物质模读取差量与真实记录时间变化由同一局部势矩和传播块控制',
 'C19 参考态':'716原完整参考需控制波包加权势矩，固定图热矩存在不等于跨图一致',
 'C20 尺度映射':'716体积归一、同一波包和原跳跃共同进入条件界；二维三维都适用，未选维数',
 'C22 来源反作用':'716固定物理波包时有几何仪器导数，须与原来源及热准备响应同时保留'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 716物理波包、原参考和来源的同一合同

- 保原全部32CAR／节点、动态标量、完整Gauss与原跳跃。读取只选规范平凡sterile模，没有删带荷物质。
- 598旧束缚势使读取差量范数受κ+C(AP_f+B0)^(1/4)控制。相同P_f控制其四阶矩和记录均方时间变化；这不是传播速度证明。
- 给定物理体积的模式归一产生r_f=max|ϕ|²／Σw|ϕ|²，局部参考矩可避免要求扩张系统有有限总势能。
- 背景非均匀变化时，模式导数由体积变化的加权方差固定。原来源、仪器改变和热准备协方差均需纳入。
- 波包、体积形式、实际原参考和跳跃仍待共同匹配；没有完成连续相互作用物质、统一热矩或自主因果装置。
- 数值为原96模式点系数及物理求积，不是完整热态或动态量子几何模拟。旧空间及699范围不变。

## 本轮合并与下一项

C03／C10／C19／C20／C22共用同一模式、局部势矩与原质量／跳跃；不能分别拟合记录、时间和来源。

接[717](round717_drafts/STATUS.md)：核实际原读取后非平稳态及有限时间过程，复用623—625、632—637及667；不继续扫描波包或范数，不把有限图域保持当统一物理连续界。
'''
write('unified_physics_condition_ledger_716.md',ledger)
write('round716_drafts/research_note_716_draft.md',(HERE/'research_note_716.md').read_text('utf8'))
write('round717_drafts/STATUS.md','''# 第717轮入口：原记录后的参考、传播与共同资源

接[716](../research_note_716.md)、[全账](../unified_physics_condition_ledger_716.md)。同一物理波包的局部势矩、原跳跃和体积匹配共同控制读取差量、时间变化及指定几何导数；实际跨图参考及传播还未实现。

1. 优先检验实际读取后非平稳态。原CAR仪器与全部玻色乘法量对易，立即时的局部势矩应原样保留；但其后真实H演化不能假定仍为原Gibbs。先核此对象身份，再核有限时间控制。
2. 复用623—625原多时响应、632—637来源守恒／区域边界、667已有过程运输。固定图域结论不重复计为新进展，不把非选择记录后态重新热化。
3. 联合条件仍为实际状态、读取、原质量／跳跃、几何来源及尺度，不退回自由波包精度优化。若需额外传播或参考假设，说明同时补哪个部门及排除什么模型。
4. 若转向局域探针，只在核实Fewster—Verch原定理前提及与原场的实际映射后使用；有限时内部耦合、资源和总来源先行，具体认知硬件设计后置。

633瞬时分布式测量限制、634来源拼接限制及699固定候选反例保留。旧空间382—386、425、522—524直接复用，不恢复已删假设。
''')
write('round716_drafts/literature_scope_audit.json',json.dumps(dict(
    source='https://arxiv.org/abs/1810.06512',authors='Christopher J. Fewster and Rainer Verch',
    read='Abstract and bibliographic record; no uninspected realization theorem imported.',
    scope='Local system-probe coupling and scattering induce operations; causal factorization is a separate condition.',
    own_proof='Original CAR cross-block algebra, Jensen bound, volume derivative and inherited Kubo identity.',
    inherited='598 confinement; 623 domains; 633 free smooth mode and causality limitation; 714 thermal response.',
    not_claimed='Continuum existence, full Gibbs numerical integral, dimension selection, autonomous instrument.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round716_drafts/scope_and_dedup_review.md','''# 716范围与去重审查

入口三模恒等式、598束缚势、712一般热谱矩、633连续自由场结果均复用。本轮增量是实际波包的局部势矩同时控制读取、原时间变化及几何响应。

逐纤维矩阵界与标量P_f对易，允许正常相关态的二阶及四阶矩估计。全H仍含动态玻色和全部带荷CAR；数值只测其点系数，不冒充完整热态。

f配置无关；w是背景几何而非动态F。否则玻色对易性失效。质量在声明Einstein变量／固定φ／标架下不随空间体积变；更一般变化需补质量、节点及标架导数。

几何遗漏项的实际残差非零，原热准备协方差另存。统一矩只给充分预算和尾控制，不能产生传播速度或证明共同连续极限。

旧空间定理保持原准确版本：384已删Lipschitz；386与425替代；523共同实现及524局域UV不重证。完整物理对象身份仍开放。
''')
names=('research_note_716.md','joint_smooth_mode_contract.py',
       'joint_smooth_mode_contract_results.json','unified_physics_condition_ledger_716.md')
write('round716_drafts/final_review.txt',
      'Main-agent review; no independent review. Full CAR signs, mass norms, local potential moments, '
      'stationary time bound, geometry normalization derivative and Kubo scope checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared716 report evidence, ledger and717 interface.')
