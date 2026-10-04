"""Preserve round710's explicit extension, scope audit and continuation."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_709.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：710实际荷变化、共同来源与拓扑相位边界\n'+rest
ledger=ledger.replace('2026-10-02。接[708全账](unified_physics_condition_ledger_708.md)，回填[709报告](research_note_709.md)。[结果](joint_charge_quantum_coarse_results.json)、[核验](research_round_709_checks.json)。',
    '2026-10-03。接[709全账](unified_physics_condition_ledger_709.md)，回填[710报告](research_note_710.md)。[结果](joint_charge_changing_vertex_results.json)、[核验](research_round_710_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**710当前增量：** 在原CAR／Gauss对象上新增指定局部QQQL乘积，实际改变总荷而保原跨区跳跃、Gauss及立即读口预算；有新的完整热参考。残余粗化保这些更新，全U(1)粗化不保。荷源、相位和体积响应必须共同匹配；只给系数命名为拓扑角仍可在该有限扩展中消去，不构成反常生成。新增系数、味张量和EFT范围明确记为输入，原H与699失败范围保留。\n\n'+marker)
updates={
    'C01 量子对象':'710以原CAR不变量给实际有界扩展；新残余条件期望保新增荷更新',
    'C02 区域组合':'710局部荷源与原Wilson边流同时存在；保原Gauss切口配对，不逐区套全局中心商',
    'C03 事件记录':'710新增项与原标量读口对易，立即能量预算算子相同，但使用新态且概率可变',
    'C04 内部演化':'710新增有界B给同一域自伴H1；原H0仍单列，H1的总荷有真实非零导数',
    'C16 反常测度':'710允许的QQQL项不等于瞬子生成；单系数相位在有限扩展可消去，独立拓扑贡献仍缺',
    'C17 质量机制':'710保全部原Dirac／Majorana质量；新增项次数4n_g，不以其荷规律产生原质量参数',
    'C18 参数':'710有拓扑菜单时新增相位给eta-theta_w不变量；一代受限参数族多一个圆，系数未导出',
    'C19 参考态':'710有界扰动给新完整Gauss Gibbs及热迹比较；不宣称原Gibbs原样保持',
    'C20 尺度映射':'710沿664单元归一得到w^(1-2n_g)；固定连续系数不保统一裸UV界，需实际匹配',
    'C21 主体存储':'710保更新合同排除删荷相干的旧粗化；参数常量无外部驱动，总能以新H计',
    'C22 来源反作用':'710同一顶点强制荷源、相位和体积混合源身份；未给量子度规或完整动态高阶域'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            c=row.split('|');c[2]+='；'+val;rows[i]='|'.join(c);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 710新增输入与同时减少的自由选择

- 保留原H0，在相同CAR／Gauss对象上声明H1=H0+B。QQQL代表及每代一因子的乘积、幅度和相位是额外模型输入，不是认知推导或已计算的瞬子顶点。
- 非零局部不变量有36个单项式、原真空范数平方72；新增项有界。因此原无界玻色及规范动力学仍在，新H和新Gibbs有同一严格存在证明。
- 一个残余荷条件期望同时保新H、新热态及中性原记录；全U(1)夹断会删真实荷源。原区域流与局部荷源共同存在，不能逐区域删除参考相干。
- 连续拓扑菜单与显式QQQL相位共同给eta-theta_w；原629三角结论只属于原无显式重子破坏菜单。有限新H只有一个共同相位时，该相位仍可通过全局重相位消去。新增荷变化不是拓扑机制的充分证据。
- 同一单元归一给w^(1-2n_g)，固定连续系数时的几何源为6(1-2n_g)B，混合相位源同步。不能另外独立拟合体积响应。固定有限图成立不保证UV极限。

## 本轮合并与下一项

C01／C03／C04／C19在新增作用下共同实现，C02保原区域配对，C16／C18严格区分系数相位与独立拓扑贡献，C20／C22给实际来源匹配。没有选出群、维数、代数目、作用强度或量子GR。

接[711](round711_drafts/STATUS.md)：回查原规范构形和629拓扑部门，先核独立拓扑资料是否能在原操作／拼接中承载，以及相同物理来源能否辨认；不追加自由拟合相位以假装生成。空间旧接口、四分支及699范围继续保留。
'''
write('unified_physics_condition_ledger_710.md',ledger)
write('round710_drafts/research_note_710_draft.md',(HERE/'research_note_710.md').read_text('utf8'))
write('round711_drafts/STATUS.md','''# 第711轮入口：独立拓扑贡献与同一内部过程

接[710](../research_note_710.md)、[全账](../unified_physics_condition_ledger_710.md)。明确的QQQL扩展可以在原CAR／Gauss对象上承担荷变化、原记录及新热参考；但单个共同系数相位可消去，不能等同于629独立拓扑权重。

1. 先回查原有限图规范构形、区域拼接、连续全局群和拓扑权重；核已有对象中是否有独立可辨的拓扑载体，避免重复旧大规范／平坦接缝问题。
2. 保留原完整过程的实际物理来源，区分单个有效复系数、两份可干涉的贡献、拓扑历史／边界资料。若需新变量、限制构形或新幅度，明确列为输入，不用theta名称替代构造。
3. 保继续更新的认知要求只在声明的任务、态、精度及时间范围中验收；H0、710有界扩展、原手征候选和连续匹配仍分开。不要把710的存在证明当成反常生成，也不再追加任意QQQL系数凑轮次。
4. 旧空间382—386、425、522—524直接复用；384消去的Lipschitz不恢复，386和425是替代桥。原699仅否定指定反射正性接法，统一目标未完成。
''')
write('round710_drafts/literature_scope_audit.json',json.dumps(dict(date='2026-10-03',sources=[
    dict(url='https://arxiv.org/html/1405.0486',read_scope='Introduction Eq(1), Qqqql tensor',
        reused='standard gauge and spin contraction',not_claimed='original discovery, UV matching, measured baryon violation'),
    dict(url='https://arxiv.org/pdf/1402.6340',read_scope='Eqs(7)-(14), explicit B+L interaction and weak-angle interference',
        reused='independent phase must be counted jointly with topology',not_claimed='derived instanton amplitude, observed CP asymmetry, current experimental limits')],
    original_phase_matrix='629 exact A and C reused; phase sign follows629',
    local_normalization='664 psi=c/sqrt(w), w=epsilon^3 exp(6 sigma) reused'),ensure_ascii=False,indent=2)+'\n')
write('round710_drafts/scope_and_dedup_review.md','''# 710主代理去重与范围审查

上次目标轮709完成并发布710入口，属于进展。此次先读导航、最新报告、结果、入口与进程；无活动Python。QQQL收缩为成熟工具，629反常、664单元归一、665有界接触扩展、709条件期望及710入口中心商均复用。

新内容是原完整对象中的特定非零荷变化扩展、共同Gauss／热参考／读口／区域流、来源与几何系数的绑定，以及单系数相位无法提供独立拓扑资料的限定证明。不同分支不能混签。

逐项检查：创建／湮灭顺序和费米符号；36项及880个全部输出；原32位标的旁观模式；SL2只验张量不当CAR酉boost；多代偶乘积不当唯一瞬子张量；新H与新Gibbs替代而非原样保持；有界扰动固定正背景；源恒等式的符号；全局Z_n_g与局部边界的区别；相位商使用原629约定且严格限制参数族；单共同相位消去时准备与仪器同步；固定单位lapse的空间体积源不当四维Weyl或GR；裸系数增长不当UV普遍反证。

未安装新运行时、未做图像检查、未创建应用任务、未更改目标。没有独立代理审查。数值只校准实际系数／CAR身份，不冒充完整无界热谱或连续极限。旧空间接口逐项保留，699失败范围未扩大。
''')
review='710主代理最终审查完成；没有独立代理审核。\n'
for name in ('research_note_710.md','joint_charge_changing_vertex.py','joint_charge_changing_vertex_results.json','unified_physics_condition_ledger_710.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round710_drafts/final_review.txt',review)
print('Prepared710 ledger, review and711 continuation.')
