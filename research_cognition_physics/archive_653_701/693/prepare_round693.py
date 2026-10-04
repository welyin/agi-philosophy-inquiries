"""Prepare693 evidence and publication; preserve all frozen prior artifacts."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    keys=sorted(mapping,key=len,reverse=True)
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in keys]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_692.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：692原辅助轨道与全部物理来源的精确有限求积',
                      '# 联合条件总账：693完整规范平均的近零谱与来源误差')
ledger=ledger.replace('接[691全账](unified_physics_condition_ledger_691.md)，回填[692报告](research_note_692.md)。[结果](joint_auxiliary_orbit_quadrature_results.json)、[核验](research_round_692_checks.json)。',
    '接[692全账](unified_physics_condition_ledger_692.md)，回填[693报告](research_note_693.md)。[结果](joint_gauge_sublevel_control_results.json)、[核验](research_round_693_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 693原全Haar平均的定量控制

- 原2空间点×2AP时间片、两传播方向保持不变。单位时间边界上，任意空间规范链路都有见证谱隙2−√3；这不是对全部配置添加谱隙。
- 原SU3／SU2／U1的Cayley坐标和614荷模块给每群元22次正分母，256维原Wilson行列式次数至多22528。Remez界给原双Haar近零谱概率的统一幂上界，指数1/1306624，不另加谱密度或横截性。
- Pfaffian扰动界覆盖奇异端点及全部固定阶原物理来源；669热矩接成完整配置／S9／双Haar的误差式。677／686的定性收敛直接继承，不作为新定理重复领取。
- 同一个上界控制692九点调节表达趋向精确物理泛函；有限软投影不满足精确九点恒等式，二者明确分开。
- 当前界极保守，热矩系数也未数值求出。尚无实际完整积分的正负号证书、指定非零归一、动态RP、H_F身份或共同连续／量子GR。

## 用户指定的空间旧接口：直接复用与实际剩余映射

|历史|直接继承|当前仍须接到原共同对象的内容|
|---|---|---|
|[382](research_note_382.md)|完整反向qubit方向接口给三维上界；上下界合并的条件性三维|真实完整位置边界、端口身份和方向效果不能由任意Bloch标签代替|
|[383](research_note_383.md)|连续自由反向对合足够，无须先与标准对径共轭|当前实际端口上的两次回归和无固定点|
|[384](research_note_384.md)|下界中的额外Lipschitz已消去；可逆全族极小作用够用|有限尺度认证仍须有明确变化模，不能重加已消去的定性假设|
|[386](research_note_386.md)|渐近位移连接真实邻域，已有有限尺度维数证书|将实际端点／缩放对象映入该定理，保持覆盖及误差量词|
|[425](research_note_425.md)|一致半幅、成本收缩等已声明条件给光滑坐标|这是与386可选的坐标路线，不把两套充分条件强行叠加|
|[522](research_note_522.md)|原热参考资源下界|保持共同／相对参考模及原热态定义|
|[523](research_note_523.md)|实际方向仪器与条件性三维已有共同实现|它与当前h、s、CAR／Gauss物质过程的同一映射；E辅助字段不是实际s记录|

本轮近零谱概率指数是有限矩阵积分的模量，不是空间维数、光滑度或跨尺度一致性。未重证这些空间结果，也未解决该表右列的对象匹配。

## 本轮合并与下一项

C01物理来源、C14完整Gauss、C19未归一泛函接成同一有限盒误差控制；C07—C09旧空间与参考结果原样复用。非零归一及正性仍为独立门槛。

接[694](round694_drafts/STATUS.md)：核原来源在实际谱穿越附近是否额外消去，或存在不被来源抑制的跳变。不得靠裸谱常数优化、层数扫描或无误差随机积分替代实际物理二次型。目标不改。
'''
write('unified_physics_condition_ledger_693.md',ledger)
mapping={'joint_auxiliary_orbit_quadrature':'joint_gauge_sublevel_control',
         '691':'692','692':'693','693':'694','3289':'3291','3291':'3293',
         '1385':'1388','1388':'1391','2641':'2655','2655':'2668'}
verify=replace((HERE/'verify_round692.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_neutral_symmetry_reduction as prior_model',
                      'import joint_auxiliary_orbit_quadrature as prior_model')
verify=verify.replace('Verify693 original auxiliary orbit reduction and complete-source quadrature.',
                      'Verify693 original Haar sublevel bound and full-source error control.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=('unified_physics_condition_ledger_693.md','round693_drafts/research_note_693_draft.md',
       'round693_drafts/final_review.txt','round694_drafts/STATUS.md',
       'round693_drafts/finite_character_entry.py','round693_drafts/finite_character_entry_results.json',
       'round693_drafts/finite_character_entry.md','round693_drafts/check_and_publish_entry.py',
       'round693_drafts/entry_checks.json','round693_drafts/primary_source_audit.json')
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==16")
write('verify_round693.py',verify)
pub=replace((HERE/'publish_round692.py').read_text('utf8'),mapping|{'## 338.':'## 339.','## 243.':'## 244.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第693轮完成：** [完整规范平均的近零谱与原来源误差]({p}research_note_693.md)'
         '给出对全部原空间规范背景统一的近零谱Haar幂界，接成原完整物理来源及九点调节表达的误差式。'
         '两组、十六式通过，最新693／3293，1391份编号科学文件、2668份保护证据。'
         '[核验]({p}research_round_693_checks.json)、[全条件账]({p}unified_physics_condition_ledger_693.md)。'
         '固定有限盒，界很松、热矩系数未求值；不构成实际符号、RP或连续时空证书。')
order=('**当前执行顺序（693后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[694临界谱附近的原物理来源]({p}round694_drafts/STATUS.md)，'
       '核实际来源的消去或跳变，不做裸谱常数优化；旧空间382—386、425、522—523逐项复用，目标不改。')
"""+pub[b:]
pub=pub.replace('完整物理泛函的辅助轨道与九点精确求积','完整规范平均的近零谱与原来源误差')
pub=pub.replace('完整物理泛函的辅助轨道与精确求积','完整规范平均的近零谱与来源控制')
pub=pub.replace('旧空间合同复用，辅助轨道代表不是物理参考','旧空间合同复用，有限谱模量不是空间尺度')
write('publish_round693.py',pub)
write('postcheck_round693.py',replace((HERE/'postcheck_round692.py').read_text('utf8'),mapping))
with (HERE/'round693_drafts/research_note_693_draft.md').open('xb') as f:
    f.write((HERE/'research_note_693.md').read_bytes())
review='''693 primary-agent proof/code/scope review; no independent agent.
Previous completed692 and executed693 entry: progress. Active goal unchanged.
Original2space/2APtime/two-direction box, all16 internal and all auxiliary channels retained.
Unit time seam Xs=I-Us gives unitary Us and uniform witness gap2-sqrt3 for all spatial links.
This is not a global hard gap. General full H norm<=3 inherited677.
Cayley Haar Uj, SUj determinant removal and original quotient all preserve correct Haar pushforward.
Density normalization and Remez theorem checked in primary papers; derived tail uses marginal angles, not independent angles.
Original614 left modules clear with degree22 positive denominator; n256 K4 M56 D22528.
Polynomial witness uniform in spatial links; probability power1/1306624 is explicit but extraordinarily weak.
Pfaffian derivative bound valid for complex skew matrices including singular endpoints; no propagator inverse.
Original source bounds polynomial in original scalar weight;669 heat moments finite but not numerically evaluated.
All q/S9/doubleHaar retained. Full-source approximation only at fixed box, positive heat time and finite source order.
Finite soft epsilon not projector: no675 factorization or692 exact-degree claim at finite a,L.
Nine-node soft expression separately defined; uniform source bound makes it converge to692 exact target.
No new La² condition. No normalized-ratio bound without nonzero partition function.
Old382/383/384/386/425/522/523 reviewed and mapped in ledger; deleted Lipschitz not restored.
386 and425 remain alternative sufficient coordinate routes.523 joint instrument model not re-proved.
Current h/s/CAR/Gauss mapping still missing; E is not s record.
Two new groups and16 display equations; no images; no RP/HF/continuum/GR completion.
Next694 actual-source critical behavior; no precision-only repetition.
'''
for name in ('research_note_693.md','joint_gauge_sublevel_control.py','joint_gauge_sublevel_control_results.json','unified_physics_condition_ledger_693.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round693_drafts/final_review.txt',review)
print('693 publication preparation completed')
