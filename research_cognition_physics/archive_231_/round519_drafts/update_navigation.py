"""Snapshot and compare seven mutable navigation documents before updating."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checks=json.loads((HERE/'research_round_519_checks.json').read_text('utf8'))
assert checks['all_reported_checks_passed']
assert checks['fresh_tests']==dict(run=4,failures=0,errors=0)
assert checks['previous_protected_evidence_hashes_verified']==996
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
snapshot=HERE/'round519_drafts/navigation_before_round519'
snapshot.mkdir(exist_ok=False)
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
            'archive_231_/' if p.parent==RESEARCH else '')
    banner=(f'**第519轮完成：** [边界换边与可分辨返回预测]({prefix}research_note_519.md)'
        '在518同一实际来源与原持续H下，证明外部NNI进入返回信号四阶首项；'
        '末端角色已足以给修正预测，无需额外局部量子摘要。'
        '全规模／未知图参考下给有限误差小于非零信号的证书，孤立星预测有严格差距。'
        '时间与来源率极端保守，读者、计时和实现仍输入，未生成宏观坐标。'
        '4项、12式及独立终审通过；累计2546项、868份编号科学文件、1000份保护证据。三维与GR仍未完成。')
    correction=(f'**518后范围修正（不增轮次）：** [已知局部态及证书分辨能力]({prefix}record_determined_predictor_review.md)'
        '收窄旧导航“仍需局部量子摘要”的判断：518特定good准备已固定局部纯态，经典记录足以确定该预测输入；'
        '通用E₂容差却也允许零点击预测，不能据此签收定位。4项未编号诊断单列，旧992份证据不改，'
        '新增4份证据随519完整继承。')
    first,rest=body.split('\n',1)
    body=first+'\n\n'+banner+'\n\n'+correction+'\n'+rest
    if p==RESEARCH/'research_direction.md':
        old='最新科学轮次与检查数为518／2542'
        assert body.count(old)==1
        body=body.replace(old,'最新科学轮次与检查数为519／2546')
        body=body.replace('完成231—518轮。','完成231—519轮。',1)
        body+='''

**519后的当前接续：** 518特定准备不再要求额外获取局部量子摘要；519已给包含边界作用的非零信号证书。后继返回有实际来源的宏观端点、跨尺度定位菜单及三维选择，不以更多Taylor阶、局部系数或率常数当作新轮。用户的宏观有效空间方向继续采用；感知错觉不直接证明微观几何异常，理想操作精确／极限语义保持原范围。

**519冻结与额外证据：** [科学检查](archive_231_/research_round_519_checks.json)、[复算入口](archive_231_/verify_boundary_corrected_return_round.py)、[未编号修正核验](archive_231_/record_determined_predictor_checks.json)。编号2546不包含新增4项未编号诊断。完整继承996份保护证据，再增3份科学文件和1份终审稿，共1000份；未生成三维或GR，不结项。
'''
    if p==RESEARCH/'RESEARCH_STATE.md':
        body=body.replace('已完成第231—518轮。','已完成第231—519轮。',1)
        body+='''

**519核验入口：** [实际返回与修正预测](archive_231_/boundary_corrected_return.py)、[保存结果](archive_231_/boundary_corrected_return_results.json)、[科学检查](archive_231_/research_round_519_checks.json)、[科学复算](archive_231_/verify_boundary_corrected_return_round.py)、[整合复核入口](archive_231_/verify_round519_integration.py)。全部原证据保留；未编号来源修正单列，不能把它算作第519的四项科学检查。
'''
    if p==HERE/'spatial_premise_closure_audit.md':
        body+='''

### 156.2 范围修正：518特定来源已决定局部纯态

上一节对未知局部摘要的警告在一般来源下有效，但对518过宽。树无环迫A内只有三边星，单激发末测w迫局部数据为已知基态，所以good条件态的A因子纯且与补因子／旧参考解耦。按实际经典记录CR边缘准备数学预测输入，即由498／499给E₂＋2√q；投影比较可改变R边缘，未假设不扰动。此直接推论保留为[未编号审查](record_determined_predictor_review.md)，不另造轮次。

继承E₂过宽：叶根真实点击始终≤min(1,s²)，而E₂处处不小于该界，故零点击预测器也合格。不能从这份容差声称已辨别传播。四项诊断、9式、独立完整复算通过；新4份证据在旧992份外继承。原笔记与本节前的历史表述保持原字节。

## 157. 第519轮：原边界作用与可分辨的实际返回

本回合先读导航、518和保存结果，发现并关闭上述摘要范围错误，再独立核验真实信号。已知局部态不能令跨界动力学消失：末端w为内部点时，沿a−u−w−d换接外支再跳跃，与w→u→a同为二阶振幅；其四阶概率系数由1／4变为3／4。逆向从u新增的唯一邻居恢复输入图，故作用效果是good支撑上的标量，覆盖全部未知图相干／参考。

[519](research_note_519.md)复用498嵌套交换子界，保留一般相干输入的五阶余项；将实际518来源两次温和比较计入，给单位耦合下s＝1／(64K⁵)、t＝s⁴／(64L)的全规模有限证书。修正预测用实际记录和末端角色，完整根CR误差≤s⁴／16，真实点击≥3s⁴／16；I≥2时实际点击比孤立星至少高s⁴／32。没有对某条记录后选择，也没有省掉源成功率t⁸。

四项检查覆盖整数全good作用、2520图反演无碰撞、精确有理预算、实际来源与二维旧参考诊断。代码、结果与12式最终稿均获独立审阅并完整复算，正式稿与终审稿原字节相同。累计2546项、868份编号科学文件；992原文件＋4未编号修正＋本轮4份科学／终审文件，共1000份保护证据。所有旧轮、失败反作用及冻稿保留，没有图像检验或新应用任务。

这轮补的是实际预测的分辨能力，不是宏观空间。读取、计时及极小误差实现仍是具体模型输入；不会把它们未经证明升级为认知原则。接下来检验有来源的共同端点／尺度及其几何组织，不把局部信号时间幂、三度角色或三个Bloch分量当空间维数。

### 157.1 已接续的去重与下一项边界

已回看491／493／517的实际统计来源、504的均值Hamiltonian障碍、506联合角色代码及宏观猜想审查。519不能消除这些旧结论，也不自动让细端口历史成为Markov过程。继续计算同一返回菜单更高Taylor阶，或将已知角色组成一个人为三坐标拟合，均不作为520。

后继须先给一个由实际过程选择的复合端点与尺度，明确其定位／传播报告和比较标准，再检验共同关系能否稳定串接及是否选择三维。微观载体无需各有光滑坐标，也不要求位置预测全部私人未来。当前完成了独立有效增量，尚无完整三维／GR结论，目标继续，不触发结项或阻塞判定。
'''
    if p==HERE/'three_dimensional_four_conditions_review.md':
        body+='''

## 62. 519：可分辨的局部返回仍未签收三维条件

518特定good来源的局部纯态可由记录确定，上一节“列表不能替代A内量子态”的说法按[范围修正](record_determined_predictor_review.md)收窄。[519](research_note_519.md)进一步纳入末端外支NNI，给未知图／参考及实际来源下误差小于信号的返回预测，不再仅依靠允许零输出的宽容差。

四项三维条件仍未关闭：局部紧的位置结构、一致半幅、保组合重定向及完整方向合同均不从该单个返回菜单直接得到。节点角色提供动力学系数，不是空间维数。编号519／2546，保护证据1000；三维与GR开放，阶段不收尾。
'''
    planned[p]=body.replace('\n',nl).encode(enc)
    label=('project_README.md' if p==ROOT/'README.md' else
           'research_'+p.name if p.parent==RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:f.write(raw)
for p,raw in baseline.items():assert p.read_bytes()==raw,str(p)
for p,data in planned.items():
    assert p.read_bytes()==baseline[p],str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},f,indent=2,ensure_ascii=False)
print(json.dumps({'navigation_files_updated':len(planned),'snapshots_preserved':len(baseline)}))
