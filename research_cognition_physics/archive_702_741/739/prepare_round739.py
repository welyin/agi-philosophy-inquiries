"""Prepare739 common covariant response, finite-band proof and publication."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def remap(text,values):
    return re.sub('|'.join(re.escape(k) for k in sorted(values,key=len,reverse=True)),
                  lambda m:values[m.group()],text)


ledger=(HERE/'unified_physics_condition_ledger_738.md').read_text('utf8')
ledger='# 联合条件总账：739同一中性—几何方程与线性约束\n'+ledger.split('\n',1)[1]
ledger=ledger.replace('接[737全账](unified_physics_condition_ledger_737.md)，回填[738报告](research_note_738.md)。[结果](joint_mixed_neutral_response_results.json)、[核验](research_round_738_checks.json)。',
                      '接[738全账](unified_physics_condition_ledger_738.md)，回填[739报告](research_note_739.md)。[结果](joint_covariant_response_closure_results.json)、[核验](research_round_739_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**739当前增量：** 同一匹配常真空的两中性、共形和剪切响应，连同允许局部作用jet，组成协变12分量方程。实际互补非退化、零过去相容源及固定有限空间频段下有模坐标规范的短时线性退迟解，线性引力约束共同保持；不扩成全内部规范、任意尺度、非线性或稳定性结论。\n\n'+marker)
updates={'C04':'739限定中性—几何方程的有限频段短时退迟闭合；有效分支及非线性仍需验收',
         'C19':'739同真空的质量与剪切共用同一参考，真实装置来源相容性另核',
         'C22':'739局部jet、非局部交叉和全度规Ward共同接通；约束保持限声明线性部门'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 739共同退迟方程与有限频段

- 原632匹配常真空上，738中性矩阵谱与631剪切通式使用相同质量；全部接触和允许局部作用jet共同保持。
- 12分量协变符号保四个坐标规范零方向；重构8分量离壳来源块，不把五剪切当五个物理引力子。
- 互补主部非退化时，局部项真实落入737合同，不再单独假设其导数阶数已经合格。
- 非零空间动量通过精确阈值位移接入原谱；常数矩阵不变，记忆余项在每个固定有限K内有共同短时界。
- 波动逆积分与逆对数共同给零过去相容来源的短时线性解；四个线性引力约束由同一Ward身份传播。
- 不重复732外源经典发展，也不把当前反馈发展扩为730变化背景的非线性自洽。
- K是适用范围，不是新认知阈值。局部系数可令可证T极短，数学存在不自动满足共同EFT有效窗口。

## 本轮合并与下一项

C04／C19／C22在限定连续部门中合成一份退迟方程；统一目标仍开放。

接[740](round740_drafts/STATUS.md)：核同一总符号的物理分支、极点及601共同有效阶次，之后再推进任意空间频率、实际动态参考及非线性。旧空间、604、649／699保持。
'''
write('unified_physics_condition_ledger_739.md',ledger)
write('round739_drafts/research_note_739_draft.md',(HERE/'research_note_739.md').read_text('utf8'))
write('round740_drafts/STATUS.md','''# 第740轮入口：同一总反馈的极点与共同有效分支

接[739](../research_note_739.md)、[条件账](../unified_physics_condition_ledger_739.md)。限定中性—几何部门的固定有限空间频段已有相容源的短时线性退迟解及引力约束传播。

1. 回查601、630—632、732、736—739及成熟半经典稳定性结果，不重复一般高阶runaway或纯谱正性。
2. 计算同一总符号，保632实际势Hessian、738原中性混合核和739引力混合；不能用纯I核无零点推断总反馈无增长极点。
3. 区分把一圈有效作用精确重求和、严格共同圈阶展开、以及低能降阶。前者的高频解不自动成为原完整物理预测。
4. 先给具体原参数和处方下的可复算证据，再区分具体分支失败、给定尺度外失效与统一计划的否定。
5. 有限局部系数与有效能标未由认知原则决定，保持显式；不靠扫系数找一个好结果后宣布唯一。
6. 之后接730实际参考、装置来源及非线性初值。旧空间、604、649／699范围与统一目标保持。
''')
write('round739_drafts/literature_scope_audit.json',json.dumps(dict(sources=[dict(
    url='https://arxiv.org/html/gr-qc/0209075',authors='Anderson, Molina-Paris, Mottola',
    locations='Appendix A equations55-60; local renormalization around equation30; scalar-example stability scope',
    use='Covariant tensor projectors and common local subtraction; no importing scalar-field stability into the original fermionic/mass response.',
    boundary='Finite-band short-time response is distinct from gauge-invariant long-time stability and from nonlinear semiclassical existence.')],
    new='Original local jets and complete neutral/metric symbols, exact spatial-threshold shift, unchanged logarithmic constant and uniform finite-band causal closure.',
    inherited='583 local frame mixing,601 order reduction,631 shear,632 matching,732 Ward propagation mechanism,735 chosen conserving prescription,737-738 response kernels.',
    excluded='All internal gauge dynamics, arbitrary Cauchy data, UV-uniform estimate, actual changing reference, nonlinear spacetime or stable physical branch.'),ensure_ascii=False,indent=2)+'\n')
write('round739_drafts/scope_and_dedup_review.md','''# 739范围、证明与去重审查

前一目标轮次737—738完成并出版，归类为进展。当前先核最新导航、报告、结果和进程，无活跃Python。主代理审查，无子代理或图像检验。

- 732已有冻结相对源的一阶经典发展；本轮非局部响应在左侧，不把同一外源定理换名重算。
- 全部谱使用632同一背景；剪切不复用630 q★的数值。局部总势的零值和零梯度由632已声明匹配条件保证。
- 独立Fierz—Pauli和线性Riemann/Ricci/R给局部曲率算符，再与投影公式比较。
- Euclidean符号作解析双线性延拓，代码不用Hermitian转置替换普通转置。局部α、β不是实时间稳定性符号结论。
- Γ_E二阶负连通核的符号与630一致，质量和剪切不同减除阶的号分别保留。
- 全12分量Ward及局部交叉保持；内规范矢量和角向未被宣称已共同求解。
- 非零空间动量的正测度位移是精确身份；静态F(k²)必须保留，才能得到常数B不变。
- 余核点态界在固定K内统一，可进入C_t L²_k的Volterra估计，不只逐模式存在。
- 时间投影中的逆波动算子使用零过去J_k，不增加任意齐次解；P2约束保持避免依赖随频率变化的偏振基定义演化。
- 零过去相容来源给线性引力约束保持；不替代任意高阶初值、实际装置制备或非线性可积。
- 数值q仅含纯非局部核。完整局部项可证窗口与物理EFT窗口是否重合尚待下一轮。
- 旧空间、604、649／699及全部冻结文件保持；目标未改，不设置定时任务。
''')
main=('research_note_739.md','joint_covariant_response_closure.py','joint_covariant_response_closure_results.json','unified_physics_condition_ledger_739.md')
write('round739_drafts/final_review.txt','Primary-agent review only. Common stationary neutral/metric retarded linear system and finite-band constraints; no UV-uniform or nonlinear conclusion.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'739':'740','738':'739','737':'738','3415':'3417','3413':'3415',
         '1526':'1529','3336':'3345','3327':'3336','384':'385','289':'290',
         'joint_mixed_neutral_response':'joint_covariant_response_closure'}
publication=remap((HERE/'publish_round738.py').read_text('utf8'),mapping)
publication=re.sub(r'^summary=.*$',"summary='**第739轮完成：** [同一中性—几何退迟方程与线性约束]({p}research_note_739.md)原局部jet与完整谱组成协变方程；互补非退化、零过去相容源和固定有限空间频段下，短时线性退迟解与引力约束共同成立。两组、十八式通过，最新739／3417，1529份编号科学文件、3345份保护证据。[核验]({p}research_round_739_checks.json)、[条件账]({p}unified_physics_condition_ledger_739.md)。内部规范、有效分支、任意尺度及非线性仍开放。'",publication,flags=re.M)
publication=re.sub(r'^order=.*$',"order='**当前执行顺序（739后，优先于下方历史安排）：** 接[740总反馈与共同有效分支]({p}round740_drafts/STATUS.md)，核同一总符号的极点和601有效阶次；纯谱正性及短时可解不保证物理稳定。随后接实际动态参考。旧空间、604、649／699及统一目标保持。'",publication,flags=re.M)
publication=publication.replace('实际中性混合质量谱与完整矩阵因果逆','同一中性—几何退迟方程与有限空间频段的约束闭合')
publication=publication.replace('旧空间合同保持，矩阵响应不替代完整约束','旧空间合同保持，线性约束不替代非线性自洽')
write('publish_round739.py',publication)
write('postcheck_round739.py',remap((HERE/'postcheck_round738.py').read_text('utf8'),mapping))
verification=remap((HERE/'verify_round738.py').read_text('utf8'),mapping)
verification=verification.replace("checks['display_formulas']==16","checks['display_formulas']==18")
verification=verification.replace('original_unequal_mass_projectors_and_common_scaling_projection_checked=True','original_local_curvature_and_spatial_spectral_shift_checked=True')
write('verify_round739.py',verification)
print('Prepared739: common symbols, finite-band proof, ledger, review and next scope.')
