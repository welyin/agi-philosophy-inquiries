"""Prepare678 publication with preservation of all prior evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_677.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：677全谱调节与完整物理来源的共同极限',
    '# 联合条件总账：678完整来源的局部辅助表示与共同归一')
ledger=ledger.replace('接[676全账](unified_physics_condition_ledger_676.md)，回填[677报告](research_note_677.md)。[结果](joint_rational_physical_limit_results.json)、[核验](research_round_677_checks.json)。',
    '接[677全账](unified_physics_condition_ledger_677.md)，回填[678报告](research_note_678.md)。[结果](joint_local_source_lift_results.json)、[核验](research_round_678_checks.json)。')
ledger=ledger.replace('677全谱调节下原物理来源与完整平均的共同极限|',
    '677全谱调节下原物理来源与完整平均的共同极限；678全部来源与原质量的有限范围体表示及归一抵消|')
ledger=ledger.replace('677完整未归一来源在固定有限盒与正热时收敛；',
    '677完整未归一来源在固定有限盒与正热时收敛；678保留全部来源的局部辅助表示及相位；')
ledger=ledger.replace('任意费米来源反射、动态正性、指定归一与连续局域性仍缺',
    '678有限体测度与全部物理来源的表示身份已证，体因子可局部抵消；任意费米来源反射、动态正性、非零物理归一与连续局域性仍缺')
ledger=ledger.replace('动态量子反馈、一般几何求导与共同重整化|',
    '678原体行列式与观测字典的有限背景响应须共同抵消；动态量子反馈、一般几何求导与共同重整化|')
ledger=ledger.replace('这不移除Wilson H谱隙、固定index、动态正性和整体归一要求。',
    '670当时按无零Wilson谱片书写；673后来已用矩形字典跨不同index，并以几乎处处谱与热矩建立有限积分。不能重新要求固定index或统一硬谱隙；动态正性、非零物理归一及连续局域性仍需另证。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 678局部约束体精确保留原非零质量、全部物理Grassmann来源与相位，不只一个有理传播子。辅助Q逆的zz块为零，所以没有额外的来源接触指数。
- 体因子det K由同一原Wilson背景确定；补偿费米K†与正定交换变量K†K精确抵消它。局部有限范围与玻色高斯收敛不是反射正性。
- 先积分全部辅助变量，再使用677原完整平均极限；不对维数随L变化的未积分体测度直接套支配收敛。
- 677末段索引曾列C26，但无新的宏观／流体结果；本轮明确不推进C26。旧空间、热参考及实际仪器条件仍按既有范围继承。

## 本轮合并与下一项

678连接C01、C14、C16、C17、C20、C22：原物理来源、质量、规范背景、局部辅助表示和归一响应共同确定。排除的是忽略体因子／来源变化的表示，不是原认知原则或整个物理候选；新增辅助变量不认作新粒子或主体。

当前候选的有限局部表示已经落实，但一般动态规范正性、非零物理归一、原H_F与物理时间、实际记录、共同连续极限及量子GR仍未完成。原固定图过程、连续物质、经典动态几何和手征辅助测度四个分支保持分别记账。

接[679](round679_drafts/STATUS.md)：回到同一物理半区来源，复用661自由Weyl反射、668质量与674标量实性，核非平坦动态背景的二／四来源反射字典及有限调节／原极限区别。完成体表示后不继续层数优化；目标及旧空间合同不改。
'''
write('unified_physics_condition_ledger_678.md',ledger)
mapping={'joint_rational_physical_limit':'joint_local_source_lift',
         '676':'677','677':'678','678':'679','3259':'3261','3261':'3263',
         '1340':'1343','1343':'1346','2460':'2473','2473':'2485'}
verify=replace((HERE/'verify_round677.py').read_text('utf8'),mapping)
verify=verify.replace('Verify full-spectrum physical-source limit and preserved historical scope.',
    'Verify complete local source lift, normalization and preserved evidence.')
verify=verify.replace('import joint_nonflat_seed_support as prior_model',
    'import joint_rational_physical_limit as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_678.md','round678_drafts/research_note_678_draft.md',
           'round678_drafts/final_review.txt','round679_drafts/STATUS.md',
           'round678_drafts/local_block_probe.py','round678_drafts/local_block_probe_results.json',
           'round678_drafts/local_block_entry.md','round678_drafts/check_entry.py',
           'round678_drafts/entry_checks.json')
"""+verify[b:]
write('verify_round678.py',verify)
publish=replace((HERE/'publish_round677.py').read_text('utf8'),mapping|{'## 323.':'## 324.','## 228.':'## 229.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第678轮完成：** [完整物理来源的局部体表示与共同归一]({p}research_note_678.md)'
         '同一原Wilson约束体精确保留全部物理来源、原质量与复相位；'
         '体行列式由局部补偿费米和收敛玻色积分抵消，原背景响应同步。'
         '两组、十六式通过，最新678／3263，1346份编号科学文件、2485份保护证据。'
         '[核验]({p}research_round_678_checks.json)、[全条件账]({p}unified_physics_condition_ledger_678.md)。'
         '正性、非零物理归一、原H_F时间、连续及量子GR仍开放；旧空间接口复用。')
order=('**当前执行顺序（678后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[679原物理来源的反射合同]({p}round679_drafts/STATUS.md)，'
       '核动态规范背景及有限调节／原极限区别，不继续层数优化；目标不改。')
'''+publish[b:]
publish=publish.replace('全谱有理表示与完整物理来源平均的共同极限','完整物理来源的局部体表示与共同归一')
publish=publish.replace('完整平均的受控极限与旧空间接口复用','完整局部表示与实际物理过程的区别')
publish=publish.replace('有限调节极限不替代实际空间与物理过程','局部辅助体不增加空间维数或物理主体')
write('publish_round678.py',publish)
write('postcheck_round678.py',replace((HERE/'postcheck_round677.py').read_text('utf8'),mapping))
write('round678_drafts/research_note_678_draft.md',(HERE/'research_note_678.md').read_text('utf8'))
review='''678 primary-agent code/proof/scope review; no independent agent.
Prior677 published and678 entry executed: progress. No running old Python.
Read navigation/state/report/results and related659-661,665,668,673-677.
Old382-386,425,522-523 inherited; no Lipschitz reinstatement or route conflation.
K recurrence same original Wilson X, full16 channels and gauge links.
R=Aplus^-1 Aminus=B^-1 T B, order checked; K invertible for0<a<1.
detK=(detAplus)^L detB det(I+R^L); positive here via accretivity/continuity.
Psi,chi original; z,eta auxiliary constraint only. Exact Grassmann shift
Jacobian1; eta saturates every shifted z, so all pure physical sources retained.
Berezin sign nu=(-1)^(m(m+1)/2); original m multiple64 gives+1.
Local source Y reduces exactly to677 Phi, kinetic factor1/2 retained.
Mass zz block included. Q^-1 zz=0 means no missing source-source contact.
Schur action and source match; no inverse of physical weight, mass or D needed.
Extra ghost integral gives nu detKdagger, positive complex Gaussian gives
1/det(Kdagger K); product with numerator gives original signed weight.
Compensation works for complex nonzero K, not an absolute-Pf replacement.
Kdagger K positive ensures finite Gaussian convergence, not reflection positivity.
Finite gauge configuration compact supplies fixed-regulator lower bound only.
Local finite-range on declared auxiliary cycle; no physical space derivation.
All gauge fields, M, phi, Y and sources transported; original doubleHaar preserved.
Derivative of determinant is retained, regulator background source cancels.
No declaration of new particles, physical clocks or cognitive architecture.
Auxiliary dimension changes withL: integrate first, apply677 afterwards.
No arbitrary source/limit interchange, continuum, original HF or quantumGR claim.
Full1536/1538/1540 Pf with original variables, phases and2/4 sources verified.
Log Pf uses615 ordering to avoid overflow; direct512 Pf cross-check.
Flux full source/action identity checked without division by original zero target.
677 C26 list had no new macro result; current ledger preserves C26 gap.
670 former index scope superseded by673 rectangular/integrability result.
Next679 actual physical source reflection, not layer-count optimization.
2 groups,16 equations. All frozen history retained.
'''
for name in ('research_note_678.md','joint_local_source_lift.py','joint_local_source_lift_results.json','unified_physics_condition_ledger_678.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round678_drafts/final_review.txt',review)
print('678 preparation completed')

