"""Prepare696 publication; no historical or frozen entry artifacts changed."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)
def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_695.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：695原辅助逐配置非负性的严格反例',
                      '# 联合条件总账：696完整辅助平均后的固定接缝边界核')
ledger=ledger.replace('接[694全账](unified_physics_condition_ledger_694.md)，回填[695报告](research_note_695.md)。[结果](joint_auxiliary_sign_certificate_results.json)、[核验](research_round_695_checks.json)。',
    '接[695全账](unified_physics_condition_ledger_695.md)，回填[696报告](research_note_696.md)。[结果](joint_dynamic_auxiliary_integral_results.json)、[核验](research_round_696_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 696辅助平均本身不能修复固定接缝历史核

- 原一空间周期自环、两AP时间片、16通道及256维N，半区空间链路分别为I与原弱中心。不同半区组合有中心时空小回路；不含不同空间点之间的传播，不替代695四时空节点背景。
- 同一原辅助Pfaffian相对参考精确为D_m八次方；完整两个S9积分由有限有理矩算完。三份原标量值分别属于Q、Q(√2)、Q(√5)，保留完整相位与物理prefactor。
- 固定两时间接缝为I的完整辅助边界核有明确有理负测试向量(1,−1/5)，范数严格小于−5×10⁻¹⁰；非对角／对角几何均值之比严格大于28。排除只做S9平均就必修复该裸核的接法。
- 原时间Haar与H_b尚未纳入，此结果不是原物理RP反例，也不否定671静态正性。入口120维局部重排筛查仅保留为方法边界，没有独立新增一轮。
- 两个半区取规范不变标量点及中心空间链路时，完整S9平均允许沿673同步输送消去一份边界Haar，保留最终Ω=ab。原H_b因子成为kτ(q_i,q_j)²；剩余全原群holonomy平均未完成。
- 旧C07—C09／C20直接复用382—386、425、522—523；不重列已消去的Lipschitz或共同实现，不把外幂次数、自环或辅助E当位置尺度、空间或实际s记录。

## 本轮合并与下一项

C01、C16、C19的实际来源现有一份完整辅助积分基准及严格固定接缝负核，须由完整规范闭合和原H_b共同验收。不能把三份正标量积分当作正历史过程。

接[697](round697_drafts/STATUS.md)：核指定中心半区的实际全群holonomy来源及完整积分或受控界。复用673换元与675必要条件，不重复固定接缝负证书，不把子群平均改名原群。原物理RP、指定归一、H_F同一性、实际空间传播、共同连续及量子GR仍开放；目标不改。
'''
write('unified_physics_condition_ledger_696.md',ledger)
mapping={'joint_auxiliary_sign_certificate':'joint_dynamic_auxiliary_integral',
         '694':'695','695':'696','696':'697','3295':'3297','3297':'3299',
         '1394':'1397','1397':'1400','2681':'2700','2700':'2720'}
verify=replace((HERE/'verify_round695.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_critical_source_certificate as prior_model','import joint_auxiliary_sign_certificate as prior_model')
verify=verify.replace('Verify696 exact opposite-sign original auxiliary configurations.',
                      'Verify696 full auxiliary integrals and strict fixed-seam negative kernel.')
names=('unified_physics_condition_ledger_696.md','round696_drafts/research_note_696_draft.md',
       'round696_drafts/final_review.txt','round697_drafts/STATUS.md',
       'round696_drafts/exterior_average_entry.py','round696_drafts/exterior_average_entry_results.json',
       'round696_drafts/exterior_average_entry.md','round696_drafts/check_and_publish_entry.py',
       'round696_drafts/entry_checks.json','round696_drafts/prepare_entry_publication.py',
       'round696_drafts/dynamic_center_probe.py','round696_drafts/dynamic_center_probe_first.py',
       'round696_drafts/joint_dynamic_auxiliary_integral_first.py',
       'round696_drafts/joint_dynamic_auxiliary_integral_first_results.json',
       'round696_drafts/integral_first_saved_result.json','round696_drafts/integral_kernel_first_results.json',
       'round696_drafts/integral_source_audit.json')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
write('verify_round696.py',verify)
pub=replace((HERE/'publish_round695.py').read_text('utf8'),mapping|{'## 341.':'## 342.','## 246.':'## 247.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第696轮完成：** [完整辅助平均后的固定接缝负核]({p}research_note_696.md)'
         '原中心空间历史的全部S9积分精确求出；固定时间接缝的2×2标量核仍有严格负方向。'
         '两组、十六式通过，最新696／3299，1400份编号科学文件、2720份保护证据。'
         '[核验]({p}research_round_696_checks.json)、[全条件账]({p}unified_physics_condition_ledger_696.md)。'
         '时间Haar与H_b尚未纳入，不构成原物理RP或H_F反例；单空间自环不替代实际空间传播。')
order=('**当前执行顺序（696后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[697中心边界的完整规范闭合]({p}round697_drafts/STATUS.md)，'
       '复用673同步输送，保留最终全群holonomy及原H_b；'
       '旧空间接口逐项复用，目标不改。')
"""+pub[b:]
pub=pub.replace('原辅助逐配置非负性的严格反例','完整辅助平均后的固定接缝负核')
pub=pub.replace('原辅助符号与完整物理平均的边界','完整辅助平均与尚待完成的规范闭合')
pub=pub.replace('旧空间合同复用，辅助E不是实际方向仪器','旧空间合同复用，单空间自环不等于空间生成')
write('publish_round696.py',pub)
write('postcheck_round696.py',replace((HERE/'postcheck_round695.py').read_text('utf8'),mapping))
with (HERE/'round696_drafts/research_note_696_draft.md').open('xb') as f:f.write((HERE/'research_note_696.md').read_bytes())
review='''696 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed695 and executed696 entry: progress. Navigation and saved results read.
Get-Process found no live Python. Reviewed653/657/671/673/675 and the original physical factor.
Entry120-dimensional indefinite realignment is only screening, not a new round or physical RP proof.
Specified original finite volume: one periodic spatial loop, two AP times, all16 channels and N256.
Spatial weak-centre histories have a nontrivial offdiagonal plaquette but no intersite spatial propagation.
Actual small Wilson blocks depend on mean m and scalar shift delta; signs remain opposite.
Both actual64-mode auxiliary Pfaffian factors retained; no square-root determinant phase chosen.
Accidental Spin6 x Spin4 stabilizer checked with original21 integer generators, not new gauge group.
Central Clifford element Q has +/-i multiplicity8. Each restricted quadratic determinant has multiplicity4.
This gives two identical D_m^4 factors and normalized auxiliary D_m^8 at all sphere endpoints by continuity.
Original physical chiral factor P_m positive; reference phase follows nonzero real matrix homotopy to free.
Finite complete sphere moments computed in Fraction for m0,1,2. Free m0 agrees with inherited653.
Original256 Pfaffians, compressed32 determinants and independent complete four-variable quadrature agree.
Eigenframe phase normalized by same nonzero reference, not absolute value; variable-shadow bug caught before save.
All first computations retained. Expanded result field distinguishes offdiagonal plaquette from diagonal cases.
Rational sqrt2/sqrt5 enclosures give strict negative test(1,-1/5), bound<-5e-10, offdiagonal ratio>28.
No fixed E remains, but temporal seams remain fixed. Not a full Haar/Hb physical RP counterexample.
673 boundary covariance removes one Haar only at gauge-fixed central half-configurations after full E averaging.
Remaining holonomy integral and full Hb factor not computed or discarded.675 full-average failure test not invoked.
Old382/383/384/386/425/522/523 reused, no new dimension claim; removed Lipschitz not restored.
Two new groups,16 displayed equations,no images; next697 complete group closure, not repeating bare-kernel test.
'''
for name in ('research_note_696.md','joint_dynamic_auxiliary_integral.py','joint_dynamic_auxiliary_integral_results.json','unified_physics_condition_ledger_696.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round696_drafts/final_review.txt',review)
print('696 publication preparation completed')
