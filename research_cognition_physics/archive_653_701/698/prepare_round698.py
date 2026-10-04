"""Prepare698 checked publication; preserve all prior evidence and probes."""
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

names=('unified_physics_condition_ledger_698.md','round698_drafts/research_note_698_draft.md',
    'round698_drafts/final_review.txt','round699_drafts/STATUS.md',
    'round698_drafts/polynomial_error_entry.py','round698_drafts/polynomial_error_entry_results.json',
    'round698_drafts/polynomial_error_entry.md','round698_drafts/check_and_publish_entry.py',
    'round698_drafts/entry_checks.json','round698_drafts/prepare_entry_publication.py',
    'round698_drafts/polynomial_error_entry_first.py',
    'round698_drafts/sphere_majorant_probe.py','round698_drafts/sphere_majorant_probe_results.json',
    'round698_drafts/block_majorant_probe.py','round698_drafts/block_majorant_probe_results.json',
    'round698_drafts/polynomial_majorant_probe.py','round698_drafts/polynomial_majorant_probe_first.py',
    'round698_drafts/joint_static_haar_majorant_first.py','round698_drafts/joint_static_haar_majorant_first_results.json', 'round698_drafts/research_note_698_before_link_fix.md')
protected=2738+3+len(names);assert protected==2761
ledger=(HERE/'unified_physics_condition_ledger_697.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：697完整规范平均的统一谱隙与数值候选',
                      '# 联合条件总账：698原静态中心项的完整积分上界')
ledger=ledger.replace('接[696全账](unified_physics_condition_ledger_696.md)，回填[697报告](research_note_697.md)。[结果](joint_full_holonomy_reduction_results.json)、[核验](research_round_697_checks.json)。',
    '接[697全账](unified_physics_condition_ledger_697.md)，回填[698报告](research_note_698.md)。[结果](joint_static_haar_majorant_results.json)、[核验](research_round_698_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 698原静态中心项的完整积分已有严格上界

- 同一697原256维来源精确降维后，以奇异值算术—几何均值给辅助行列式上界。实际十个T矩阵的支撑把20维球面二次型分为五个4×4块；完整两个S9的16次矩精确收缩。
- 对原中心极分解证明二次模界、跨片系数界和物理因子凸弦界，得到保持全部原群变量的正有理多项式上界。它是积分的上界，不是替换原测度或假设原逐配置正性。
- 原Weyl常数项由44阶网格、11个经试除确认的素数和整数界唯一恢复；48阶及新素数独立核同一值。完整B11满足|B11|≤U≈2.4548996203337804×10⁻¹⁰<2.46×10⁻¹⁰；U不是B11精确值。
- 既有B00精确值加本轮U，将一个充分负核判据集中到尚未证明的B10≥727/10¹⁶；向量(1,−3/10000)的条件性上界<−3.24×10⁻²⁰。697抽样估计不能签收这个下界。
- 入口157阶误差合同保留；本轮交付另一条可实际积分的上界。完整Bᴳ负号、原H_b有限τ、H_F身份及共同连续仍未完成。
- C07—C09按697 §1.1继承：384额外Lipschitz已消去；386／425为替代桥；522—523共同实现已完成，当前h、s、CAR／Gauss对应仍待建立。E不是s记录。

## 本轮合并与下一项

C01、C16、C19的同一有限候选取得完整原测度下的静态项严格控制，负核测试留下一个明确的非对角下界任务。没有新认知公理或空间假设，也未得到总体目标的反例。

接[699](round699_drafts/STATUS.md)：保原来源求B10的精确收缩或足够强的解析下界，复用有限Laurent结构、完整球面矩与群常数项。若下界不成立，记录真实结果；不能以扩大抽样替代证书。只有实际得到完整负号后才接675原H_b条件。目标与四分支保持不变。
'''
write('unified_physics_condition_ledger_698.md',ledger)
mapping={'joint_full_holonomy_reduction':'joint_static_haar_majorant',
    '696':'697','697':'698','698':'699','3299':'3301','3301':'3303',
    '1400':'1403','1403':'1406','2720':'2738','2738':str(protected)}
verify=replace((HERE/'verify_round697.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_dynamic_auxiliary_integral as prior_model','import joint_full_holonomy_reduction as prior_model')
verify=verify.replace('assert prior_model.run()==core.read(prior_model.TARGET)',
    "assert prior_model.algebra_checks()==core.read(prior_model.TARGET)['algebra']")
verify=verify.replace('uniform-gap original source reduction and full-measure numerical candidate.',
                      'exact full-Haar and sphere upper bound for the original static centre entry.')
verify=verify.replace("assert text['display_formulas']==16","assert text['display_formulas']==18")
a=verify.index('    names=');z=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[z:]
write('verify_round698.py',verify)
pub=replace((HERE/'publish_round697.py').read_text('utf8'),mapping|{'## 343.':'## 344.','## 248.':'## 249.'})
a=pub.index('summary=');z=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第698轮完成：** [原静态中心项的完整积分严格上界]({p}research_note_698.md)'
         '整个原群和两个S9平均后|B11|≤U<2.46×10⁻¹⁰，整数常数项已认证。'
         '负核判据仍缺B10的严格下界，未宣称完整RP失败。'
         '两组、十八式通过，最新698／3303，1406份编号科学文件、2761份保护证据。'
         '[核验]({p}research_round_698_checks.json)、[全条件账]({p}unified_physics_condition_ledger_698.md)。'
         '旧空间接口直接复用，目标及四分支不变。')
order=('**当前执行顺序（698后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[699原非对角完整积分的下界]({p}round699_drafts/STATUS.md)，'
       '核B10≥727/10¹⁶是否成立，保留原群、球面与H_b；'
       '数值估计不代替证书，旧空间合同继续复用。')
'''+pub[z:]
pub=pub.replace('完整规范平均的统一谱隙与数值候选','原静态中心项的完整积分严格上界')
pub=pub.replace('完整规范平均的解析约化与数值边界','原静态项的完整积分上界及剩余非对角条件')
pub=pub.replace('旧空间合同复用，全群谱隙不等于连续空间','旧空间合同复用，积分界不新增空间条件')
write('publish_round698.py',pub)
write('postcheck_round698.py',replace((HERE/'postcheck_round697.py').read_text('utf8'),mapping))
with (HERE/'round698_drafts/research_note_698_draft.md').open('xb') as f:f.write((HERE/'research_note_698.md').read_bytes())
review='''698 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed697 and executed698 entry: progress; navigation/results read.
No live Python found at entry; sequential research sessions completed normally.
Original697 finite centre candidate, all16 channels, actual N256 and original phase preserved.
Frobenius determinant majorant is an inequality for the source, not an absolute-value replacement measure.
Actual T support pairs give five4x4 real blocks. Gaussian radial extraction integrates both completeS9.
Generating coefficients are nonnegative in a,b,C, so componentwise majorants are legitimate.
Beta-square inequality certified by exact positive Bernstein polynomial after removing roots.
Physical odd upper uses convex chord and rational sqrt5 upper56/25.
Color diagonal bound retains odd trace and trace square; cross upper follows support triangle and Young4/5.
Weak cross square1/80 and diagonal1/4 are exact bounds from original even-odd support.
Finite compact sphere polynomial has81 positive rational coefficients, checked against the full20x20 formula.
Weyl full rank4 integration retained.44^4 root grid exact by coordinate degrees43,43,43,42.
Weyl chamber reduction counts each nonzero orbit12 times; zeros on repeated eigenvalues are exact.
Common denominator and strict positive integer bound ensure CRT unique recovery,11 trial-division primes.
Independent48-grid/new-prime check passes; int64 multiplication and sums bounded before use.
Even physical average53/536870912 computed independently by integer Laurent products, not assumed.
Only B11 upper is certified. B10>=727/10^16 remains unproved;697 MonteCarlo cannot replace it.
Conditional vector(1,-3/10000) negative upper<-3.24e-20 computed in Fraction, not a completed RP counterexample.
First broad bounds and all exploratory data kept;157-degree entry unchanged, no repeat round credit.
Space382-386/425/522-523 inherited; no restored Lipschitz,386/425 alternatives,E is not s.
Two new verification groups,18 equations,no images; next699 actual offdiagonal lower bound, goal unchanged.
'''
for name in ('research_note_698.md','joint_static_haar_majorant.py','joint_static_haar_majorant_results.json','unified_physics_condition_ledger_698.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round698_drafts/final_review.txt',review)
print('698 prepared; protected evidence',protected)
