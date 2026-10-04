"""Prepare689 publication only; preserve all688 and inherited evidence."""
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
ledger=(HERE/'unified_physics_condition_ledger_688.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：688完整Gauss支撑与CAR边缘态','# 联合条件总账：689原物理Y与CAR投影的共同观测')
ledger=ledger.replace('接[687全账](unified_physics_condition_ledger_687.md)，回填[688报告](research_note_688.md)。[结果](joint_gauss_support_marginal_results.json)、[核验](research_round_688_checks.json)。','接[688全账](unified_physics_condition_ledger_688.md)，回填[689报告](research_note_689.md)。[结果](joint_physical_car_projector_results.json)、[核验](research_round_689_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 689实际观测字典的新增连接

- 原653纯时间框中K=(I−Q)/2、A=(I+Q)/2；实际w/barw收缩C=I−2(I−Q)^−1。等时接触系数由原矩阵确定，不能任意略去。646及679的一般接触边界直接继承。
- 同片原Y多项式O_H=det((I+H)/2)exp[−barw(H−I)(H+I)^−1w]精确匹配Γ(H)插入。表面分母由补子式消去，实际块行列式对K及H+I同时奇异仍正则。
- 多个时间片给原有序词det(I−H_diag Q)，局部群源只作用CAR，不改变辅助群链路。普通同片Grassmann积不当作算符乘法；不同片保持原时间序。
- Haar平均O_H给原规范不变偶Y多项式O_Pi，其原完整条件分支期望恰为688的pN。由此关闭688的这个具体原观测缺口；并非仅配分函数的重建记号。
- 不同时间片重复该投影的期望仍为pN。比较态仍限CAR单独零H参考；未判定完整H_F态或动态Q0。
- 空间非平凡时，不能自动把同一符号仍称为CAR投影。E非s记录，四分支、旧空间、654脱耦范围不变。

## 本轮合并与下一项

C01、C14、C19、C22在原纯时间分支中增加了实际观测、规范支撑和状态的同一字典。继承653的时间转移，不重复计定理；原完整H_b、空间传播和一般动态Q0仍须共同验收。

接[690](round690_drafts/STATUS.md)：把明确物理子代数带回673—675完整泛函，先验空间背景中的同一算符含义，再保全部原平均核正性。停止纯时间链重证，目标及任务不变。
'''
write('unified_physics_condition_ledger_689.md',ledger)
mapping={'joint_gauss_support_marginal':'joint_physical_car_projector','687':'688','688':'689','689':'690',
         '3281':'3283','3283':'3285','1373':'1376','1376':'1379','2596':'2609','2609':'2618'}
verify=replace((HERE/'verify_round688.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_full_holonomy_positive_measure as prior_model','import joint_gauss_support_marginal as prior_model')
verify=verify.replace('Verify689 exact Gauss multiplicities, CAR marginal, and corrected653 attribution.','Verify689 actual Weyl/CAR source and projected observable dictionary.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_689.md','round689_drafts/research_note_689_draft.md','round689_drafts/final_review.txt',
       'round690_drafts/STATUS.md','round689_drafts/physical_car_source_probe.py','round689_drafts/physical_car_source_probe_results.json']
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace(",'round688_drafts/attribution_correction_653.md'",",'round687_drafts/attribution_correction_653.md'")
write('verify_round689.py',verify)
pub=replace((HERE/'publish_round688.py').read_text('utf8'),mapping|{'## 334.':'## 335.','## 239.':'## 240.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第689轮完成：** [原物理Weyl观测与CAR规范投影]({p}research_note_689.md)'
         '由原矩阵固定等时接触项，构造同片物理Y多项式，精确读出688的Gauss支撑概率；保非交换时间序及零权重多项式。'
         '两组、十六式通过，最新689／3285，1379份编号科学文件、2618份保护证据。'
         '[核验]({p}research_round_689_checks.json)、[全条件账]({p}unified_physics_condition_ledger_689.md)。'
         '653时间链直接继承，原空间动态Q0、H_F身份、共同连续与量子GR仍开放；旧空间复用。')
order=('**当前执行顺序（689后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[690原物理投影与完整动态泛函]({p}round690_drafts/STATUS.md)，'
       '先核非平凡空间中的算符身份，再保全部原平均检验；不重做纯时间链，目标不改。')
"""+pub[b:]
pub=pub.replace('完整Gauss支撑与CAR边缘态','原物理Weyl观测与CAR规范投影').replace('原联合Gauss约束与子系统状态','原物理观测与同一Gauss支撑').replace('旧空间合同复用，修正继承成果归属','旧空间合同复用，物理观测字典不替代空间')
write('publish_round689.py',pub)
write('postcheck_round689.py',replace((HERE/'postcheck_round688.py').read_text('utf8'),mapping))
write('round689_drafts/research_note_689_draft.md','')
(HERE/'round689_drafts/research_note_689_draft.md').write_bytes((HERE/'research_note_689.md').read_bytes())
review='''689 primary-agent mathematical/source/code review; no independent agent.
Read688 navigation/results,653/654/646/661/670/679 and relevant prior-note keyword matches.
Existing653 temporal spectrum/transfer/Gauss positivity inherited, not new theorem.
Actual original gamma4-adapted spinor projections produce K=(I-Q)/2 and A=(I+Q)/2.
Plus-exponent ordered covariance is -AK^-1; contact sign checked at identity (Ctt=0,n=1/2).
All Grassmann coefficients follow formal Gaussian identity, not inferred from sample orders.
O_H coefficients polynomial by Cauchy-Binet/Jacobi complementary minors, including detD=0.
Original unnormalized insertion block determinant does not divideK or its weight.
Same-time32-mode sources actCAR only, auxiliary links not incorrectly transformed.
Multiple time insertions reduce algebraically to I-HdiagQ, original noncommuting cyclic word.
Same-time ordinary exterior product not identified with CAR multiplication; ordering explicit.
Local gauge covariance tested with independent actualG elements at each time site.
Haar O_H gives actual originalY singlet projector expectationpN, same rational values as688.
Two distinct time projection insertions analytically satisfy PiF^2=PiF.
Initial entry incorrectly required absolute unnormalized determinant>1; fixed to nonsingular inserted matrix.
Full prior HF has other charge carriers; CAR-only comparison not substituted for it.
No extension to generic spatial Q0/RP, no actuals instrument, no new space or GR theorem.
Old382-386/425/522-523 inherited with removed hypotheses not reintroduced.
Two groups/16 equations. Next690 must test actual spatial/dynamic object identity, not repeat puretime.
'''
for name in ('research_note_689.md','joint_physical_car_projector.py','joint_physical_car_projector_results.json','unified_physics_condition_ledger_689.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round689_drafts/final_review.txt',review)
print('689 publication preparation completed')
