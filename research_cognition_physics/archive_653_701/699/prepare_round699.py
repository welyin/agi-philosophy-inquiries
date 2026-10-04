"""Prepare checked 699 publication; keep all previous evidence immutable."""
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


names=('unified_physics_condition_ledger_699.md','round699_drafts/research_note_699_draft.md',
    'round699_drafts/final_review.txt','round700_drafts/STATUS.md',
    'round699_drafts/free_domination_entry.py','round699_drafts/free_domination_entry_results.json',
    'round699_drafts/free_domination_entry.md','round699_drafts/check_and_publish_entry.py',
    'round699_drafts/entry_checks.json','round699_drafts/prepare_entry_publication.py',
    'round699_drafts/clifford_moment_probe.py','round699_drafts/clifford_moment_probe_results.json',
    'round699_drafts/clifford_moment_certificate.py',
    'round699_drafts/offdiagonal_exact_sphere.py','round699_drafts/offdiagonal_exact_sphere_results.json',
    'round699_drafts/offdiagonal_exact_group.py','round699_drafts/offdiagonal_exact_group_results.json')
protected=2761+3+len(names);assert protected==2781
ledger=(HERE/'unified_physics_condition_ledger_698.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：698原静态中心项的完整积分上界',
                      '# 联合条件总账：699完整积分的严格负方向与候选边界')
ledger=ledger.replace('接[697全账](unified_physics_condition_ledger_697.md)，回填[698报告](research_note_698.md)。[结果](joint_static_haar_majorant_results.json)、[核验](research_round_698_checks.json)。',
    '接[698全账](unified_physics_condition_ledger_698.md)，回填[699报告](research_note_699.md)。[结果](joint_offdiagonal_haar_certificate_results.json)、[核验](research_round_699_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**699当前判定优先：** 全部原群和两个S9的标量边缘已有严格负方向；675将其接成固定候选族的充分大有限τ反例。它关闭“该固定overlap权重与原H_b对全部τ都正”的接法，不关闭原正H_F、另行共同缩放的路线或连续统一目标。下方旧轮次过程记录中的“仍待证”保留其历史语义，当前状态以本段和699结尾为准。\n\n## 当前共同对象及仍存在的分支')
ledger=ledger.replace('一般动态规范正性、指定归一、原全图物理时间／CAR态与连续极限|',
    '699已排除此固定无质量候选对所有τ的动态RP；共同时间／参数路径、指定归一、原全图物理时间／CAR态与连续极限仍须另接|')
updates={
    'C01 量子对象':'699在真实E无关物理标量子代数获得完整积分负方向，故此候选族不能对全部τ给正态；不否定原H_F态',
    'C04 内部演化':'699固定overlap／质量参数而增大原H_b热时的族失败；物理共同时间轨道尚未识别，不能把τ随意等同原H_F时间',
    'C16 反常测度':'699精确B10与698严格B11上界认证全平均负核；675给充分大有限τ失败，不宣称全部τ或连续极限失败',
    'C19 参考态':'699原H_b正基态和谱隙直接复用，负二次型不能由该固定缝合在大τ获得正物理参考；未计算数值τ阈值'}
rows=ledger.splitlines()
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):rows[i]=row[:-1]+'；'+value+'|'
ledger='\n'.join(rows)+'\n'
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 699原完整平均后的负核已经严格认证

- 697原256维带相位来源、全16通道及中心历史不变。原Clifford代数中精确四次湮灭恒等式给32维辅助行列式等于95项多项式e4的四次方；解析相位在h=I颜色端点固定，不替换为绝对值。
- 两个完整S9用Wick／Dirichlet矩积掉，最终45项Q(√2)完全齐次多项式。自由Clifford极限精确校准；h=I结果精确回到696独立值，旧B00另经同一算法校准。
- 原全商群Weyl常数项次数(19,19,19,18)，20阶格点及6个素数唯一恢复两个整数；24阶及新素数核算。B10精确值约7.313249396801595×10⁻¹⁴，有理下界严格大于727/10¹⁶。
- 继承原B00和698的|B11|上界，向量(1,−3/10000)的完整二次型≤−2.9193954230389916×10⁻¹⁹；负号由有理运算认证，不依赖697抽样。
- 直接应用675，保原完整H_b的固定参数候选在充分大的有限τ下不满足物理RP。热核阈值未数值求出；不外推到全部τ、原H_F、共同连续极限或认知公理。单自环不替代真实空间传播。
- C07—C09的旧接口按699 §1.1逐项继承，没有新空间条件；384额外Lipschitz不恢复，386／425不叠加，523共同实现不再列为空白。E不是s记录。

## 本轮合并与下一项

C01、C04、C16、C19从“原完整积分符号未判定”推进到一条限定明确的候选反例：固定overlap辅助权重和原H_b任意热时不能自动组成共同正过程。原正H_F、指定连续物质、给定作用经典几何及辅助候选四分支继续区分。

接[700](round700_drafts/STATUS.md)：回查已有共同时间与来源字典，核对候选失败究竟限制哪条物理参数路径，寻找可保全部物理观测的替代连接。停止同一已判定积分的优化；先合并共同条件，再讨论底层认知设计。统一目标不变。
'''
write('unified_physics_condition_ledger_699.md',ledger)
mapping={'joint_static_haar_majorant':'joint_offdiagonal_haar_certificate',
    '697':'698','698':'699','699':'700','3301':'3303','3303':'3305',
    '1403':'1406','1406':'1409','2738':'2761','2761':str(protected)}
verify=replace((HERE/'verify_round698.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_full_holonomy_reduction as prior_model','import joint_static_haar_majorant as prior_model')
verify=verify.replace("assert prior_model.algebra_checks()==core.read(prior_model.TARGET)['algebra']",
    'assert prior_model.run()==core.read(prior_model.TARGET)')
verify=verify.replace('exact full-Haar and sphere upper bound for the original static centre entry.',
    'exact offdiagonal full-Haar and sphere integral and strict negative Gram.')
verify=verify.replace("assert text['display_formulas']==18","assert text['display_formulas']==20")
a=verify.index('    names=');z=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[z:]
write('verify_round699.py',verify)
pub=replace((HERE/'publish_round698.py').read_text('utf8'),mapping|{'## 344.':'## 345.','## 249.':'## 250.'})
a=pub.index('summary=');z=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第699轮完成：** [原完整积分的严格负方向与候选边界]({p}research_note_699.md)'
         '原B10精确求出，结合698上界，完整群及两个S9平均后有严格负二次型。'
         '由675可知固定候选与原H_b在充分大有限τ下RP失败；不否定原H_F或统一目标。'
         '两组、二十式通过，最新699／3305，1409份编号科学文件、2781份保护证据。'
         '[核验]({p}research_round_699_checks.json)、[全条件账]({p}unified_physics_condition_ledger_699.md)。'
         '旧空间接口逐项复用，四分支继续区分。')
order=('**当前执行顺序（699后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[700候选的共同时间与替代接口]({p}round700_drafts/STATUS.md)，'
       '停止优化已判定积分，核物理时间与参数路径、同一来源的替代映射；'
       '不把有限候选失败外推为全部时间、连续极限或认知原则失败。')
'''+pub[z:]
pub=pub.replace('原静态中心项的完整积分严格上界','原完整积分的严格负方向与候选边界')
pub=pub.replace('原静态项的完整积分上界及剩余非对角条件','完整平均负核与固定候选族的适用边界')
pub=pub.replace('旧空间合同复用，积分界不新增空间条件','旧空间合同逐项复用，候选负核不否定空间接口')
write('publish_round699.py',pub)
write('postcheck_round699.py',replace((HERE/'postcheck_round698.py').read_text('utf8'),mapping))
with (HERE/'round699_drafts/research_note_699_draft.md').open('xb') as f:f.write((HERE/'research_note_699.md').read_bytes())
review='''699 primary-agent proof/code/scope review; no independent agent.
Read latest navigation698/entry699, original675/673/591 and user-specified old space reports.
No Python process returned by Get-Process. CIM process inspection was denied; no dependence on it.
Original256 source retains all16 channels, phase, entire quotient Haar and twoS9.
Exact six-vector Clifford algebra plus weakR; R-odd trace zero because at most3 weak vectors.
Exact quartic annihilator, normalized trace recurrence and Newton identities give charpoly Q^8.
Source reality and connected real-analytic domain fix detS'=e4^4 at positive color endpoint.
No absolute-value change to the original measure; positivity only this offdiagonal centre subfamily.
First sphere uses full Wick/radial moments; second sphere Dirichlet1^5 moments;45 terms degree8.
Independent free U calibration exact; identity holonomy recovers696 complete sphere fraction.
Original T Clifford and weak-volume dictionary checked; six original256 matrix comparisons pass.
Complete original Weyl CT degree19,19,19,18;20grid exact; ordered chamber orbit12 cancels divisor.
Qsqrt2 conjugate corresponds to orthogonal U, contractive S and bounded physical factor.
Both real embeddings bounded1 imply integer coefficient boundsD; CRT6 proven primes modulus>2D.
Independent24grid/new-prime check; all int64 operations within proven limits.
Exact B10 lower plus inherited B00 and certified698 U11 imply negative real quadratic form.
Original Hb large-tau corollary uses675,not a new independent round or changed Casimir model.
Finite spatial self-loop preserves591 connected elliptic/confining hypotheses; original positive parameters retained.
Finite tau threshold not computed; all tau failure,originalHF failure,continuum failure not claimed.
No cognition axiom refutation; four branches remain distinct. Old space contracts unchanged.
384 extra Lipschitz removed;386/425 alternatives;522 absolute difference task;523 joint model exists.
Current h,s,CAR/Gauss mapping remains; auxiliary E is not s record. No images/new task/dependencies.
Two verification groups,20 equations; next700 shared time/candidate alternative, goal unchanged.
'''
for name in ('research_note_699.md','joint_offdiagonal_haar_certificate.py','joint_offdiagonal_haar_certificate_results.json','unified_physics_condition_ledger_699.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round699_drafts/final_review.txt',review)
print('699 prepared; protected evidence',protected)
