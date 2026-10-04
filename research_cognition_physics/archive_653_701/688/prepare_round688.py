"""Prepare688 with an explicit653 attribution correction; preserve frozen files."""
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
ledger=(HERE/'unified_physics_condition_ledger_687.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：687完整群holonomy正型与精确积分','# 联合条件总账：688完整Gauss支撑与CAR边缘态')
ledger=ledger.replace('接[686全账](unified_physics_condition_ledger_686.md)，回填[687报告](research_note_687.md)。[结果](joint_full_holonomy_positive_measure_results.json)、[核验](research_round_687_checks.json)。','接[687全账](unified_physics_condition_ledger_687.md)，回填[688报告](research_note_688.md)。[结果](joint_gauss_support_marginal_results.json)、[核验](research_round_688_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 688优先归属补正及新支撑

[归属补正](round687_drafts/attribution_correction_653.md)优先于上方继承的687历史描述：完整群正型、实际任意纯时间链、球谐谱及扩展Gauss正投影已在653完成。687的新增证据限定为全局多项式帧、显式统一下界及完整群精确积分；688原时间入口仅作旧结果复算，不新计科学组。

- 新算原32CAR与各球谐层的Gauss重数mℓ=(228248,1511416,4996416,11074440,18658304,25771104,30876104,33661096,34779760)。完整原群Weyl常数项用两素数精确认证；独立恢复687Haar值。
- 辅助自身单态数为偶ℓ取1、奇ℓ取0；原CAR单态数n0=228248。该计数不是基本粒子数或空间维数。
- 653正重建态施加联合Gauss后，CAR自身单态概率pN=n0Σ偶μℓ^N/Σmℓμℓ^N。与限定CAR单独零H参考ΠF/n0的迹距离严格为1−pN，任意有限正整数N都非零。
- 该投影是重建CAR的规范不变偶算符，不需测量E；但原Weyl Y等时复合插入字典尚须证明。辅助E不是实际s记录。
- 原完整H_F含可补偿CAR电荷的其他自由度，不等于CAR单独参考。新差异不是原H_F不等价或Q0负性的证明。
- 654已有连续极限脱耦直接继承；旧空间按679表复用，384额外Lipschitz不恢复，386与425保持替代，523共同实现已存在而当前对象映射仍缺。

## 本轮合并与下一项

688连接C01、C14、C16、C19、C20的条件分支支撑与状态；C22全部原Y来源识别单列，不能用配分函数同一代替观测同一。四分支、原动态Q0、H_F身份、共同连续及量子GR的状态不变。

接[689](round689_drafts/STATUS.md)：回查646／660—661／670／679，核纯时间原物理Weyl全部来源及等时复合观测，明确接触、排序与规范合同。只新增尚未接通的对象映射；不重复时间谱或质量拉回。目标不改。
'''
write('unified_physics_condition_ledger_688.md',ledger)
mapping={'joint_full_holonomy_positive_measure':'joint_gauss_support_marginal','686':'687','687':'688','688':'689','3279':'3281','3281':'3283','1370':'1373','1373':'1376','2586':'2596','2596':'2609'}
verify=replace((HERE/'verify_round687.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_spectral_bulk_matching as prior_model','import joint_full_holonomy_positive_measure as prior_model')
verify=verify.replace("text['display_formulas']==18","text['display_formulas']==16")
verify=verify.replace('c in (10,9)','c in (10,13,9)')
verify=verify.replace('Verify strict-time sources, actual gauge covariance and scoped compensation reflection failure.','Verify688 exact Gauss multiplicities, CAR marginal, and corrected653 attribution.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_688.md','round688_drafts/research_note_688_draft.md','round688_drafts/final_review.txt','round689_drafts/STATUS.md',
       'round688_drafts/temporal_sphere_probe.py','round688_drafts/temporal_sphere_probe_results.json','round688_drafts/temporal_sphere_entry.md',
       'round688_drafts/check_entry.py','round688_drafts/entry_checks.json','round687_drafts/attribution_correction_653.md']
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace("'round689_drafts/STATUS.md'):","'round689_drafts/STATUS.md','round687_drafts/attribution_correction_653.md'):")
write('verify_round688.py',verify)
pub=replace((HERE/'publish_round687.py').read_text('utf8'),mapping|{'## 333.':'## 334.','## 238.':'## 239.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第688轮完成：** [完整Gauss支撑与CAR边缘态]({p}research_note_688.md)'
         '精确计算原各球谐层的Gauss重数，证明联合投影后CAR边缘态在有限时间格数下仍不同于限定的CAR单独参考。'
         '完整群正型与时间拼接继承653；[687归属补正]({p}round687_drafts/attribution_correction_653.md)优先。'
         '两组、十六式通过，最新688／3283，1376份编号科学文件、2609份保护证据。'
         '[核验]({p}research_round_688_checks.json)、[全条件账]({p}unified_physics_condition_ledger_688.md)。'
         '原Y观测、动态Q0、H_F身份、共同连续及量子GR仍开放；旧空间复用。')
order=('**当前执行顺序（688后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[689原Weyl来源与CAR等时观测]({p}round689_drafts/STATUS.md)，'
       '回查已有接触项及物理字典，只核新增对象匹配；目标不改。')
"""+pub[b:]
pub=pub.replace('完整群holonomy、球面正核与精确积分','完整Gauss支撑与CAR边缘态').replace('原整个群的平坦部门与精确非零投影','原联合Gauss约束与子系统状态').replace('旧空间合同复用，内部群正核不作时空','旧空间合同复用，修正继承成果归属')
write('publish_round688.py',pub)
write('postcheck_round688.py',replace((HERE/'postcheck_round687.py').read_text('utf8'),mapping))
write('round688_drafts/research_note_688_draft.md',(HERE/'research_note_688.md').read_text('utf8'))
(HERE/'round688_drafts/research_note_688_draft.md').write_bytes((HERE/'research_note_688.md').read_bytes())
review='''688 primary-agent proof/code/scope review; no independent agent.
Baseline687 navigation/result/hash receipts inspected. No concurrent Python job found.
Read653/654/661 and679 space interface table; user task and goal unchanged.
Discovered653 already proves full group positivity, actual arbitrary-N chain and full spectrum.
Written687 attribution correction preserves frozen687 files and overrides novelty claims.
688 temporal entry retained as inherited-result reproduction, not new scientific group.
New actual F=Lambda(R16 tensor C2) character, no unexplained F tensor F* substitution.
Harmonic character h_l-h_(l-2) uses actual10 real representation.
Full original S(U3xU2) Weyl density/divisor12 retained; all Gauss compensation sectors included.
Degree bounds19,19,19,18 justify exact20-root extraction; CRT exceeds positive multiplicity bound.
Original653 rational spectrum read and verified exactly; recomputation does not rename old theorem.
Two primes/exact-order roots/int64 overflow bounds checked by reused algorithm.
Auxiliary invariant count agrees with analytic two-radial-variable harmonic decomposition.
Exact sum m_l mu_l/2^32 independently recovers687 fullHaar integral.
PiF times jointGauss equals PiF tensor PiK, so CAR singlet block is scalarPiF.
Orthogonal block positivity proves exact trace distance1-pN, not just lower bound.
m1>0 and mu1>0 gives strict difference at every finite positiveN.
654 auxiliary decoupling inherited; no new scale optimization or continuum theorem claimed.
PiF is even gauge-invariant reconstructed CAR operator; physical Y composite dictionary still open.
Conditional CAR-only zeroH reference not replaced by original fullHF (other fields can compensate charges).
E not real s record; four branches unchanged; no general originalQ0 or GR conclusion.
Space382-386/425/522-523 inherited; no removed Lipschitz restored; alternative routes not stacked.
Initial code expected old spectrum key degree, corrected to actual l before successful result save.
Two groups/16 equations; next689 targets original all-source equal-time dictionary, avoiding oldmass route.
'''
for name in ('research_note_688.md','joint_gauss_support_marginal.py','joint_gauss_support_marginal_results.json','unified_physics_condition_ledger_688.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round688_drafts/final_review.txt',review)
print('688 publication preparation completed')
