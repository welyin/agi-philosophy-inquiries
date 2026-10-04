"""Prepare690 from verified689 without overwriting prior evidence."""
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
ledger=(HERE/'unified_physics_condition_ledger_689.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：689原物理Y与CAR投影的共同观测','# 联合条件总账：690空间传播、物理来源与守恒电荷合同')
ledger=ledger.replace('接[688全账](unified_physics_condition_ledger_688.md)，回填[689报告](research_note_689.md)。[结果](joint_physical_car_projector_results.json)、[核验](research_round_689_checks.json)。','接[689全账](unified_physics_condition_ledger_689.md)，回填[690报告](research_note_690.md)。[结果](joint_spatial_charge_dictionary_results.json)、[核验](research_round_690_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 690空间传播与完整原平均的限制

- 689同片来源多项式保留，但一般来源为det(D_H K−B_H A)。L=K+A可逆片的Qeff=(A−K)L^−1不必是单步移位；保detL和原帧方向，不从重写推正时间过程。
- 原ΛevenV的中性单态受全部G恒等作用；在原λ=0基线，实际物理ν关联与全部规范／标量／E背景无关。670／673分解给I0(fν)=ων(fν)I0(1)，完整H_b、S9和双Haar平均保留。
- 原两点空间／两AP时间的π动量收缩为aJ，a=√5−2。占据均值1/2、两时相关(5−2√5)/2；对全中性四模的电荷／逆电荷，原完整平均因子为[1−(1−a²)sin²(θ/2)]²。
- θ=π给161−72√5，非1。因此原来源、689同片符号的同一守恒CAR电荷识别及非零正规化不能同时成立。限定无中性数改变顶点的匹配分支，不泛称所有原相互作用都保该电荷。
- 这是原规范不变物理观测的完整平均身份，不是固定E或不完整Gauss抽样；相应中性小反射Gram反而严格正，不构成Q0负范数反例。
- 扩大过程、改识别真实守恒电荷、受控尺度匹配仍可能；612固定格距全时障碍、646有限物理窗口直接继承，不重复计证明。
- E不是s，旧空间与四分支保持；没有恢复384已消去输入，没有把386和425叠为必要条件。

## 本轮合并与下一项

C01、C04、C14、C19、C22的同一源／守恒电荷合同进一步受限；C20仍允许已有受控尺度路线。原完整Q0正性、H_F身份、量子几何和全目标未完成。

接[691](round691_drafts/STATUS.md)：对接成熟overlap守恒流和切口来源，核真实电荷、原来源与物理正时间支撑，记录非局部项或新增承载自由度。不要再以行列式改写自动签收同一CAR过程。目标不改。
'''
write('unified_physics_condition_ledger_690.md',ledger)
mapping={'joint_physical_car_projector':'joint_spatial_charge_dictionary','688':'689','689':'690','690':'691',
         '3283':'3285','3285':'3287','1376':'1379','1379':'1382','2609':'2618','2618':'2628'}
verify=replace((HERE/'verify_round689.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_gauss_support_marginal as prior_model','import joint_physical_car_projector as prior_model')
verify=verify.replace('Verify690 actual Weyl/CAR source and projected observable dictionary.','Verify690 spatial physical source and full-average conserved charge mismatch.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_690.md','round690_drafts/research_note_690_draft.md','round690_drafts/final_review.txt',
       'round691_drafts/STATUS.md','round690_drafts/spatial_car_source_probe.py','round690_drafts/spatial_car_source_probe_results.json',
       'round690_drafts/spatial_car_source_probe_initial_results.json']
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round690.py',verify)
pub=replace((HERE/'publish_round689.py').read_text('utf8'),mapping|{'## 335.':'## 336.','## 240.':'## 241.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第690轮完成：** [空间传播、原物理电荷与完整平均]({p}research_note_690.md)'
         '原中性通道给完整平均后的精确电荷配对161−72√5，排除将689同片符号无改动当作同一守恒CAR奇偶；不判定Q0反射负性。'
         '两组、十六式通过，最新690／3287，1382份编号科学文件、2628份保护证据。'
         '[核验]({p}research_round_690_checks.json)、[全条件账]({p}unified_physics_condition_ledger_690.md)。'
         '扩大过程、守恒量重识别及受控尺度路线仍开放；旧空间与全目标不变。')
order=('**当前执行顺序（690后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[691守恒流、时间切口与物理观测]({p}round691_drafts/STATUS.md)，'
       '核真实守恒作用、同一来源与正时间支撑；继承612／646，停止纯时间链重证。')
"""+pub[b:]
pub=pub.replace('原物理Weyl观测与CAR规范投影','空间传播、原物理电荷与完整平均').replace('原物理观测与同一Gauss支撑','空间传播后的同一物理电荷合同').replace('旧空间合同复用，物理观测字典不替代空间','旧空间合同复用，守恒量对象须共同匹配')
write('publish_round690.py',pub)
write('postcheck_round690.py',replace((HERE/'postcheck_round689.py').read_text('utf8'),mapping))
write('round690_drafts/research_note_690_draft.md','')
(HERE/'round690_drafts/research_note_690_draft.md').write_bytes((HERE/'research_note_690.md').read_bytes())
review='''690 primary-agent scope/proof/code review; no independent agent.
Previous goal turn688/689: progress; inspected current navigation,689 proof/results/postcheck.
No live Python process. Get-Process empty result is not a waiting job.
Read612,646,653,661,670-675,679-680; no old spectral/static-RP theorem re-counted.
Same original670 256-dimensional kernel, all16 channels and independent originalG links.
Gamma4-adapted physicalY matches689 via fixed spin basis; AP sign fixed before result save.
General source det(DK-BA) is polynomial; Qeff chart requires squareK and invertibleL only there.
Original actual neutral even-exterior vacuum invariant under everyG; no charged probe substitution.
At lambda0, physical neutral Weyl block is Gaussian/free for all actual backgrounds.
670/673 factorization gives raw I0(fnu)=omega_nu(fnu)I0(1), including charged zero weights.
FullHb/S9/doubleHaar averaging preserved analytically, not approximated by seed samples.
Two-site/two-time exact covarianceaJ, a=sqrt5-2, both spin and both momentum modes retained.
Q(sqrt5) exact arithmetic verifies neutral parity pair161-72sqrt5 and occupation mismatch.
Strong contract explicitly identifies local symbols with conserved neutralCAR charge/parity.
No claim that every projection is conserved or all fullHF interactions preserve neutralnumber.
Small neutral reflected Gram positive: dictionary failure not RP counterexample.
Nonzero normalized candidate cannot hide mismatch in remaining average; zeroZ not an accepted state.
Initial random scale violated old kernel-check gap>0.5; kept old bound and weak nontrivial actual links.
Initial occupation result preserved before adding charge-pair evidence; no frozen result overwritten.
Old space384 removed condition not restored;386/425 alternatives;523 existence inherited.
Four branches/fullgoal unchanged; next691 true current/cut dictionary, no arbitrary CAR refit.
'''
for name in ('research_note_690.md','joint_spatial_charge_dictionary.py','joint_spatial_charge_dictionary_results.json','unified_physics_condition_ledger_690.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round690_drafts/final_review.txt',review)
print('690 publication preparation completed')
