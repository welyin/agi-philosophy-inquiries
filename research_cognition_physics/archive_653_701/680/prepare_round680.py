"""Prepare680 publication; preserve historical notes and all current entry evidence."""
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

ledger=(HERE/'unified_physics_condition_ledger_679.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：679原物理来源反射与共同Hermitian型',
    '# 联合条件总账：680完整质量的正性约化与严格归一边界')
ledger=ledger.replace('接[678全账](unified_physics_condition_ledger_678.md)，回填[679报告](research_note_679.md)。[结果](joint_physical_source_reflection_results.json)、[核验](research_round_679_checks.json)。',
    '接[679全账](unified_physics_condition_ledger_679.md)，回填[680报告](research_note_680.md)。[结果](joint_mass_reflection_congruence_results.json)、[核验](research_round_680_checks.json)。')
ledger=ledger.replace('679全部物理来源反射及完整Gauss平均的Hermitian型|',
    '679全部物理来源反射及完整Gauss平均的Hermitian型；680完整质量RP与同基线无质量RP等价|')
ledger=ledger.replace('679有限软调节及原极限的全部物理来源反射已证；',
    '679有限软调节及原极限的全部物理来源反射已证；680原半区质量为保Gauss的可逆线性乘法，非保单位代数自同构；')
ledger=ledger.replace('动态正性、非零物理归一与连续局域性仍缺|',
    '680固定原H_b下全部允许有限局部质量的RP只需判定原无质量完整对象；该动态正性、指定非零归一与连续局域性仍缺|')
ledger=ledger.replace('677原有限质量来源导数与完整平均的调节极限相容；动态过程',
    '677原有限质量来源导数与完整平均的调节极限相容；680完整物理代数对质量指数及其逆封闭，停止质量调参寻求全代数正性；动态过程')
ledger=ledger.replace('679物理全多项式反射型Hermitian仍不证明正性或指定归一非零；',
    '679物理全多项式反射型Hermitian；680精确外代数反例表明RP、Z0=1及质量乘法可逆仍不保证指定Zlambda非零；原模型指定归一仍需单独证明；')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 680逐项核清原物理时间半区、局部质量Gauss不变性、正反质量指数对多项式增长代数的封闭；完整S9／双Haar／原H_b平均全部保留。因此全多项式RP在固定共同基线下与无质量RP双向等价。
- Rlambda是可逆线性乘法，不保单位／乘法，不自动输送物理时间；不能由此宣布质量无物理效应或各质量模型等于原H_F。
- 有限精确Grassmann例给全部实lambda的RP与可逆质量乘法，同时Zlambda=(1-lambda)^2。该例只反驳一般逻辑蕴涵，不是原16通道模型反例。指定归一继续作为独立验收。
- 正标量权重、低次Gram变好或质量调参均不能替代全代数正性；原Q0完整平均的正负性仍开放。672固定边界负方向不得未经全平均升级。
- 668一般指数插入及669热矩直接复用，不计新定理；本轮新接的是完整动态Gauss物理代数、双向约化及归一独立性。旧空间各项仍按679复用表。

## 本轮合并与下一项

680合并C01、C14、C16、C17、C19：在固定原H_b及共同背景定义下，允许的全部有限局部质量不再需要分别检验RP，只需判定原无质量完整泛函；非零归一仍独立。改变质量同时改变H_b的物理参数族需要另行匹配，不能冒用此固定基线等价。

无质量完整动态规范RP、原H_F／时间与记录身份、受控共同连续极限、量子GR与现实预测仍未完成。四个物理分支保持分别记账；不新增认知公理、不改变目标。

接[681](round681_drafts/STATUS.md)：核成熟动态规范手征构造能否提供同一无质量物理泛函、全部来源与正时间过程；新文献仅按实际维数、背景、流与边界范围接入。停止质量调参或重复小Gram，不用替代模型的成功冒充原对象身份。
'''
assert '680完整质量' in ledger
write('unified_physics_condition_ledger_680.md',ledger)
mapping={'joint_physical_source_reflection':'joint_mass_reflection_congruence',
    '678':'679','679':'680','680':'681','3263':'3265','3265':'3267',
    '1346':'1349','1349':'1352','2485':'2501','2501':'2513'}
verify=replace((HERE/'verify_round679.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_local_source_lift as prior_model',
    'import joint_physical_source_reflection as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_680.md',
       'round680_drafts/research_note_680_draft.md','round680_drafts/final_review.txt',
       'round681_drafts/STATUS.md']
names+=['round680_drafts/'+s for s in ('mass_congruence_probe.py',
    'mass_congruence_probe_results.json','mass_congruence_entry.md','check_entry.py','entry_checks.json')]
assert len(names)==9
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==14")
verify=verify.replace('Verify physical-source reflection, singular extension and preserved evidence.',
    'Verify full original mass congruence and exact normalization boundary.')
write('verify_round680.py',verify)
publish=replace((HERE/'publish_round679.py').read_text('utf8'),
    mapping|{'## 325.':'## 326.','## 230.':'## 231.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第680轮完成：** [完整物理代数中的质量约化与严格归一边界]({p}research_note_680.md)'
         '固定原H_b下，完整Gauss物理代数的质量插入可逆，全部允许局部质量的RP等价于无质量RP；'
         '精确Grassmann反例保留指定归一的独立门槛。两组、十四式通过，最新680／3267，1352份编号科学文件、2513份保护证据。'
         '[核验]({p}research_round_680_checks.json)、[全条件账]({p}unified_physics_condition_ledger_680.md)。'
         '原无质量动态正性、指定归一、H_F身份及共同连续／量子GR仍开放；旧空间接口复用。')
order=('**当前执行顺序（680后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[681原无质量泛函与动态规范构造]({p}round681_drafts/STATUS.md)，'
       '核成熟构造的实际背景、来源与正时间映射，不再质量调参；目标不改。')
"""+publish[b:]
publish=publish.replace('原物理来源的反射、接触抵消与共同内积','完整物理代数中的质量约化与严格归一边界')
publish=publish.replace('来源反射闭合与旧空间条件的逐项复用','质量约化不替代原过程与空间接口')
publish=publish.replace('不恢复旧下界已消去条件，不叠加替代坐标路线','旧方向与热参考合同保持，归一门槛独立')
write('publish_round680.py',publish)
write('postcheck_round680.py',replace((HERE/'postcheck_round679.py').read_text('utf8'),mapping))
write('round680_drafts/research_note_680_draft.md',(HERE/'research_note_680.md').read_text('utf8'))
review='''680 primary-agent code/proof/scope review; no independent agent.
Previous turn679 published and680 entry executed: progress.
Current nav/state/679/results/680 entry inspected; no running Python.
668 forward half-mass perturbation and669 growth moments inherited, not new theorem.
673 full dynamic doubleHaar, original Hb, S9 and physical Y retained.
Mass split is on physical half-algebra, not nonlocal bare Phi support.
Original local mass Gauss scalar; both signs preserve half support and polynomial growth.
Finite Grassmann exponential gives exact inverse without physical matrix inverses.
Complete unnormalized RP congruence and equivalence holds at fixed Hb/other parameters.
Changing mass with simultaneous Hb changes is outside the same-baseline implication.
Rlambda linear invertible, not unital algebra automorphism or original time intertwiner.
No claim that mass is unobservable, models physically identical, or original HF matched.
Full negative witness transported by inverse, not necessarily retained at low source degree.
Normalization Z=Q0(R1,R1) needs separate non-null condition.
Exact rational four-generator Berezin example has all-real-lambda rank-one PSD Gram,
Z0=1, invertible mass multiplier and Z1=0; not original-model or cognition counterexample.
Indefinite exact control: positive scalar weights but transported norm remains minus2.
264 original mass terms,132 per half, all included and then removed; none deleted.
Rank-two Pf and covariance updates checked against direct full augmented Pf for0/2/4/6.
Numerical inverse restricted to resolved fixture; general proof polynomial only.
No full Haar/S9 integral sampled or numerical RP conclusion.
Finite scope, no continuum/geometry derivatives or quantumGR completion.
382-386/425/522-523 inherit679 table, no reinstated conditions or repeated proofs.
Next681 actual mature dynamic gauge/source/time mapping, no new app/automation/goal.
2 groups14 equations; all prior science/drafts retained.
'''
for name in ('research_note_680.md','joint_mass_reflection_congruence.py',
             'joint_mass_reflection_congruence_results.json','unified_physics_condition_ledger_680.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round680_drafts/final_review.txt',review)
print('680 preparation completed')
