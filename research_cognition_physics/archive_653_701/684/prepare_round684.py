"""Prepare684 review, ledger and publication scripts without altering old files."""
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
ledger=(HERE/'unified_physics_condition_ledger_683.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：683同时间片读口、共同极限与辅助反射',
    '# 联合条件总账：684辅助复制障碍与原物理子代数')
ledger=ledger.replace('接[682全账](unified_physics_condition_ledger_682.md)，回填[683报告](research_note_683.md)。[结果](joint_time_local_source_interface_results.json)、[核验](research_round_683_checks.json)。',
    '接[683全账](unified_physics_condition_ledger_683.md)，回填[684报告](research_note_684.md)。[结果](joint_auxiliary_physical_positivity_results.json)、[核验](research_round_684_checks.json)。')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '684原完整补偿扩充可保全部物理来源而有辅助负方向，完整辅助RP仅为充分路线；原物理态与辅助过程同一性仍缺|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 684适配原G的层/自旋R保持作用协变，但反线性平方为−I；整体复相位不能修复。
- 对G_p=I_p⊗G及S=A⊗R的全部有限实正交复制族，Θ²=1要求A²=−I、p偶。原自由空间零模式的反射Gram谱严格为同数±a/(4+a²)²。因此该整类直接修补失败，不排除所有反射或所有扩充。
- 两复制加入两对真正的G费米补偿后，完整辅助积分逐背景为1，全部原物理来源、S9、双Haar及H_b平均不变。背景响应也须把补偿因子一并微分。
- 在658／661已证的原自由物理正态上，该完整补偿扩充仍有明确辅助玻色负方向。故某个完整辅助表示的RP不是原物理RP的必要条件；这一具体结论不撤销动态原Q0的必要正性。
- 683的物理来源极限可以承接未来对每个物理测试的渐近下界，但该下界尚未证明，不能只凭辅助负数变小就签收正性。
- 直接物理子代数路线保留；661自由结果及680全部物理质量约化继承，不重复计为新定理。停止同类复制设计，转回原共同泛函和真正物理来源。

## 本轮合并与下一项

684连接C01、C04、C14、C16、C20、C22：具体辅助表示、反射对合、全部来源与真正物理子代数必须分别验收。排除一个完整有限复制族，并删除将该辅助全代数RP当作必要条件的误加要求。

原Q0正性、指定归一、H_F及记录身份、共同连续、量子GR和预测仍开放；四分支不混同。空间382—386、425、522—523按679表继承，不恢复已经消去的条件。

接[685](round685_drafts/STATUS.md)：不再改装同类辅助反射。核原共同极限中的未减除体权重及规范/几何来源，保留直接物理子代数或同来源正过程路线；目标不改。
'''
write('unified_physics_condition_ledger_684.md',ledger)
mapping={'joint_time_local_source_interface':'joint_auxiliary_physical_positivity',
    '682':'683','683':'684','684':'685','3271':'3273','3273':'3275',
    '1358':'1361','1361':'1364','2539':'2555','2555':'2567'}
verify=replace((HERE/'verify_round683.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_bulk_fluctuation_matching as prior_model',
    'import joint_time_local_source_interface as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_684.md','round684_drafts/research_note_684_draft.md',
       'round684_drafts/final_review.txt','round685_drafts/STATUS.md']
names+=['round684_drafts/'+p for p in ('compensation_involution_probe.py',
    'compensation_involution_probe_results.json','compensation_involution_entry.md','check_entry.py','entry_checks.json')]
assert len(names)==9
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round684.py',verify)
pub=replace((HERE/'publish_round683.py').read_text('utf8'),mapping|{'## 329.':'## 330.','## 234.':'## 235.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第684轮完成：** [辅助反射的有限复制障碍与原物理子代数]({p}research_note_684.md)'
         '排除原补偿的一整类有限实复制反射；完整补偿可精确保全部物理来源，同时保留辅助负方向。'
         '完整辅助RP仅为充分路线。两组、十六式通过，最新684／3275，1364份编号科学文件、2567份保护证据。'
         '[核验]({p}research_round_684_checks.json)、[全条件账]({p}unified_physics_condition_ledger_684.md)。'
         '原动态Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（684后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[685原共同极限中的体权重与物理来源]({p}round685_drafts/STATUS.md)，'
       '停止同类复制设计，核真实权重及来源，保留直接物理子代数路线；目标不改。')
"""+pub[b:]
pub=pub.replace('同时间片物理读口、完整来源极限与辅助反射边界','辅助反射的有限复制障碍与原物理子代数')
pub=pub.replace('同时间片读口接回原来源，辅助反射仍需核验','辅助全代数正性并非原物理必要条件')
pub=pub.replace('旧空间合同复用，不把辅助反射障碍扩为物理反例','旧空间合同复用，回到真正物理子代数')
write('publish_round684.py',pub)
write('postcheck_round684.py',replace((HERE/'postcheck_round683.py').read_text('utf8'),mapping))
write('round684_drafts/research_note_684_draft.md',(HERE/'research_note_684.md').read_text('utf8'))
review='''684 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed683 and ran684 entry: progress.
Navigation, published results, postcheck, entry and live process list inspected.
Inherited658/661 free physical positivity and normalization; no fixedE substitute.
Actual R acts covariantly on G but anti-linear square equals minus identity.
Phase repair fails since phase cancels in R R*.
All finite real orthogonal factorized copies classified; A²=-I implies even p and A^T=-A.
Actual free W0 sector gives Q=beta A tensor J, with exact balanced ±beta spectrum.
No arbitrary complex, nonfactorized or background-dependent reflection exclusion.
Two copies completed with two G Grassmann pairs; full spectator integral exactly1.
Pfaffian orientation explicit; independent Pfaffian/determinant numeric cross-check.
All original physical sources and complete averages preserved pointwise, even at zero weight.
Free fullS9 physical state remains positive; declared auxiliary observable is strictly negative.
This proves full positivity of a particular auxiliary representation is not necessary.
No negative dynamic physical Gauss witness and no change of the original candidate.
Reflection restriction and convergence statements are acceptance contracts, not RP proofs.
Source-dependent small auxiliary negatives do not prove asymptotic physical positivity.
Primary JaffeJanssens1607.07126 definitionsII.1/II.7 read; bosons not given fermionic twist.
Old space results and four original branches remain distinct; no new principle.
Next685 tests unsubtracted original body weight in actual677 limit, not copy optimization.
Two groups and16 equations; no figure inspection or goal completion.
'''
for name in ('research_note_684.md','joint_auxiliary_physical_positivity.py',
             'joint_auxiliary_physical_positivity_results.json','unified_physics_condition_ledger_684.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round684_drafts/final_review.txt',review)
print('684 preparation completed')
