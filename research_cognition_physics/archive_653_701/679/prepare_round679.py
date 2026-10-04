"""Prepare679 with immutable historical evidence and explicit old-space reuse."""
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

ledger=(HERE/'unified_physics_condition_ledger_678.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：678完整来源的局部辅助表示与共同归一',
    '# 联合条件总账：679原物理来源反射与共同Hermitian型')
ledger=ledger.replace('接[677全账](unified_physics_condition_ledger_677.md)，回填[678报告](research_note_678.md)。[结果](joint_local_source_lift_results.json)、[核验](research_round_678_checks.json)。',
    '接[678全账](unified_physics_condition_ledger_678.md)，回填[679报告](research_note_679.md)。[结果](joint_physical_source_reflection_results.json)、[核验](research_round_679_checks.json)。')
ledger=ledger.replace('678全部来源与原质量的有限范围体表示及归一抵消|',
    '678全部来源与原质量的有限范围体表示及归一抵消；679全部物理来源反射及完整Gauss平均的Hermitian型|')
ledger=ledger.replace('678保留全部来源的局部辅助表示及相位；',
    '678保留全部来源的局部辅助表示及相位；679有限软调节及原极限的全部物理来源反射已证；')
ledger=ledger.replace('任意费米来源反射、动态正性、非零物理归一与连续局域性仍缺',
    '679全部原物理Weyl来源反射、奇异延伸及完整物理型Hermitian已证；动态正性、非零物理归一与连续局域性仍缺')
ledger=ledger.replace('一般Grassmann来源合同仍须另外证明。',
    '该轮留下的一般物理Grassmann来源合同由679方程项及接触抵消补齐；不等同任意裸镜像变量反射。')
ledger=ledger.replace('原热参考路线。','原热参考路线。')
ledger=ledger.replace('675精确辅助平均后仍保留完整双Haar；群选择',
    '675精确辅助平均后仍保留完整双Haar；679完整物理来源型Hermitian；群选择')
ledger=ledger.replace('总权重正性和指定零点仍待排除；',
    '679物理全多项式反射型Hermitian仍不证明正性或指定归一非零；')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 679的反射观测与原观测之差为SA；对应来源接触项精确抵消，所以所有物理高斯矩一致。条件只需原配对与gamma5-Hermiticity，不需有限软投影幂等。
- 增广Pfaffian的多项式延伸包含奇异核和零权重；原质量是同一物理代数中的有限指数，不需逆质量。裸N_lambda及裸Phi的简单反射身份一般不成立，661支撑障碍保留。
- 原H_b／双Haar／S9的完整平均继承673、674，故原物理半区型Hermitian。没有证明RP、非零归一、原H_F时间或共同连续极限。
- 依用户提醒再次回查382、383、384、386、425、522、523，逐项复用表见679 §1。已完成的旧合同不重新列为空白；实际端点、参考、仪器与当前物质过程的同一实现仍需接通。

## 本轮合并与下一项

679连接C01、C14、C16、C17、C19、C20、C22：原物理来源、原质量、完整规范平均和共同Hermitian型相容。一般动态规范RP、指定非零归一、原H_F时间与记录、共同连续极限及量子GR仍未完成。原有限图、连续物质、经典几何、手征辅助候选四分支继续分别记账。

接[680](round680_drafts/STATUS.md)：回查668—669，将半区质量插入准确接到完整Gauss候选和物理观测代数，检验可逆性及RP等价范围。若仅重述已有指数公式不新增轮次，不以质量调参、有限Gram样本或第五方向正性代替实际物理正性。目标不改，认知设计后置。
'''
assert '679原物理来源反射' in ledger
write('unified_physics_condition_ledger_679.md',ledger)
mapping={'joint_local_source_lift':'joint_physical_source_reflection',
    '677':'678','678':'679','679':'680','3261':'3263','3263':'3265',
    '1343':'1346','1346':'1349','2473':'2485','2485':'2501'}
verify=replace((HERE/'verify_round678.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_rational_physical_limit as prior_model',
    'import joint_local_source_lift as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_679.md',
       'round679_drafts/research_note_679_draft.md','round679_drafts/final_review.txt',
       'round680_drafts/STATUS.md']
names += ['round679_drafts/'+s for s in (
    'physical_reflection_probe.py','physical_reflection_probe_results.json',
    'physical_reflection_entry.md','check_entry.py','entry_checks.json',
    'reflection_algebra_probe.py','reflection_algebra_probe_results.json',
    'source_shift_probe.py','source_shift_probe_results.json')]
assert len(names)==13
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
verify=verify.replace('Verify complete local source lift, normalization and preserved evidence.',
    'Verify physical-source reflection, singular extension and preserved evidence.')
write('verify_round679.py',verify)
publish=replace((HERE/'publish_round678.py').read_text('utf8'),
    mapping|{'## 324.':'## 325.','## 229.':'## 230.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第679轮完成：** [原物理来源的反射、接触抵消与共同内积]({p}research_note_679.md)'
         '方程项变换与来源接触抵消给全部原物理Weyl来源反射，含有限软调节、奇异质量及零权重；'
         '完整Gauss物理型为Hermitian。两组、十六式通过，最新679／3265，1349份编号科学文件、2501份保护证据。'
         '[核验]({p}research_round_679_checks.json)、[全条件账]({p}unified_physics_condition_ledger_679.md)。'
         '正性、非零归一、原H_F时间及共同连续／量子GR仍开放；旧空间接口逐项复用。')
order=('**当前执行顺序（679后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[680完整质量插入与正性约化]({p}round680_drafts/STATUS.md)，'
       '回查已有半区插入，核完整Gauss物理代数的可逆合同；目标不改。')
"""+publish[b:]
publish=publish.replace('完整物理来源的局部体表示与共同归一','原物理来源的反射、接触抵消与共同内积')
publish=publish.replace('完整局部表示与实际物理过程的区别','来源反射闭合与旧空间条件的逐项复用')
publish=publish.replace('局部辅助体不增加空间维数或物理主体','不恢复旧下界已消去条件，不叠加替代坐标路线')
write('publish_round679.py',publish)
write('postcheck_round679.py',replace((HERE/'postcheck_round678.py').read_text('utf8'),mapping))
write('round679_drafts/research_note_679_draft.md',(HERE/'research_note_679.md').read_text('utf8'))
review='''679 primary-agent proof/code/scope review; no independent agent.
Read current navigation,678,679 entry/results; no running Python.
Re-read382/383/384/386/425/522/523. Reuse explicit table, no dimension reproof.
No Lipschitz reinstatement;386 and425 alternative;523 joint realization retained.
Arbitrary Hermitian epsilon algebra does not claim Wilson realization.
Same original M unitary/skew/chiral, Mbar=Mdagger=-Mstar.
Phi simplification uses Pv=Qminus+gamma5D, no idempotence.
Bare N0 reflection sign and Berezin determinant1 checked.
Tilde Phi differs from Phi star; S independent of D and source.
Contact S F^T-F S^T-S A S^T=0 is exact block identity.
Explicit common block expression in note7 reviewed, source reversal retained.
Gaussian inverse used only on dense set. D=tI gives detN0=t^(2n).
Polynomial extension includes singular N0, variable ranks, zero weights.
Full mass is original physical polynomial exponential; conjugate minus sign.
Not claiming bare Nlambda congruence or time support of Phi.
All source coefficients follow augmented Pfaffian, not low-order sampling.
Original full group, variable phi/E and original16channels retained.
0/2/4/6 source probes: relative only resolved, absolute diagnostics for singular.
Full average reflection follows673/674 invariance and669/677 integrability.
Physical half-algebra form Hermitian, NOT PSD or strictly normalized.
Finite soft included;677 fixed finite box positive tau limit reused, not new theorem.
No uniform Wilson gap, geometry derivative or original HF/time identity claimed.
No new cognitive axiom or revised task/goal; design deferred.
Next680 audits invertible mass insertion in complete Gauss algebra,
not repetition of668 elementary perturbation formula or mass tuning.
2 groups,16 equations. All historical files and drafts preserved.
'''
for name in ('research_note_679.md','joint_physical_source_reflection.py',
             'joint_physical_source_reflection_results.json','unified_physics_condition_ledger_679.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round679_drafts/final_review.txt',review)
print('679 preparation completed')
