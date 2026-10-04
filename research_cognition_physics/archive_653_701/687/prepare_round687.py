"""Prepare687 publication from the latest686 baseline."""
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
ledger=(HERE/'unified_physics_condition_ledger_686.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：686谱体匹配与原完整来源','# 联合条件总账：687完整群holonomy正型与精确积分')
ledger=ledger.replace('接[685全账](unified_physics_condition_ledger_685.md)，回填[686报告](research_note_686.md)。[结果](joint_spectral_bulk_matching_results.json)、[核验](research_round_686_checks.json)。','接[686全账](unified_physics_condition_ledger_686.md)，回填[687报告](research_note_687.md)。[结果](joint_full_holonomy_positive_measure_results.json)、[核验](research_round_687_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 687在615同一一格、AP时间、空间平凡、φ=0部门推广到整个原非阿贝尔商群；全局谱帧是R†的多项式，不选群平方根，不丢Jacobian。
- 实际Pf等于[(1+E·O(g)E)/2]^8，全S9平均为五个cos²本征值的h8/495；原SU5载体约束给统一下界25^-8。此处不需原群属于Spin9固定向量分支。
- 原Weyl权重同一字典给|det((I+R)/2)|²，保留其零点；辅助非零不等于每个配置的物理权重非零。
- 全权重是有限酉表示在正密度中的矩阵系数，故为群正型；完整群Haar是正不变投影，权重至少2^-40。
- 保完整Weyl密度的秩四Laurent常数项用20阶根和两素数CRT精确计算，原全S9／单holonomy群Haar权重为261746352167/17416264183971840；无随机积分误差。
- 群正型及非零Haar不是物理时间RP：原673的多空间点、完整H_b／双边界、质量、物理Y来源仍未被替代。特征Hilbert空间只是证书，未认作原H_F或实际记录。
- 旧空间按679表继承。内部群参数不当坐标，384旧条件不恢复，386与425保持替代路线。

## 本轮合并与下一项

687连接C01、C14、C16、C19、C22：原整个规范群、全辅助测度、Weyl权重、正型内积与群投影已在同一平坦一格部门明确成立。认知动机不替代该部门的物理输入。

原动态Q0、指定动态归一、H_F与实际记录身份、共同连续、量子GR和预测仍开放；四分支不混同。静态群Gram不升级为时间过程。

接[688](round688_drafts/STATUS.md)：用实际两个时间片和独立原链路检验球面转移及holonomy表达，保AP闭合及长度归一；不另造正核，不重复同类辅助调参。目标不改。
'''
write('unified_physics_condition_ledger_687.md',ledger)
mapping={'joint_spectral_bulk_matching':'joint_full_holonomy_positive_measure','685':'686','686':'687','687':'688','3277':'3279','3279':'3281','1367':'1370','1370':'1373','2579':'2586','2586':'2596'}
verify=replace((HERE/'verify_round686.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_body_weight_source_matching as prior_model','import joint_spectral_bulk_matching as prior_model')
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==18")
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_687.md','round687_drafts/research_note_687_draft.md','round687_drafts/final_review.txt','round688_drafts/STATUS.md',
       'round687_drafts/full_holonomy_probe.py','round687_drafts/full_holonomy_probe_results.json','round687_drafts/full_holonomy_entry.md']
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round687.py',verify)
pub=replace((HERE/'publish_round686.py').read_text('utf8'),mapping|{'## 332.':'## 333.','## 237.':'## 238.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第687轮完成：** [完整群holonomy、球面正核与精确积分]({p}research_note_687.md)'
         '原平坦一格部门的全S9积分具有整个原群上的正型证书，完整非阿贝尔Haar权重精确为261746352167/17416264183971840。'
         '两组、十八式通过，最新687／3281，1373份编号科学文件、2596份保护证据。'
         '[核验]({p}research_round_687_checks.json)、[全条件账]({p}unified_physics_condition_ledger_687.md)。'
         '群正型不等于原动态物理RP；H_F身份、共同连续及量子GR仍开放，旧空间复用。')
order=('**当前执行顺序（687后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[688实际两时间片的holonomy与球面转移]({p}round688_drafts/STATUS.md)，'
       '由原矩阵检验时间拼接及长度归一，保独立原链路，不另造正核；目标不改。')
""" +pub[b:]
pub=pub.replace('谱体匹配、全部原来源与完整平均','完整群holonomy、球面正核与精确积分').replace('谱体匹配接回原完整平均','原整个群的平坦部门与精确非零投影').replace('旧空间合同复用，谱匹配不等于空间极限','旧空间合同复用，内部群正核不作时空')
write('publish_round687.py',pub)
write('postcheck_round687.py',replace((HERE/'postcheck_round686.py').read_text('utf8'),mapping))
write('round687_drafts/research_note_687_draft.md',(HERE/'research_note_687.md').read_text('utf8'))
review='''687 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed685/686: progress. Current navigation and686 receipts inspected.
No live Python process found; Get-Process empty result is not a running job.
Read615,657,659,672-675,679-680; old space table inherited.
Primary Kikukawa1710.11618v3 text inspected; no general gauge positivity theorem borrowed.
Same615 one-site AP-time massless/zero-phi diagnostic extended to original full group.
Global frames polynomial in R†; no square-root branch or discarded Jacobian.
Actual Pfaffian identity checked, not absolute-value replacement.
Original real vector representation law retained and tested with noncommuting elements.
Full S9 integral via exact moments; uniform positivity uses odd-dimensional SU5 determinant.
Weyl determinant has physical zero weights; no normalized ratio there.
Positive-type proof uses actual coefficient with positive density and finite unitary representation.
Its Hilbert space is a certificate, not an identified physical-time CAR state.
Entire group Haar positivity uses invariant projector; lower bound2^-40 exact.
S(U3×U2) rank4 Haar includes Weyl density and divisor12, accounts originalZ6.
Laurent degree bounds19,19,19,18 justify N20 constant extraction.
Two primes verified by trial division, roots exact order20, integer overflow bounds checked.
CRT modulus exceeds12*4^24*495; positivity and upper bound1 certify unique integer.
Independent20/23 floating quadratures check value, do not establish exactness.
No complete two-time673 H_b/Gauss physical RP, no general spatial fields or nonzero masses.
18 equations/two groups. Full goal unchanged. Next688 must derive actual temporal transfer.
'''
for name in ('research_note_687.md','joint_full_holonomy_positive_measure.py','joint_full_holonomy_positive_measure_results.json','unified_physics_condition_ledger_687.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round687_drafts/final_review.txt',review)
print('687 preparation completed')

