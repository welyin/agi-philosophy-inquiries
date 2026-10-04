"""Prepare686 ledger and publication without changing older evidence."""
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
ledger=(HERE/'unified_physics_condition_ledger_685.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：685原体权重、非零支撑与共同来源','# 联合条件总账：686谱体匹配与原完整来源')
ledger=ledger.replace('接[684全账](unified_physics_condition_ledger_684.md)，回填[685报告](research_note_685.md)。[结果](joint_body_weight_source_matching_results.json)、[核验](research_round_685_checks.json)。','接[685全账](unified_physics_condition_ledger_685.md)，回填[686报告](research_note_686.md)。[结果](joint_spectral_bulk_matching_results.json)、[核验](research_round_686_checks.json)。')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|','686新增明确非局部谱匹配，原全部未归一来源在完整平均下收敛；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('动态量子反馈、一般几何求导与共同重整化|','685排除常数／曲率型体匹配，686精确谱匹配保全部固定阶来源及原平均；动态量子反馈、一般几何求导与共同重整化|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 686明示S_bulk=(L+1)log det B+L Tr log(I+|H_a|)，原det K减除后残余R在[1,2^n]；可逆H上原共同极限给R→1，k零模处给2^k。
- 673原双Haar零集为零测度、669热矩、677来源多项式界直接复用，得到全部固定阶来源在完整原H_b／双Haar／S9平均下的加权L1收敛。无需硬谱隙、固定谱秩或La²→0。
- 有限质量来源不改变H和控制测度，故其固定阶导数在紧参数集上一致收敛。一般规范／几何导数不在该结论内。
- 该匹配一般非局部且可有零谱尖点；固定盒2^n不提供共同体积／连续界。正标量不等于RP，指定归一及原H_F身份仍缺。
- 615原非零holonomy族证明只减线性领先项，在La²→1的原允许序列上仍留背景相关常数。精确方案与截断方案分开，不能把后者的额外速率要求强加前者。
- 678有限精确补偿与686共同极限匹配各有明确范围；停止重复辅助复制或层数优化，回到原真正物理反射型。
- 旧空间按679表、685 §1继承；没有恢复384已消去条件，没有把386和425合并为必要条件。

## 本轮合并与下一项

686连接C01、C04、C14、C16、C20、C22：明示谱体匹配后，全部原未归一来源和完整平均共同回到原候选；这关闭该匹配方案的来源连接，没有提供正物理时间。

原Q0正性、指定归一、H_F及实际记录身份、共同连续、量子GR和预测仍开放；四分支不混同。空间382—386、425、522—523的旧定理继续复用。

接[687](round687_drafts/STATUS.md)：复用679／680，回查672—675，优先检验真正物理子代数中完整Gauss／S9平均的正性证书或有误差界的负见证；不重复固定E、固定边界或辅助负方向。目标不改。
'''
write('unified_physics_condition_ledger_686.md',ledger)
mapping={'joint_body_weight_source_matching':'joint_spectral_bulk_matching','684':'685','685':'686','686':'687','3275':'3277','3277':'3279','1364':'1367','1367':'1370','2567':'2579','2579':'2586'}
verify=replace((HERE/'verify_round685.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_auxiliary_physical_positivity as prior_model','import joint_body_weight_source_matching as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_686.md','round686_drafts/research_note_686_draft.md','round686_drafts/final_review.txt','round687_drafts/STATUS.md']
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round686.py',verify)
pub=replace((HERE/'publish_round685.py').read_text('utf8'),mapping|{'## 331.':'## 332.','## 236.':'## 237.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第686轮完成：** [谱体匹配、全部原来源与完整平均]({p}research_note_686.md)'
         '明确非局部谱匹配在原调节极限下保全部固定阶来源和完整平均，无需新增硬谱隙或La²条件。'
         '两组、十六式通过，最新686／3279，1370份编号科学文件、2586份保护证据。'
         '[核验]({p}research_round_686_checks.json)、[全条件账]({p}unified_physics_condition_ledger_686.md)。'
         '原物理RP、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（686后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[687原真正物理反射型]({p}round687_drafts/STATUS.md)，'
       '回查672—675，核完整Gauss／S9物理证书或误差受控反例；停止同类辅助表示调参，目标不改。')
""" +pub[b:]
pub=pub.replace('原体权重、非零物理支撑与匹配边界','谱体匹配、全部原来源与完整平均').replace('原体权重接到实际非零物理支撑','谱体匹配接回原完整平均').replace('旧空间合同复用，匹配原权重与来源','旧空间合同复用，谱匹配不等于空间极限')
write('publish_round686.py',pub)
write('postcheck_round686.py',replace((HERE/'postcheck_round685.py').read_text('utf8'),mapping))
write('round686_drafts/research_note_686_draft.md',(HERE/'research_note_686.md').read_text('utf8'))
review='''686 primary-agent proof/code/scope review; no independent agent.
Previous685 publication and postcheck completed; no frozen files edited.
Uses original K, same physical Phi/masses and fixed Berezin orientation.
Exact spectral bulk removal leaves1<=R<=2^n, not exactly1 at finite regulator.
Original Wilson zero sets handled through673 Haar-null theorem, not discarded samples.
Exact k zero modes give residual2^k; no uniform convergence assertion.
Original joint a→0,aL→infinity retained; no new hard gap or La² restriction.
Full original H_b/double Haar/S9 weighted L1 proof uses669/677 envelope.
All fixed-order physical and finite mass sources included; no generic geometry derivatives.
Normalization only conditional on specified nonzero target, not presumed.
Leading-only subtraction counterexample uses615 actual nonzero window and background-dependent action.
Nonlocal spectral matching declared as added input, potentially nonsmooth at zero.
2^n bound not uniform in volume; auxiliary regulator is not physical spatial continuum.
Numerics keep original256 spin/internal nonflat fixture and0/2/4 sources.
At a=.001953125, relative errors6.7% exceeded5% check; added one next dyadic sample.
Final errors1.64%,1.61%,1.60%; no claim of exact physical accuracy from finite samples.
No full Haar or S9 integration numerically approximated; complete-average theorem analytic.
No positive auxiliary process, physical RP, originalHF identification or GR theorem.
Old space results inherited without restoring removed assumptions.
Next687 returns to full physical reflection form; no repeated auxiliary-copy or regulator optimization.
'''
for name in ('research_note_686.md','joint_spectral_bulk_matching.py','joint_spectral_bulk_matching_results.json','unified_physics_condition_ledger_686.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round686_drafts/final_review.txt',review)
print('686 preparation completed')

