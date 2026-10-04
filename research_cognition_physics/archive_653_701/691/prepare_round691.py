"""Prepare691 from verified690; all prior scientific files remain frozen."""
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
ledger=(HERE/'unified_physics_condition_ledger_690.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：690空间传播、物理来源与守恒电荷合同','# 联合条件总账：691中性酉对称与完整正性的分次约化')
ledger=ledger.replace('接[689全账](unified_physics_condition_ledger_689.md)，回填[690报告](research_note_690.md)。[结果](joint_spatial_charge_dictionary_results.json)、[核验](research_round_690_checks.json)。','接[690全账](unified_physics_condition_ledger_690.md)，回填[691报告](research_note_691.md)。[结果](joint_neutral_symmetry_reduction_results.json)、[核验](research_round_691_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 691原中性对称、完整分离与奇Gauss测试

- λ=0时I0(aν bR)=ων(aν)R0(bR)，R0定义为原泛函限制，保留中性归一、全部辅助16通道、H_b、S9及双Haar，不删T(E)行列。
- 全局中性相位严格保原Y正半区、Gauss和反射，保持完整Q0；在已正的中性反射商空间上得到强连续酉作用和整数电荷分解。完整Q0未正前不能称全Hilbert酉重建。
- 热表示的对称实现固定参考向量，不自动是可观测代数中的左电荷算符，也不是690同片符号乘法；局部电荷和实际s仪器仍缺。
- 中性非零Gram配对使分次重排符号消失，完整RP等价于剩余R0的RP。即使只要求全模型偶代数RP，剩余奇Gauss测试也必须通过，因原中性奇场可作正范数伙伴。
- 原u^c d^c d^c颜色epsilon复合为非零奇Gauss不变量，已核整个原商群；中性乘它为合法偶测试，原自由Wick范数精确相同。不是新增粒子或新增物种。
- 680同一H_b下质量RP合同可与本约化合并；非零质量不自动保中性相位对称，各指定归一仍单列。612／646、四分支和全部旧空间条件保持。
- 691切口入口的成熟局部投影流及Ward源复用，不当作全局电荷算符证明；不再在已正中性部门寻找原动态负性。

## 本轮合并与下一项

C01、C04、C14、C19、C22获得共同的原对称／状态／分次正性约化，C03内部实际操作未完成。真正未决者是保留完整辅助和规范平均的R0，以及原过程／尺度／量子几何身份。

接[692](round692_drafts/STATUS.md)：将奇Gauss复合B及中性偶配对接回675—680完整原边界泛函；保全部来源余子式、E与H_b／Gauss积分，不以固定E或带电Gram替代。目标不改。
'''
write('unified_physics_condition_ledger_691.md',ledger)
mapping={'joint_spatial_charge_dictionary':'joint_neutral_symmetry_reduction','689':'690','690':'691','691':'692',
         '3285':'3287','3287':'3289','1379':'1382','1382':'1385','2618':'2628','2628':'2641'}
verify=replace((HERE/'verify_round690.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_physical_car_projector as prior_model','import joint_spatial_charge_dictionary as prior_model')
verify=verify.replace('Verify691 spatial physical source and full-average conserved charge mismatch.','Verify691 actual neutral symmetry, graded RP reduction and odd Gauss probe.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_691.md','round691_drafts/research_note_691_draft.md','round691_drafts/final_review.txt',
       'round692_drafts/STATUS.md','round691_drafts/projected_cut_current_probe.py','round691_drafts/projected_cut_current_probe_results.json',
       'round691_drafts/projected_cut_current_entry.md','round691_drafts/check_entry.py','round691_drafts/entry_checks.json',
       'round691_drafts/neutral_symmetry_initial_results.json']
verify=verify[:a]+'''    for name,digest in core.read(HERE/'round691_drafts/entry_checks.json')['artifact_hashes'].items():
        assert core.digest(HERE/'round691_drafts'/name)==digest,name
'''+ '    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round691.py',verify)
pub=replace((HERE/'publish_round690.py').read_text('utf8'),mapping|{'## 336.':'## 337.','## 241.':'## 242.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第691轮完成：** [中性相位、完整正性约化与奇Gauss测试]({p}research_note_691.md)'
         '原中性相位在已正反射商上有酉实现；完整无质量RP精确约化到保留全部辅助／规范平均的剩余泛函，偶全代数仍需奇Gauss测试。'
         '两组、十六式通过，最新691／3289，1385份编号科学文件、2641份保护证据。'
         '[核验]({p}research_round_691_checks.json)、[全条件账]({p}unified_physics_condition_ledger_691.md)。'
         '酉对称不冒充局部电荷或仪器；剩余动态RP、H_F、共同尺度及量子GR仍开放。')
order=('**当前执行顺序（691后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[692原奇Gauss复合与完整泛函]({p}round692_drafts/STATUS.md)，'
       '保留全部余子式与原平均，复用675—680，不继续修补已分离的中性部门；目标不改。')
"""+pub[b:]
pub=pub.replace('空间传播、原物理电荷与完整平均','中性相位、完整正性约化与奇Gauss测试').replace('空间传播后的同一物理电荷合同','完整物理正性的中性分离与奇Gauss条件').replace('旧空间合同复用，守恒量对象须共同匹配','旧空间合同复用，酉对称不冒充测量仪器')
write('publish_round691.py',pub)
write('postcheck_round691.py',replace((HERE/'postcheck_round690.py').read_text('utf8'),mapping))
write('round691_drafts/research_note_691_draft.md','')
(HERE/'round691_drafts/research_note_691_draft.md').write_bytes((HERE/'research_note_691.md').read_bytes())
review='''691 primary-agent proof/code/scope review; no independent agent.
Previous goal turn690 plus executed691 entry: progress. Current navigation/results/receipts inspected.
No live Python process found; no duplicate job started.
Read661,670,673,679-680 and historical symmetry/odd-Gauss searches.
690 actual neutral free factor extends to all retained-source polynomials and complete average.
Rest functional defined as originalrestriction; all auxiliary16 channels and normalizations retained.
Globalphase preserves sourcehalf,Gauss,Theta and fullHermitianQ; unitarity only on a positive quotient.
Neutral quotient positive inherited freeRP. Strongcontinuity follows finite charge polynomials and unitarybound.
No fulltwo-time-to-HF reconstruction inferred; translation intertwinement only where defined/reconstructed.
Thermal symmetry implementer fixesOmega; not automatically left observable charge or internalinstrument.
Graded reordering sign vanishes in nonzero neutral matrix elements; tensor/Schur RP equivalence exact.
Even fullRP needs oddrestRP via existing neutralodd positive-norm partner, not newprobe species.
Actual B=epsilon u^c d^c d^c nonzero and invariant under original fullG, integer hypercharge0.
Actual source8-point even partner has same norm as odd6-point baryon, preserves Grassmann order.
MassRP congruence680 combines but does not transport neutralU1 symmetry to Majorana/mixedmass.
No originalQ0 sign conclusion or quantumGR closure. Oldspace constraints and fourbranches unchanged.
Full256 neutral diagonal Wick plus120 offdiagonals checked, charge sector counts derived.
Initial random mixedsources mostlyzero; preserved initialresult and replaced with nonzero balanced sources.
Mature localcut Ward entry preserved and reproved as scoped diagnostic, not physicalcharge proof.
Two groups/16 equations. Next692 real oddGauss source in complete dynamic averaged candidate.
'''
for name in ('research_note_691.md','joint_neutral_symmetry_reduction.py','joint_neutral_symmetry_reduction_results.json','unified_physics_condition_ledger_691.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round691_drafts/final_review.txt',review)
print('691 publication preparation completed')
