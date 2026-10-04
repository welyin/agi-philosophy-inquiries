"""Prepare685 ledger and publication without overwriting prior evidence."""
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
ledger=(HERE/'unified_physics_condition_ledger_684.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：684辅助复制障碍与原物理子代数','# 联合条件总账：685原体权重、非零支撑与共同来源')
ledger=ledger.replace('接[683全账](unified_physics_condition_ledger_683.md)，回填[684报告](research_note_684.md)。[结果](joint_auxiliary_physical_positivity_results.json)、[核验](research_round_684_checks.json)。','接[684全账](unified_physics_condition_ledger_684.md)，回填[685报告](research_note_685.md)。[结果](joint_body_weight_source_matching_results.json)、[核验](research_round_685_checks.json)。')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 685原K的精确体行列式在677共同极限含Tr X及Tr|H|两部分；归一到aL的领先差不需La²→0或统一硬谱隙。
- 原2×2无自环图Tr X常数；615一格周期诊断有绕回自链路，Tr X随holonomy变化，不能混用。
- 直接复用615全S9非零窗口，证明原未减除体因子影响真正非零物理支撑；有限软系数经多项式极限接通，并非固定E替代。
- 同背景归一费米来源比会消去公共体因子；背景权重和来源响应仍改变。纯体积或曲率及其导数项不能抵消平坦holonomy差。
- 一格族的明确来源导数由闭式证明；没有把677推广为一般规范／几何导数收敛。
- 精确非局部谱匹配留下[1,2^n]的残余，下一项核全部来源和原完整平均；所有非局部匹配、其他理论和原Q0正性未被排除。
- 旧空间逐项继承679及685 §1，384已消去条件不恢复，386与425为不同充分路线；无新增空间维数定理。

## 本轮合并与下一项

685连接C01、C04、C14、C16、C20、C22：辅助体、原非零物理支撑、背景权重与来源必须匹配。615旧积分被复用，新连接是它与678原K及677极限的共同验收。

原Q0正性、指定归一、H_F及记录身份、共同连续、量子GR和预测仍开放；四分支不混同。空间382—386、425、522—523按679表继承。

接[686](round686_drafts/STATUS.md)：精确谱体匹配残余、全部固定阶物理来源与原H_b／双Haar／S9平均；不继续辅助复制，不改目标。
'''
write('unified_physics_condition_ledger_685.md',ledger)
mapping={'joint_auxiliary_physical_positivity':'joint_body_weight_source_matching','683':'684','684':'685','685':'686','3273':'3275','3275':'3277','1361':'1364','1364':'1367','2555':'2567','2567':'2579'}
verify=replace((HERE/'verify_round684.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_time_local_source_interface as prior_model','import joint_auxiliary_physical_positivity as prior_model')
a=verify.index('    names=');b=verify.index('    preserved=',a)
names=['unified_physics_condition_ledger_685.md','round685_drafts/research_note_685_draft.md','round685_drafts/final_review.txt','round686_drafts/STATUS.md']
names+=['round685_drafts/'+p for p in ('body_weight_limit_probe.py','body_weight_limit_probe_results.json','body_weight_limit_entry.md','check_entry.py','entry_checks.json')]
assert len(names)==9
verify=verify[:a]+'    names='+repr(tuple(names))+'\n'+verify[b:]
write('verify_round685.py',verify)
pub=replace((HERE/'publish_round684.py').read_text('utf8'),mapping|{'## 330.':'## 331.','## 235.':'## 236.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+"""summary=('**第685轮完成：** [原体权重、非零物理支撑与匹配边界]({p}research_note_685.md)'
         '复用615全S9积分，证明未减除体因子改变原非零物理支撑，曲率或体积常数匹配不足。'
         '两组、十六式通过，最新685／3277，1367份编号科学文件、2579份保护证据。'
         '[核验]({p}research_round_685_checks.json)、[全条件账]({p}unified_physics_condition_ledger_685.md)。'
         '原Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（685后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[686谱体匹配与原完整物理来源]({p}round686_drafts/STATUS.md)，'
       '核残余界、原完整平均与来源；目标不改。')
""" +pub[b:]
pub=pub.replace('辅助反射的有限复制障碍与原物理子代数','原体权重、非零物理支撑与匹配边界').replace('辅助全代数正性并非原物理必要条件','原体权重接到实际非零物理支撑').replace('旧空间合同复用，回到真正物理子代数','旧空间合同复用，匹配原权重与来源')
write('publish_round685.py',pub)
write('postcheck_round685.py',replace((HERE/'postcheck_round684.py').read_text('utf8'),mapping))
write('round685_drafts/research_note_685_draft.md',(HERE/'research_note_685.md').read_text('utf8'))
review='''685 primary-agent proof/code/scope review; no independent agent.
Latest684 receipt, entry and postcheck read. No concurrent Python processes.
Original determinant checked directly and with exact spectral identity.
Leading normalized limit does not impose La²→0 or hard spectral gap.
Original one-site615 and two-by-two677 models explicitly distinguished.
Self-links in615 make Tr X variable; trace norm is instead constant64.
615 full S9 integral and original Weyl determinant reused; not a new integral proof.
Nonzero target support plus polynomial regulator convergence imply finite-soft nonzero support eventually.
Finite soft S9 not numerically integrated; exact target never mislabeled finite regulator.
Conditional physical-source ratios cancel common factor; background weights do not.
No full Haar normalized-model difference or physical RP verdict asserted.
Curvature/volume-only matching excluded, not all local lattice or nonlocal counterterms.
Explicit one-site derivative not a general source-derivative interchange theorem.
Exploratory absolute-error assertions .06/.25 were too restrictive for different response magnitudes.
Final numeric criterion is stated1% relative error plus decreasing errors; analytic limits unchanged.
Primary FurmanShamir hep-lat/9405004v2 sections1-3,6 read; does not prove original K positivity.
16 equations, two groups; old space results and full goal unchanged.
Next686 exact nonlocal spectral matching and full average, no auxiliary-copy redesign.
'''
for name in ('research_note_685.md','joint_body_weight_source_matching.py','joint_body_weight_source_matching_results.json','unified_physics_condition_ledger_685.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round685_drafts/final_review.txt',review)
print('685 preparation completed')

