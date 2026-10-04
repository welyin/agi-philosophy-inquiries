"""Prepare704 publication without modifying frozen scientific evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k)
          for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

names=('unified_physics_condition_ledger_704.md','round704_drafts/research_note_704_draft.md',
       'round704_drafts/final_review.txt','round705_drafts/STATUS.md',
       'round704_drafts/own_gibbs_compression_entry.py','round704_drafts/own_gibbs_compression_entry_results.json',
       'round704_drafts/own_gibbs_compression_entry.md','round704_drafts/check_and_publish_entry.py',
       'round704_drafts/entry_checks.json','round704_drafts/prepare_entry_publication.py',
       'round704_drafts/literature_scope_audit.json','round704_drafts/initial_preparation_diagnostic.py',
       'round704_drafts/initial_preparation_diagnostic_results.json')
protected=2839+3+len(names);assert protected==2855
ledger=(HERE/'unified_physics_condition_ledger_703.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：703共同几何来源、热响应与原记录',
                      '# 联合条件总账：704准备几何与真实记录历史的共同响应')
ledger=ledger.replace('接[702全账](unified_physics_condition_ledger_702.md)，回填[703报告](research_note_703.md)。[结果](joint_geometry_thermal_limit_results.json)、[核验](research_round_703_checks.json)。',
    '接[703全账](unified_physics_condition_ledger_703.md)，回填[704报告](research_note_704.md)。[结果](joint_preparation_history_limit_results.json)、[核验](research_round_704_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**704当前增量：** 改用625已声明的Galerkin／CP近似族，其自身准备Gibbs态、真实演化和实际记录的总阶数二来源系数共同收敛。能量夹权准备态导数与625弱二阶形式接通，包含准备／动力混合项；这不是703对数转移的动态导数证明。固定图有限过程接口已关闭，空间局域近似、共同连续和动态量子反馈仍开放。\n\n## 当前共同对象及仍存在的分支')
updates={
    'C01 量子对象':'704同一625有限CP族同时保自身几何准备态与真实记录来源',
    'C03 事件记录':'704原CP仪器含尾项，准备／动力混合系数与实际概率归一共同收敛',
    'C04 内部演化':'704同一Galerkin Hamiltonian承担准备及随后真实演化，原Hessian保留',
    'C09 参考系统':'704随准备背景变化的自身Gibbs态有A夹权迹范数的0／1／2阶导数极限',
    'C19 参考态':'704自身有限Gibbs态的準备变化与动态来源共同匹配，实际热化机制仍未导出',
    'C20 尺度映射':'704固定图同一谱族的六类来源系数共同极限；空间局域区域合同另核',
    'C22 来源反作用':'704准备／动力混合项直接由原同一族确定，不能固定初态后忽略'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|原固定图完整量子过程|'):
        cells=row.split('|');cells[2]+='；704自身准备几何与真实历史的共同二阶系数';rows[i]='|'.join(cells)
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            cells=row.split('|');cells[2]+='；'+value;rows[i]='|'.join(cells);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 704自身准备态与原真实历史共用同一来源极限

- 固定A=H₀+c≥1与P_R，h_R=P_RH(x)P_R，初态为该有限族自己的Gibbs态；随后执行真实动态及625原CP记录，尾项不删。
- 623原算符／形式相对界给小系数邻域的Neumann统一控制。热迹界由有限矩阵Lie乘积与Schatten Holder证明，不用指数算符单调性。
- 实热态紧性、复热时及参数Vitali/Cauchy、半热时因子分解共同给Aρ_R^(j)A的迹范数收敛。j=0,1,2；夹权不冒充Tr(A²|ρ′|)。
- 625已经证明真实记录动力来源的共同弱二阶形式。本轮将其A拉回成统一有界、弱算符收敛的代表，与准备态迹类导数配对，关闭总阶数二的混合响应缺口。
- 任意原C²系数路径用链式法则，不新增解析历史公理。固定来源菜单、时间窗与记录次数，不承诺全部源函数空间的一致余项。
- 256维原中性诊断混合CTP总系数约19.44649600i；冻结准备会遗漏该项。相同两支的概率导数仍精确归一；复响应不当概率。
- 224维夹权二阶准备误差仍约.08217，明确不能从态值误差小推断所有导数同样精确。数值只核有限诊断，原全图定理由解析承担。
- 本族与703正转移对数族的有限对象不同；不声称取得后者的真实来源导数，也不改699反例。
- 空间382—386、425、522—524全部按既有范围复用。384已删Lipschitz，386／425替代，523共同实现及524局域UV控制不重新列缺；原空间参考与h／s／CAR／Gauss身份仍待接。

## 本轮合并与下一项

C01／C03／C04／C09／C19／C20／C22在原完整固定图中共用一个自身热准备、真实记录及二阶来源过程。数学谱截止不是物理空间分辨率，尾重置不是免费装置，背景及参数输入未消失。

接[705](round705_drafts/STATUS.md)：回查617／637—644／652／667，核保实际区域操作的局部近似与本共同过程能否接通。精确切边、全谱近似和真实空间细化分别记账；不再优化本轮热导数或重做区域熵存在性。认知设计后置，目标不变。
'''
ledger=ledger.replace('準备','准备')
write('unified_physics_condition_ledger_704.md',ledger)
mapping={'joint_geometry_thermal_limit':'joint_preparation_history_limit',
         '702':'703','703':'704','704':'705','3311':'3313','3313':'3315',
         '1418':'1421','1421':'1424','2823':'2839','2839':'2855'}
verify=replace((HERE/'verify_round703.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_gibbs_preserving_transfer as prior_model',
                      'import joint_geometry_thermal_limit as prior_model')
verify=verify.replace('original geometric thermal-source jets, Hessian contacts and actual records.',
                      'joint preparation-Gibbs and actual dynamical record jets.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==18")
write('verify_round704.py',verify)
pub=replace((HERE/'publish_round703.py').read_text('utf8'),mapping|{'## 349.':'## 350.','## 254.':'## 255.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第704轮完成：** [准备几何与真实记录历史的共同响应]({p}research_note_704.md)'
         '同一625有限CP族的自身准备Gibbs态、真实演化与记录来源共同收敛到总阶数二，'
         '包括准备／动力混合项；夹权热导数接入原弱响应。'
         '两组、十八式通过，最新704／3315，1424份编号科学文件、2855份保护证据。'
         '[核验]({p}research_round_704_checks.json)、[全条件账]({p}unified_physics_condition_ledger_704.md)。'
         '固定图有限过程接口接通，空间区域／连续及动态量子反馈仍开放；旧空间与目标不变。')
order=('**当前执行顺序（704后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[705局部区域与共同尺度映射]({p}round705_drafts/STATUS.md)，'
       '回查617／637—644等，区分精确切分、全谱逼近与真实局部截断；'
       '保留四分支与699范围，停止热导数和有限谱精度优化。')
'''+pub[b:]
pub=pub.replace('共同几何来源、热响应与原记录','准备几何与真实记录历史的共同响应')
pub=pub.replace('原几何来源及热记录共用一个正极限','准备、演化与记录来源共用一个有限过程')
pub=pub.replace('旧空间接口继承，几何热来源不替代时空生成','旧空间接口继承，固定图混合响应不替代空间映射')
write('publish_round704.py',pub)
write('postcheck_round704.py',replace((HERE/'postcheck_round703.py').read_text('utf8'),mapping))
with (HERE/'round704_drafts/research_note_704_draft.md').open('xb') as f:
    f.write((HERE/'research_note_704.md').read_bytes())
review='''704 primary-agent mathematical/code/scope review; no independent agent.
Original623 common operator/form domain and625 weak second-order record response reused.
Previous completed703 and704 entry constitute progress; no goal status change.
Fixed baseline spectral projectors, constant energy shift; original Gauss and all interactions retained.
Own finite Gibbs preparation is not C_R(original Gibbs) at finite cutoff; old625 caution preserved.
Small analytic finite coefficient neighborhood has both operator and form relative bounds.
Neumann inverse order A K^-1 correct; strong convergence follows bounded V A^-1 compressions.
Real embedded heat strong convergence upgraded to trace norm by uniform A energy tails.
Finite matrix exp trace bound proved by Lie product and Holder; no operator monotonicity.
Complex time Cauchy yields K exp heat; bounded strong inverse multiplier supplies left A.
Half-time factorization gives trace norm A heat A, then normalized preparation derivatives.
Real coefficient Schwarz reflection ensures right half heat holomorphic.
No unjustified Tr(A^2 abs(rho derivative)) needed.
625 weak coefficient forms become uniformly bounded WOT-convergent A-pulled representatives.
Trace norm weighted preparation paired with WOT coefficient gives all six joint coefficients.
Only scalar dynamic response claimed; full dynamic cq second trace differentiability not claimed.
Finite source menus, fixed time/record count; no universal source-function-space Taylor remainder.
Original C2 source paths use finite coefficient chain rule, not new analytic-history axiom.
Actual CP defect retained, finite process preparation not rethermalized during evolution.
Numerical diagnostic uses original702 radial neutral sector, not full graph CAR diagonalization.
Exp jets are derivatives, include factorial mixed terms and actual Hessian.
Independent eigensolver four-point mixed difference confirms nonzero cross term.
Complex off-diagonal history coefficients not probabilities; equal histories normalize all jets.
Weighted derivative error large relative to baseline state error explicitly reported.
Initial left-weight-only diagnostic preserved; final includes stronger sandwiched weights.
No claim for703 log-transfer real-time derivatives;699 fixed-candidate failure unchanged.
382 upper,383 involution,384 removed Lipschitz;386/425 alternative bridges.
522 thermal absolute precision,523 common implementation,524 local UV probes directly reused.
Remaining actual h/s/CAR/Gauss and region/space interface not rebranded cognition axiom.
Two scientific groups and18 equations. Next705 targets region-compatible scales.
'''
for name in ('research_note_704.md','joint_preparation_history_limit.py','joint_preparation_history_limit_results.json',
             'unified_physics_condition_ledger_704.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round704_drafts/final_review.txt',review)
print('704 prepared; protected evidence',protected)

