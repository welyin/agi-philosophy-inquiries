"""Prepare700 checked publication with frozen prior evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k) for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


names=('unified_physics_condition_ledger_700.md','round700_drafts/research_note_700_draft.md',
    'round700_drafts/final_review.txt','round701_drafts/STATUS.md',
    'round700_drafts/mass_time_scope_entry.py','round700_drafts/mass_time_scope_entry_results.json',
    'round700_drafts/mass_time_scope_entry.md','round700_drafts/check_and_publish_entry.py',
    'round700_drafts/entry_checks.json','round700_drafts/prepare_entry_publication.py',
    'round700_drafts/literature_scope_audit.json')
protected=2781+3+len(names);assert protected==2795
ledger=(HERE/'unified_physics_condition_ledger_699.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：699完整积分的严格负方向与候选边界',
                      '# 联合条件总账：700共同时间、物种与空间局域性')
ledger=ledger.replace('接[698全账](unified_physics_condition_ledger_698.md)，回填[699报告](research_note_699.md)。[结果](joint_offdiagonal_haar_certificate_results.json)、[核验](research_round_699_checks.json)。',
    '接[699全账](unified_physics_condition_ledger_699.md)，回填[700报告](research_note_700.md)。[结果](joint_anisotropic_time_contract_results.json)、[核验](research_round_700_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**700新增合同：** 在明确新声明的m0=1各向异性Wilson路径中，物理速度固定ν=c a_t/a_s。空间Wilson项同步缩小会重现八个轻角点；固定r_s>1/2保一个轻角点时，实际投影一阶空间矩≥r_s/(2ν)，故固定空间格距的时间单独细化不保统一局域性。此自由必要部门检验不判断共同连续轨道的动态RP。\n\n## 当前共同对象及仍存在的分支')
updates={
    'C01 量子对象':'700原手征P在时间单独细化下有精确空间导数r_s/(2ν)，不能无成本沿用统一局部观测字典',
    'C04 内部演化':'700实际低能速度固定ν=c a_t/a_s，单独改变热时标签不代替Wilson时间缩放；共同H_F映射仍缺',
    'C15 物种手征':'700新增路径的轻角点数8/7/4/1由r_s阈值决定；r_s=ν且时间单独细化会回到604八重部门',
    'C16 反常测度':'700实际投影导数证明非一致局域性，未把统一硬谱隙重新作为有限积分定义前提',
    'C20 尺度映射':'700固定r_s>1/2时统一物理一阶核矩M*要求a_t/a_s²≥r_s/(2cM*)；此为明确任务合同，原共同细化仍开放'}
rows=ledger.splitlines()
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):rows[i]=row[:-1]+'；'+value+'|'
ledger='\n'.join(rows)+'\n'
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 700时间细化、物种数与原空间投影的共同条件

- 700入口已用680把699负见证输送到同一基线下的任意有限质量：仅调λ(τ)不修复坏τ的全代数RP。精确保全部物理反射型的正扩充或固定参数正极限同样不成立；不计新完整轮次。
- 613的空间overlap Hamiltonian及轴荷规范化路线已经研究，直接继承，不重复实验。700正式增量是另一个显式各向异性四维Wilson路径，m0=1及原16内部表示保留。
- 低能时空导数比给c=νa_s/a_t；r_s=ν时，固定a_s而a_t趋零最终使全部八空间角点变轻。一般轻角点按r_s区间为8、7、4、1；临界点另有Wilson零谱，不擅自赋值当有隙处方。
- 固定r_s>1/2保一个轻角点时，w=0处实际投影导数范数恰r_s/(2ν)，AP固定β下也按1/a_t发散。由Fourier一阶矩直接证明不能保持时间细化一致的空间指数衰减常数；不只凭610充分谱隙界失效猜测。
- 若明确要求统一物理核一阶矩M*，必须a_t/a_s²≥r_s/(2cM*)。这是指定局域性合同的必要条件，不是认知公理，也不是完整相互作用传播速度定理。
- ν=r_s=1、a_t=a_s/c的原自由共同缩放通过本轮必要检查；原自由R²≥1与低能符号是校准，不算完整动态规范、态／来源或GR证明。
- 数值使用全部16通道的1024维原结构矩阵及AP时间；Fourier有限求和只作校准。旧空间各轮继续继承，384／386／425／523不恢复已消去前提，不重复已完成实现。

## 本轮合并与下一项

C04时间、C15物种、C20尺度和C01／C16实际谱投影的局域性现在由同一ν、r_s路径共同约束。统一模型不能分别选择一个时间极限、一个物种数和另一份空间局域性界。原正H_F、指定连续物质、给定作用经典几何和辅助测度四分支仍需真正连接。

接[701](round701_drafts/STATUS.md)：核同时细化时空时物理观测及共同正性需要的统一误差，复用646、659、677、686，区分固定盒来源极限与实际时空共同极限。不要回到已判定的699积分或700角点重复。认知设计后置，目标不变。
'''
write('unified_physics_condition_ledger_700.md',ledger)
mapping={'joint_offdiagonal_haar_certificate':'joint_anisotropic_time_contract',
    '698':'699','699':'700','700':'701','3303':'3305','3305':'3307',
    '1406':'1409','1409':'1412','2761':'2781','2781':str(protected)}
verify=replace((HERE/'verify_round699.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_static_haar_majorant as prior_model','import joint_offdiagonal_haar_certificate as prior_model')
verify=verify.replace('exact offdiagonal full-Haar and sphere integral and strict negative Gram.',
    'anisotropic time, light species and actual spatial-projector locality contract.')
verify=verify.replace("assert text['display_formulas']==20","assert text['display_formulas']==16")
a=verify.index('    names=');z=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[z:]
write('verify_round700.py',verify)
pub=replace((HERE/'publish_round699.py').read_text('utf8'),mapping|{'## 345.':'## 346.','## 250.':'## 251.'})
a=pub.index('summary=');z=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第700轮完成：** [共同时间、物种与空间局域性]({p}research_note_700.md)'
         '原Wilson的声明各向异性路径中，有限速度固定时间／空间系数比；'
         '时间单独细化会重现轻角点，或失去实际投影的统一空间局域性界。'
         '完整AP时间及原16通道已核；共同细化和动态RP仍开放。'
         '两组、十六式通过，最新700／3307，1412份编号科学文件、2795份保护证据。'
         '[核验]({p}research_round_700_checks.json)、[全条件账]({p}unified_physics_condition_ledger_700.md)。'
         '旧空间接口及统一目标不变。')
order=('**当前执行顺序（700后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[701共同尺度下的物理观测与正性]({p}round701_drafts/STATUS.md)，'
       '核固定盒来源极限到实际时空共同极限的统一量词；'
       '不重复699积分或700角点，四分支继续分别记账。')
'''+pub[z:]
pub=pub.replace('原完整积分的严格负方向与候选边界','共同时间、物种与空间局域性')
pub=pub.replace('完整平均负核与固定候选族的适用边界','共同时间缩放约束物种与实际空间投影')
pub=pub.replace('旧空间合同逐项复用，候选负核不否定空间接口','旧空间合同复用，跨尺度局域性与维数定理分开')
write('publish_round700.py',pub)
write('postcheck_round700.py',replace((HERE/'postcheck_round699.py').read_text('utf8'),mapping))
with (HERE/'round700_drafts/research_note_700_draft.md').open('xb') as f:f.write((HERE/'research_note_700.md').read_bytes())
review='''700 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed699 and executed700 entry: progress. Navigation/latest reports read.
Get-Process returned no Python; no live handle was restarted.
613 CHN and Singh comparison explicitly deduplicated; no new Hamiltonian reconstruction claimed.
New input is anisotropic original four-dimensional Wilson path, nu and rs separate,m0=1 fixed.
Full original gamma and all16 channels retained; free/massless flat sector is only a necessary test.
Small momentum ratio c=nu a_s/a_t independent of residue normalization at each light corner.
Exact corner sign count8/7/4/1 excludes critical Wilson zeros; AP finite gap not mistaken for physical mass.
Fixed rs>1/2 gives w=0 at p*=acos(1-1/rs); exact projector derivative rs/(2nu).
No inference of nonlocality from failure of a sufficient gap bound alone.
AP frequency pi/Nt at fixed beta has exact gap identity and positive scaled-derivative limit.
Bloch-angle argument is volume-uniform/continuous-symbol scope,not every discrete finite box hitting p*.
Absolute Fourier first moment dominates actual derivative; full3D bound implies slice bound.
Uniform exponential constants impossible at fixed spatial spacing; each finite nu may remain local.
Physical first-moment contract gives necessary a_t/a_s² bound; not a universal cognitive axiom/LR theorem.
Joint isotropic refinement is an inherited free-symbol calibration,not interacting RP or continuum proof.
Original full1024 matrix with flat group link matches Fourier block; no different internal toy representation.
No new images/tasks/runtime or restored hard-gap requirement for finite integrability.
384 extra Lipschitz removed;386/425 alternatives;523 joint model inherited; E not s record.
Two verification groups,16 equations; next701 actual joint-scale physical observables/positivity.
'''
for name in ('research_note_700.md','joint_anisotropic_time_contract.py','joint_anisotropic_time_contract_results.json','unified_physics_condition_ledger_700.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round700_drafts/final_review.txt',review)
print('700 prepared; protected evidence',protected)
