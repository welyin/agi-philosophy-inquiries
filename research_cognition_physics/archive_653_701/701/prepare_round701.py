"""Prepare701 publication; preserve old notes, entry evidence and navigation."""
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


names=('unified_physics_condition_ledger_701.md','round701_drafts/research_note_701_draft.md',
       'round701_drafts/final_review.txt','round702_drafts/STATUS.md',
       'round701_drafts/joint_scale_entry.py','round701_drafts/joint_scale_entry_results.json',
       'round701_drafts/joint_scale_entry.md','round701_drafts/check_and_publish_entry.py',
       'round701_drafts/entry_checks.json','round701_drafts/prepare_entry_publication.py',
       'round701_drafts/literature_scope_audit.json')
protected=2795+3+len(names);assert protected==2809
ledger=(HERE/'unified_physics_condition_ledger_700.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：700共同时间、物种与空间局域性',
                      '# 联合条件总账：701观测分辨率与共同正性极限')
ledger=ledger.replace('接[699全账](unified_physics_condition_ledger_699.md)，回填[700报告](research_note_700.md)。[结果](joint_anisotropic_time_contract_results.json)、[核验](research_round_700_checks.json)。',
    '接[700全账](unified_physics_condition_ledger_700.md)，回填[701报告](research_note_701.md)。[结果](joint_flow_resolution_limit_results.json)、[核验](research_round_701_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**701当前增量：** 681自由Abelian比较分支中，完整磁场及固定有限来源多项式有联合格距／平滑极限；显式全动量误差使反射正性在极限恢复。同族的缩放放大来源却保留严格负极限。该结果要求共同列出观测映射和分辨率量词，不能转作699原手征候选的正性证明。\n\n## 当前共同对象及仍存在的分支')
updates={
    'C01 量子对象':'701固定涂抹磁场的有限多项式在同一自由Gaussian测度中共同收敛；不覆盖未重整化点态复合场',
    'C04 内部演化':'701正时间支撑与观测幅度须随尺度明示；有限流失败可与固定观测的正极限并存',
    'C09 参考系统':'524已给523模型的局域紧支撑探针及UV有限噪声，直接复用；不再次列为空白，当前原h/s映射仍缺',
    'C16 反常测度':'701只判定681自由Abelian比较族的极限RP，699原H_b/Gauss/S9负证书不被销账',
    'C20 尺度映射':'701全模式来源误差由sqrt(2rho/pi)exp(-tau²/2rho)S2与a²(S3/12+TS4/24)+a³S4/(2pi³)共同控制；变化来源可有非零负极限'}
rows=ledger.splitlines()
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):rows[i]=row[:-1]+'；'+value+'|'
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 701共同尺度必须连同观测分辨率验收

- 直接继承681实际自由磁场的完整Gaussian积分。精确热核余项给与空间频率无关的反射型误差；完整磁场张量保留两个横向极化，不只看一个模式。
- 明确的自由空间差分族、连续Hamiltonian时间及Fourier来源映射下，全部允许模式与无限尾部都有控制。S4有限的固定物理来源在a、rho同时趋零时收敛；不要求rho/a²固定。
- 普通和反射配对同时收敛，Wick多项式界接成每个固定有限来源代数的正极限。不是完整非Abelian手征构造，也不把物理时间连续输入写成从AP候选推导。
- 同一个实际磁场模式，F_rho=rho^(-1/4)[X(z sqrt(rho))-X(2z sqrt(rho))]的反射范数趋于严格负的-A(z)；z=1/2时约-0.03755638449。负号由正的二阶导数下界证明。
- 因此不能省略来源量词：固定宏观测试成功不推出全部变化测试统一成功；各有限调节器有负测试也不直接否定固定来源的连续极限。
- 本轮空间为比较模型输入。382—386、425、522—523依旧直接复用；524局域探针和UV有限噪声亦已完成，不重做；当前过程身份未补齐，E不是记录s。

## 本轮合并与下一项

C01观测、C04物理时间、C16反射正性及C20跨尺度误差现在在同一个已有比较模型中联合验收。701入口的辅助对角逼近与本轮物理来源极限严格分开；原正H_F、指定连续物质、给定作用经典几何与手征辅助候选四分支仍未自动统一。

接[702](round702_drafts/STATUS.md)：回到原H_b／Gauss／S9反射误差及来源映射，优先查原H_F正过程的受控替代连接，停止自由Gaussian系数优化。699反例继续有效，目标不变。
'''
write('unified_physics_condition_ledger_701.md',ledger)
mapping={'joint_anisotropic_time_contract':'joint_flow_resolution_limit',
         '699':'700','700':'701','701':'702','3305':'3307','3307':'3309',
         '1409':'1412','1412':'1415','2781':'2795','2795':'2809'}
verify=replace((HERE/'verify_round700.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_offdiagonal_haar_certificate as prior_model',
                      'import joint_anisotropic_time_contract as prior_model')
verify=verify.replace('anisotropic time, light species and actual spatial-projector locality contract.',
                      'joint magnetic-source continuum/RP bound and nonuniform negative limit.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
write('verify_round701.py',verify)
pub=replace((HERE/'publish_round700.py').read_text('utf8'),mapping|{'## 346.':'## 347.','## 251.':'## 252.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第701轮完成：** [观测分辨率与共同正性极限]({p}research_note_701.md)'
         '681自由磁场的全空间模式及有限来源多项式有受控正极限；'
         '同族缩放放大来源仍保留严格负范数，故物理观测映射与尺度量词必须共同验收。'
         '不修复699原手征候选。两组、十六式通过，最新701／3309，1415份编号科学文件、2809份保护证据。'
         '[核验]({p}research_round_701_checks.json)、[全条件账]({p}unified_physics_condition_ledger_701.md)。'
         '旧空间与524局域探针结果直接复用，目标不变。')
order=('**当前执行顺序（701后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[702原候选反射误差与来源映射]({p}round702_drafts/STATUS.md)，'
       '回到完整H_b／Gauss／S9和原H_F的受控连接；'
       '停止自由Gaussian系数优化，保留四分支及699反例的准确范围。')
'''+pub[b:]
pub=pub.replace('共同时间、物种与空间局域性','观测分辨率与共同正性极限')
pub=pub.replace('共同时间缩放约束物种与实际空间投影','物理来源分辨率与正性需要共同取极限')
pub=pub.replace('旧空间合同复用，跨尺度局域性与维数定理分开','旧空间及局域探针继承，来源极限不替代维数前提')
write('publish_round701.py',pub)
write('postcheck_round701.py',replace((HERE/'postcheck_round700.py').read_text('utf8'),mapping))
with (HERE/'round701_drafts/research_note_701_draft.md').open('xb') as f:
    f.write((HERE/'research_note_701.md').read_bytes())
review='''701 primary-agent proof/code/scope review; no independent agent.
Previous goal turn completed700 and executed701 entry: progress. No live Python process found.
524 locality/UV finite probe already solved; no duplicate construction performed.
Object is exactly681 declared free Abelian magnetic comparison branch, not original673 chiral measure.
Heat identity follows rho differentiation of integrated Gaussian covariance; sign entrywise not matrixwise.
All spatial modes and both transverse polarizations are retained; magnetic zero mode identically zero.
S4 sources include smooth spatially supported tests, not just fixed momentum windows.
Added spatial lattice is explicitly free/noncompact; physical time is continuous input.
Source map is declared Fourier restriction; not claimed strictly local at finite lattice spacing.
q-k bound from sine remainder; A(q) derivative norm at most4, segment avoids zero.
All omitted modes have |k|>pi/a; tail bounded by a^3 S4/(2pi^3).
Ordinary and reflection covariances converge; Wick bound handles any fixed finite polynomial family.
No unsmeared stress composite, actual instrument, autonomous control, GR or full interacting claim.
Gaussian convolution expansion has finite second-moment remainder uniform for bounded positive omega.
L'' strictly positive certifies nonzero negative scaled-source limit; tests stay at positive times.
Shrinking/renormalized source differs from fixed physical source, so no contradiction.
Finite sums are calibration, infinite-mode conclusions are analytic; no interval certificate claimed.
699 negative candidate is not called UV artefact; its actual joint-scale map is still unproved.
384 extra Lipschitz not restored;386/425 alternatives;523 and524 inherited; E not s.
Two check groups,16 equations, next702 returns to original complete-source comparison.
'''
for name in ('research_note_701.md','joint_flow_resolution_limit.py','joint_flow_resolution_limit_results.json',
             'unified_physics_condition_ledger_701.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round701_drafts/final_review.txt',review)
print('701 prepared; protected evidence',protected)
