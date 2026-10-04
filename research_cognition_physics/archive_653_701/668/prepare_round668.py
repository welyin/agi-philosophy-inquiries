"""Prepare668 evidence and navigation publication without rewriting history."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_667.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：667局部参考、空间传播与尺度接口','# 联合条件总账：668原质量、物理观测与辅助正泛函')
ledger=ledger.replace('接[666全账](unified_physics_condition_ledger_666.md)，回填[667报告](research_note_667.md)。[结果](joint_reference_spatial_process_results.json)、[核验](research_round_667_checks.json)。',
    '接[667全账](unified_physics_condition_ledger_667.md)，回填[668报告](research_note_668.md)。[结果](joint_mass_auxiliary_reflection_results.json)、[核验](research_round_668_checks.json)。')
ledger=ledger.replace('660—661自由有限盒物理观测反射|一般规范背景、原全图物理时间／CAR态、完整质量与连续极限',
    '660—661自由有限盒物理观测；668原常数质量的正泛函及严格归一|原动态标量、一般规范背景、原全图物理时间／CAR态与连续极限')
ledger=ledger.replace('原物理态与原手征辅助泛函的同一完整表示|',
    '668已接指定常数质量；原物理态与辅助泛函的同一完整表示仍缺|')
ledger=ledger.replace('667完整有限记录输送|自主实际关系instrument和动态引力反馈',
    '667完整有限记录输送|668冻结标量不足以承载原s读口；真实共同instrument和动态引力反馈')
ledger=ledger.replace('一般规范相互作用局域测度；参考排序不是连续反常',
    '668已接固定原质量与正归一；动态标量和一般规范相互作用测度仍缺，参考排序不是连续反常')
ledger=ledger.replace('667跳跃重相位同步改变Dirac质量|机制选择、一般背景极点及现实参数',
    '667跳跃／Dirac质量同步；668原全部常数质量的物理Weyl拉回|实际动态质量、机制选择、一般背景极点及现实参数')
ledger=ledger.replace('实际热态尺度极限及宇宙初态；不再要求同一空Fock全局酉',
    '668有限正泛函不等于原完整Gibbs；实际热态尺度极限及宇宙初态，不要求空Fock全局酉')
ledger=ledger.replace('667质量、跳跃及曲目标边势共同来源|动态量子反馈',
    '667质量／跳跃／曲目标边势；668质量与观测映射的同一导数|动态量子反馈')
ledger=ledger.replace('667原条件热迹／响应可区分正确字典和只改跳跃的错误重写',
    '667热迹、668质量权重／来源可排除指定错误改写')
a=ledger.index('## 下一项')
ledger=ledger[:a]+'''## 本轮合并与下一项

668在同一物理Weyl代数中加入实际全部常数Dirac／Majorana质量，C01、C16、C17、C19、C22共同成立于声明的有限自由规范背景。全多项式反射正性及任意实质量强度严格归一已经关闭；不再以“只证自由无质量”为此范围缺口。裸局部镜像字段加质量与实际物理观测加质量不等价，前者保原权重的指定接法已排除。

仍须区分常数φ、时变外部φ与原量子动态φ。新正泛函不等于原完整Gauss Gibbs，不解除612全时谱及653荷谱限制。553单尺度／运行区别保留。221—222、230作为过程／伙伴／有限误差验收，不直接外推无界场。404是667参考表示区别的历史概念依据，369半侧模包含仍需真实共同参考和包含前提。

接[669](round669_drafts/STATUS.md)：使用623原标量热半群和同一质量函数，核动态标量、反射与积分范围，再核624原s读口；不另拟无关Gaussian标量，不把Euclidean函数插入等同实际instrument。旧空间证明及替代路线按上表直接复用。认知设计后置；四分支尚未全部统一，目标不改。
'''
write('unified_physics_condition_ledger_668.md',ledger)
mapping={'joint_reference_spatial_process':'joint_mass_auxiliary_reflection',
    '666':'667','667':'668','668':'669','3239':'3241','3241':'3243',
    '1310':'1313','1313':'1316','2329':'2342','2342':'2354'}
verify=replace((HERE/'verify_round667.py').read_text('utf8'),mapping)
verify=verify.replace('Verify local reference transport, actual spatial masses and inherited geometry scope.',
    'Verify actual original mass reflection, strict finite normalization and common sources.')
verify=verify.replace('import joint_reference_mass_stabilizer as prior_model','import joint_reference_spatial_process as prior_model')
verify=verify.replace('reference_scale_probe','mass_observation_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_668.md','round668_drafts/research_note_668_draft.md',
           'round668_drafts/final_review.txt','round669_drafts/STATUS.md',
           'round668_drafts/mass_observation_entry.md','round668_drafts/mass_observation_probe.py',
           'round668_drafts/mass_observation_probe_results.json',
           'round668_drafts/mass_observation_probe_first_attempt.py',
           'round668_drafts/mass_probe_first_failure.json')
"""+verify[b:]
write('verify_round668.py',verify)
publish=replace((HERE/'publish_round667.py').read_text('utf8'),mapping|{'## 313.':'## 314.','## 218.':'## 219.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第668轮完成：** [原非零质量、物理观测与辅助测度的共同反射结构]({p}research_note_668.md)'
         '原常数质量接入同一物理Weyl观测，有限自由盒有全多项式反射正性和严格归一；来源保留观测映射导数。'
         '两组、十六式通过，最新668／3243，1316份编号科学文件、2354份保护证据。'
         '[核验]({p}research_round_668_checks.json)、[全条件账]({p}unified_physics_condition_ledger_668.md)。'
         '原动态标量、完整Gauss过程、连续及量子GR仍开放；旧空间接口直接继承。')
order=('**当前执行顺序（668后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[669原动态标量、物理质量与记录]({p}round669_drafts/STATUS.md)，'
       '核原标量热过程、反射和积分，再接原s读口；目标不改。')
"""+publish[b:]
publish=publish.replace('局部参考、空间传播与共同量子过程的相容条件','原非零质量、物理观测与辅助测度的共同反射结构')
publish=publish.replace('局部参考与原空间传播的共同连接','原质量与物理观测的共同正泛函')
publish=publish.replace('复用旧维数与坐标证明，核实际物质和尺度映射','固定原质量的反射接口与动态标量缺口')
write('publish_round668.py',publish)
write('postcheck_round668.py',replace((HERE/'postcheck_round667.py').read_text('utf8'),mapping))
write('round668_drafts/research_note_668_draft.md',(HERE/'research_note_668.md').read_text('utf8'))
review='''668 primary-agent proof/code audit; no independent agent review.
Previous goal turn completed667 and actual668 entry: progress. Navigation,
latest note/results/postcheck and live process state read; no Python running.
No goal edits, app tasks, automation, agents or images. Existing Python/NumPy.
Actual614 original mass/exterior dictionary, all16 internal species and complex
Dirac/Majorana couplings. Constant phi and unit gauge finite free box declared.
Euclidean mass prescription is additional input, not full original Gibbs identity.
Physical mass uses Phi.T P Phi, not bare candidate local mass. Work in c,e'
variables without inverse Jminus dagger v; original auxiliary zero modes handled
by polynomial identity rather than illegal ratios. Full phase/detS preserved.
Weyl and original S9 weight factorization follows exact triangular congruence.
RP uses661 original physical algebra and even onsite Vplus+Theta Vplus. It does
not transfer support through nonlocal Phi. Full polynomial proof precedes tests.
Strict partition: actual free symbol, even AP time, no zero denominator. Takagi
diagonalizes constant symmetric mass in physical Weyl sector, not auxiliary U16
symmetry. Momentum block determinant ratio squared, Pfaffian continuous from1;
product positive for all real lambda. Not extended to dynamic scalars or gauge.
Two-point reflection diagnostics are not the proof. Numerical large Pfaffian
complex residue is assessed relatively, not hidden or rounded to exact real.
Common source includes Phi derivatives and mass derivative, all phase retained.
Finite difference comparison and independent Fourier mass increment agree.
Source numerics conditional on declared E; not full bosonic stress integral.
Entry run rerun without counting as new group; rectangular helper development
failure and first code preserved. Formal science passed unchanged thresholds.
Side-chat historical clues read against221,222,230,369,404,553. Do not claim new
general process criterion or algebraic state idea. Fixed-scale identity is not
running spectral relation. Actual observable/process mapping still needed.
382-386,425,522-523 retained: no removed Lipschitz restored; alternative bridges
not stacked. Old theorems are not new science or automatic current realizations.
Two new groups,16 equations. Next669 actual dynamic original scalar and s record
mapping; distinguish Euclidean insertions from instruments. Full goal active.
'''
for name in ('research_note_668.md','joint_mass_auxiliary_reflection.py',
             'joint_mass_auxiliary_reflection_results.json','unified_physics_condition_ledger_668.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round668_drafts/final_review.txt',review)
print('668 report/review frozen; verification and publication prepared')
