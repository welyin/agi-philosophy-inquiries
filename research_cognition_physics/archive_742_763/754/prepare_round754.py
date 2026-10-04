"""Prepare754 all-color source/constraint and relative quantum-state bridge."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_753.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：754全颜色来源、内部补偿与量子态差',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[753全账](unified_physics_condition_ledger_753.md)，回填[754报告](research_note_754.md)。[结果](joint_irreducible_source_completion_results.json)、[核验](research_round_754_checks.json)。联合目标保持。',1)
updates={
'C01':'754原背景CAR存在有限秩相差的两份准自由Hadamard态，颜色来源差均值非零；仍非动态全Gauss物理态',
'C09':'754同一配置下的小来源补偿保初片参考jet非退化；实际读取和全局参考未完成',
'C10':'754固定753非平坦颜色背景后，全颜色来源及动量/能量共同接入初始几何；首阶态差响应复用732',
'C11':'754全模颜色Gauss逆与731三向反流/共形法联立，任意光滑外给源在正C合同下可完成全部初始约束',
'C14':'754固定非零颜色连接的八个旧常数核消失；非零分量均值由原规范场承担，不是取消整体Gauss',
'C22':'754颜色补偿的动量及一阶/二阶电能全保；真态差源满足共同Ward，有限强度数值仍是外给诊断'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 754全颜色来源与共同初始约束

- 固定753颜色幅度1的原连接，-D_iD_i在全部Fourier模有正下界：零模用有理LDL证明大于1/100，非零模用解析范数界大于1/4。
- 无需给颜色源施加零坐标均值。其均值由原场的i[A_i,delta E_i]承担，不借外部荷，也不把该坐标分量积分说成脱离参考的物理净荷。
- 与731原电弱逆、Gauss中性三向反流、全部电动量能量及共形约束共同求解；新颜色动量和原非零颜色电场的一阶能量不能遗漏。
- 两份原背景CAR准自由Hadamard态可相差两个光滑u_R模式投影，确有非零颜色总分量。其全部应力/电流/质量力来自同一态差，满足旧732共同Ward。
- 因此该真实来源的短时线性反作用可接入；绝对首阶应用仍以735处方及741共同阶次为条件。
- 数值J、rho为显式诊断，并未计算成该两态的真实连续应力。有限强度初始解不等于有限强度自洽量子反馈。
- 未构造全Gauss/引力量子物理态、自治准备仪器、量子几何正性或图的连续极限；旧空间、604及649/699保持。

## 本轮合并与下一项

C01/C11/C14/C10/C22通过同一背景上的真实相对来源和全初始右逆相连，C09保持小变化范围。减少的是零色背景对量子来源的八个附加均值限制；原维数、群、作用、参考态及参数输入没有因此被推出。

接[755](round755_drafts/STATUS.md)：审查共同物理状态和约束的量子接口，先复用730—741、751与649/699；不继续枚举源函数或优化有限矩阵。统一目标不改。
'''
write('unified_physics_condition_ledger_754.md',ledger)
write('round754_drafts/research_note_754_draft.md',(HERE/'research_note_754.md').read_text('utf8'))
write('round755_drafts/STATUS.md','''# 第755轮入口：共同背景、量子物理状态与实际过程

接[754](../research_note_754.md)、[条件账](../unified_physics_condition_ledger_754.md)、[范围审计](../round754_drafts/scope_and_dedup_review.md)。

1. 753/754给原作用的非平坦颜色共同背景和全初始来源右逆。真Hadamard态差可有非零颜色分量均值，并能接短时线性反作用。
2. 以上是背景CAR与经典平均约束接口，不等于动态规范/引力的全量子物理态。755优先处理这项对象身份，避免继续外给诊断源。
3. 先复用730—741同态参考与有限阶响应、751完整来源/资源、590/623—625正过程和649/699限定失败。成熟BV、Hadamard、Gauss或路径积分结果不能重复计轮。
4. 无连续背景稳定子改善候选背景，但不自动提供BRST正性、全部约束闭合或实际instrument。不得把本轮的源补偿数学右逆叫作自治制备。
5. 需要分支时沿认知观察提出相容假说组合，写清共同原因、独立输入及失败判据，不把所有修补累积成新公理。
6. 不继续任意颜色背景、源样式、低频矩阵或求解精度优化；若只有整理，保存工作报告不虚增完成轮次。
7. 旧空间382—386、425、522—523、604、649/699保持。目标、定时、应用任务不变，无图像。
''')
write('round754_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
      dict(title='A Phase Space Approach to the Conformal Construction of Non-Vacuum Initial Data Sets in General Relativity',
           url='https://arxiv.org/html/2106.15027v2',
           checked='Canonical matter data, CMC source scaling and momentum compatibility; examples do not directly prove full chiral quantum feedback.',
           use='Mature method inherited from731; original color/energy dictionary independently used.'),
      dict(title='The renormalized locally covariant Dirac field',
           url='https://arxiv.org/html/1210.4031',
           checked='Gauge/Yukawa backgrounds, Hadamard existence, Wick fields; stated stress conservation result has a restricted background range.',
           use='Tools only; full original matrix/Majorana and source Ward scope relies on inherited730/732/735, not a blanket extension of the paper.')
    ],
    inherited='731 EW/counterflow/conformal right inverse;732 relative source Ward and local linear response;735/741 absolute conditional ordering;753 actual constrained nonzero color background.',
    new='Fixed color configuration all-mode inverse, removal of old zero-mean source restriction by internal fields, full simultaneous initial completion, actual Hadamard pair source witness.',
    excluded='Numerical diagnostic stress as an actual quantum stress; exact nonlinear semiclassical feedback; autonomous state preparation; full BRST or gravity physical state.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_754.md','joint_irreducible_source_completion.py','joint_irreducible_source_completion_results.json','unified_physics_condition_ledger_754.md')
write('round754_drafts/final_review.txt',
      'Primary-agent review only. Fixed-connection Gauss inverse proved for all modes independently of full stabilizer claim. Exact rational zero-mode certificate; nonzero-mode analytic bound. All induced color momentum and electric energy retained. Quantum finite-rank Hadamard state comparison distinct from numerical prescribed-source diagnostics. No full quantum-gravity or exact feedback claim. Historical evidence preserved.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'754':'755','753':'754','752':'753','3446':'3449','3443':'3446',
         '1571':'1574','3541':'3550','3530':'3541','399':'400','304':'305',
         'joint_reference_constraint_strata':'joint_irreducible_source_completion'}
pub=remap((HERE/'publish_round753.py').read_text('utf8'),mapping)
summary='**第754轮完成：** [全颜色来源与共同约束]({p}research_note_754.md)原753背景的全模颜色Gauss逆接通任意光滑来源及完整动量/能量补偿；真实Hadamard态差可带非零颜色均值，并接条件首阶反作用。数值为外给诊断，未签收全量子物理态。三组、二十式通过，最新754／3449，1574份编号科学文件、3550份保护证据。[核验]({p}research_round_754_checks.json)、[条件账]({p}unified_physics_condition_ledger_754.md)。'
order='**当前执行顺序（754后，优先于下方历史安排）：** 接[755共同物理状态与过程]({p}round755_drafts/STATUS.md)，核背景CAR、平均来源约束与全量子状态的实际连接；复用旧正过程和限定失败，不继续诊断源优化。[范围审计]({p}round754_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('原参考、可积扰动与共同背景','全颜色来源与共同约束').replace('旧空间合同保持，共同约束背景与原参考','旧空间合同保持，量子来源与完整初始约束')
write('publish_round754.py',pub)
write('postcheck_round754.py',remap((HERE/'postcheck_round753.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round753.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_754.md','round754_drafts/research_note_754_draft.md',
       'round754_drafts/final_review.txt','round754_drafts/literature_scope_audit.json',
       'round754_drafts/scope_and_dedup_review.md','round755_drafts/STATUS.md')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace("checks['display_formulas']==18","checks['display_formulas']==20")
ver=ver.replace('original_joint_classical_background_and_tangent_scope_checked=True',
                'all_mode_color_inverse_and_joint_source_scope_checked=True')
ver=ver.replace('nonlinear_classical_initial_family_realized=True,quantum_common_geometry_realized=False',
                'nonlinear_prescribed_initial_constraints_solved=True,self_consistent_quantum_feedback_proven=False')
write('verify_round754.py',ver)
print('Prepared754 all-color joint source completion.')

