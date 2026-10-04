"""Prepare753 original-action constrained background and physical-tangent audit."""
import hashlib,json,re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)

ledger=(HERE/'unified_physics_condition_ledger_752.md').read_text('utf8')
ledger=ledger.replace(ledger.split('\n')[0],'# 联合条件总账：753原参考、可积扰动与共同背景',1)
ledger=ledger.replace(ledger.split('\n')[2],
    '2026-10-04。接[752全账](unified_physics_condition_ledger_752.md)，回填[753报告](research_note_753.md)。[结果](joint_reference_constraint_strata_results.json)、[核验](research_round_753_checks.json)。联合目标保持。',1)
updates={
'C01':'753提供无连续联合稳定子的同作用经典背景候选，未证明量子BRST正性或连续物理态',
'C09':'753原满秩X消去联合稳定子的时空部分；新非阿贝尔颜色族沿旧屏障保局部参考，非全局或实际读取',
'C10':'753完整颜色电磁能进入同一Hamiltonian约束，几何随色场初值改变；原共同主锥作用不改',
'C11':'753原零色背景有满足全一阶约束的不可积颜色切向；另有全约束精确颜色族，不是原切向的补偿',
'C14':'753原SU3自由度可激活三方向非阿贝尔连接／电场，使连续颜色背景稳定子从8维降为0；未取消规范冗余或离散中心',
'C22':'753新增源由同一原颜色场算出且完整反作用，未借外部补偿荷；量子绝对来源仍须同态接入'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 753原作用内的共同背景与切向限制

- 满秩内部四参考只消去保持全部背景的连续时空变换，原零色背景仍有su(3)稳定子。不是度规单独无Killing或全约束正则性的充分证明。
- 原连续颜色切向a1=.31T1、e1=.23T2、其他一阶变化为零，通过全一阶约束，却有二次总荷−17.685980218443，故不能保留该切向补成精确闭合初值。
- 另用原三组颜色方向T1、T2、T4，E逐方向与A对齐，Gauss和颜色动量精确为零而磁场非零。其完整能量修改原共形约束，得到同作用正解及短时经典发展。
- 对全部|λ|≤1，新增Y严格小于0.1，原ψ<3/2及z单调证明保留；同一原物质四参考仍局部满秩。
- λ非零时三组颜色电场的共同中心化子在su(3)中为零，结合旧电弱及物质参考，连续联合背景稳定子为零。全部局部规范自由与有限中心仍保持。
- 此结果改善共同微扰背景，不构成全量子物理态、BRST正性、制备或测量实现，也未证明完整引力约束伴随算子满射。
- 双方向中间候选保留为草稿，不单独计轮；完成版本使用原第三个方向，未加作用项或物种。
- 旧空间382—386、425、522—523，604、649/699及来源规范化的适用范围保持。

## 本轮合并与下一项

C09/C11/C14/C10/C22同属一份原作用约束解，C01有更明确的经典背景入口。新增初值是假说候选资料；没有把无连续稳定子升格成认知定义或必需物理公理。

接[754](round754_drafts/STATUS.md)：先复用731右逆与730—741同态来源，检验新增非平坦颜色背景能否同时承接颜色量子来源及引力约束，保共同参考；不继续枚举任意色场或优化残差。统一目标不改。
'''
write('unified_physics_condition_ledger_753.md',ledger)
write('round753_drafts/research_note_753_draft.md',(HERE/'research_note_753.md').read_text('utf8'))
write('round754_drafts/STATUS.md','''# 第754轮入口：无连续稳定子背景与同一量子来源

接[753](../research_note_753.md)、[条件账](../unified_physics_condition_ledger_753.md)、[范围审计](../round753_drafts/scope_and_dedup_review.md)。

1. 753已在原作用内构造全经典初始约束、局部四参考、无连续联合背景稳定子的颜色族；不是改作用、增物种或新增认知公理。参数族不等于自治制备或时间演化。
2. 零色背景的指定一阶切向仍被二次Gauss排除。正面颜色族有另一份允许切向，不把它说成原反例已被补偿。
3. 复用731电弱/动量/共形右逆和730—741的同态来源。固定新增非零颜色背景后，核颜色来源与所有其它来源是否能够共同补偿；不能把固定零色背景的常数荷限制无条件移来，也不能仅凭稳定子为零宣布全量子约束完成。
4. 比较颜色Gauss的实际核、同一物质来源Ward身份和引力初始约束；需要时用753完整场，而非抽象另造自旋／探针模型。
5. 原图Gauss物理态与连续背景CAR态仍不同；绝对来源规范化、BRST物理态正性、实际记录和共同尺度映射保持未完成。
6. 不继续随机颜色背景、特征值精度、窗口或装置优化。只有实质联合连接才计完整轮次，整理写工作稿。
7. 原空间382—386、425、522—523以及604、649/699保持；目标、定时与应用任务不变，不做图像。
''')
write('round753_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[dict(title='Einstein-Yang-Mills Theory: Gauge Invariant Charges and Linearization Instability',
                  url='https://arxiv.org/html/2105.11744',
                  checked='SectionIII equations36-47 and conclusion: quadratic integral consistency on a closed hypersurface. Other conserved-charge constructions assume a timelike Killing field.',
                  use='Known Taub-type obstruction only. Canonical density proof supplied for original nonstationary matter background; no unverified import of pure-gravity or full matter linearization-stability theorem.')],
    inherited='570 finite EW nonlinear Gauss;572-573 exact continuum initial/reference branch;574/728/731 stabilizers;723-725 actual records;730-741 conditional quantum sources.',
    new='Original full material rank and joint symmetry relation; actual continuum colored tangent obstruction; full nonlinear constrained non-Abelian family retaining material chart and eliminating continuous background isotropy.',
    excluded='Full constraint-adjoint surjectivity; BRST positivity; exact continuous quantum instrument; all linearized modes integrable; gauge redundancy eliminated; emergence of dimension, groups or Einstein action.'
),ensure_ascii=False,indent=2)+'\n')
main=('research_note_753.md','joint_reference_constraint_strata.py','joint_reference_constraint_strata_results.json','unified_physics_condition_ledger_753.md')
write('round753_drafts/final_review.txt',
      'Primary-agent review only. Original action retained; initial color data and geometric response explicitly new. Specified nonintegrable tangent distinct from successful family. Full material rank removes joint spacetime symmetry only. Non-Abelian three-direction family satisfies original complete classical initial constraints; old analytic coordinate proof reused with uniform energy barrier. No full BRST/quantum completion claim. Frozen prior evidence preserved.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in main)+'\n')
mapping={'753':'754','752':'753','751':'752','3443':'3446','3441':'3443',
         '1568':'1571','3530':'3541','3521':'3530','398':'399','303':'304',
         'joint_relational_interaction_cones':'joint_reference_constraint_strata'}
pub=remap((HERE/'publish_round752.py').read_text('utf8'),mapping)
summary='**第753轮完成：** [原参考、可积扰动与共同背景]({p}research_note_753.md)原连续零色背景有不可积的一阶颜色切向；另以原三方向颜色场构造全经典约束、原局部参考且无连续联合稳定子的解族。原作用不改，未签收量子物理态。三组、十八式通过，最新753／3446，1571份编号科学文件、3541份保护证据。[核验]({p}research_round_753_checks.json)、[条件账]({p}unified_physics_condition_ledger_753.md)。'
order='**当前执行顺序（753后，优先于下方历史安排）：** 接[754共同背景与量子来源]({p}round754_drafts/STATUS.md)，复用731右逆和730—741同态来源，核颜色来源、引力约束及参考的共同接入；不继续任意色场扫描。[范围审计]({p}round753_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('关系定位、实际交互与共同光锥','原参考、可积扰动与共同背景').replace('旧空间合同保持，原参考交互与传播相容审计','旧空间合同保持，共同约束背景与原参考')
write('publish_round753.py',pub)
write('postcheck_round753.py',remap((HERE/'postcheck_round752.py').read_text('utf8'),mapping))
ver=remap((HERE/'verify_round752.py').read_text('utf8'),mapping)
start=ver.index('    names=(');end=ver.index('    new=',start)
names=('unified_physics_condition_ledger_753.md','round753_drafts/research_note_753_draft.md',
       'round753_drafts/final_review.txt','round753_drafts/literature_scope_audit.json',
       'round753_drafts/scope_and_dedup_review.md','round754_drafts/STATUS.md',
       'round753_drafts/initial_two_direction_source.txt','round753_drafts/initial_two_direction_results.json')
ver=ver[:start]+'    names='+repr(names)+'\n'+ver[end:]
ver=ver.replace('==(2,0,0)','==(3,0,0)').replace('fresh_tests=dict(run=2','fresh_tests=dict(run=3')
ver=ver.replace("checks['display_formulas']==16","checks['display_formulas']==18")
ver=ver.replace('original_reference_interaction_and_common_cone_contract_checked=True',
                'original_joint_classical_background_and_tangent_scope_checked=True')
ver=ver.replace('nonlinear_common_geometry_realized=False',
                'nonlinear_classical_initial_family_realized=True,quantum_common_geometry_realized=False')
write('verify_round753.py',ver)
print('Prepared753 original-action joint background and tangent audit.')
