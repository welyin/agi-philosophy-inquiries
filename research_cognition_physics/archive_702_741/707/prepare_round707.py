"""Prepare707 nonlinear scalar blocking publication; preserve all frozen evidence."""
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

names=('unified_physics_condition_ledger_707.md','round707_drafts/research_note_707_draft.md',
       'round707_drafts/final_review.txt','round708_drafts/STATUS.md',
       'round707_drafts/spatial_map_entry.py','round707_drafts/spatial_map_entry_results.json',
       'round707_drafts/spatial_map_entry.md','round707_drafts/check_and_publish_entry.py',
       'round707_drafts/entry_checks.json','round707_drafts/prepare_entry_publication.py',
       'round707_drafts/literature_scope_audit.json','round707_drafts/initial_geodesic_block_results.json')
protected=2885+3+len(names);assert protected==2900
ledger=(HERE/'unified_physics_condition_ledger_706.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：706保边界局部近似与真实记录来源',
    '# 联合条件总账：707曲目标空间分块与原共同过程')
ledger=ledger.replace('接[705全账](unified_physics_condition_ledger_705.md)，回填[706报告](research_note_706.md)。[结果](joint_local_history_source_limit_results.json)、[核验](research_round_706_checks.json)。',
    '接[706全账](unified_physics_condition_ledger_706.md)，回填[707报告](research_note_707.md)。[结果](joint_geodesic_spatial_block_results.json)、[核验](research_round_707_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**707当前增量：** 原H⁵相邻标量作全域测地中点／相对向量变换，保留相对S⁴、全部CAR和环路，只迹掉规范平凡r，得到真正保Gauss的空间分块通道。原热准备及几何导数共同输送；新中点instrument有完整H精确注能范数。全部Dirac／Majorana平均块共享cosh(r/√6)并有差分耦合，不能直接换成裸中点质量。仍非裸单节点H、嵌套尺度闭合或连续极限。\n\n## 当前共同对象及仍存在的分支')
updates={
'C01 量子对象':'707原H⁵全域酉字典后迹掉规范平凡r，正常CP且保物理Gauss支持',
'C02 区域组合':'707保完整相对S⁴及全部CAR／环路的相邻标量块，不等于独立裸节点',
'C03 事件记录':'707声明新中点instrument可精确下降到粗态；多时历史保原内部记忆',
'C14 规范群':'707共同根作用在相对方向而不在r；保角向群载体使粗Gauss支持成立',
'C17 质量机制':'707原全部Dirac／Majorana在同一块内有共同cosh幅度及差分模式耦合',
'C19 参考态':'707原Gauss热态经实际空间通道输送，保总Z；热表示算符不当等待生成元',
'C20 尺度映射':'707超出640二次分支的完整非线性标量分块；未证明跨层组合或原裸形式闭合',
'C21 主体存储':'707新中点读口原完整H注能范数=hbar²M(1/w1+1/w2)/128；装置实现仍缺',
'C22 来源反作用':'707固定空间字典保准备几何0／1／2阶迹类导数；总归一及所有相互作用不丢'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|原固定图完整量子过程|'):
        c=row.split('|');c[2]+='；707实际非线性标量分块的Gauss态、原热来源、记录预算及全质量字典';rows[i]='|'.join(c)
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            c=row.split('|');c[2]+='；'+value;rows[i]='|'.join(c);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 707非线性空间分块的共同对象

- 617／637切口不改变物理格距；707入口原串联电系数因子4及四路径字符余项仅为范围审计，不计完整轮次，也不外推为完整空间细化反例。
- 原574目标是曲率-1/6的H⁵。642树参考保全物质后，两节点可全域变换为中点m和相对z；原Higgs零点不删除。
- Jacobi场给完整10维拉回度量和J(r)=32[sinh(2r/sqrt6)/(2r/sqrt6)]^4，极坐标还含r^4。标架联络及所有旧树电项保留；未把动能改为平直平均。
- 等节点体积时，标量径向半密度转换生成16/R²+8/[R²sinh²(2r/R)]项。其余动力学不独立，域沿原过程输送。
- 保相对S⁴、全部原CAR和环路，只有r是根规范群的平凡因子。物理Hilbert空间确切分为Hc,phys张量Hr，故偏迹真正保Gauss，不仅保协变密度。
- 同一原H及Gauss热态定义Pi=Tr_r exp(-beta H')，Tr Pi为原总Z；准备几何一二阶迹类导数与空间通道交换。-log Pi/beta只是热表示，未识别为自主时间生成元。
- 新仪器测s(中点)，不是旧节点记录的粗分组。其Klein公式只依赖两端s和F，原链路依赖恰消去，完整非对角电项不注能。
- 原标量双交换子给精确注能范数hbar²M(1/w1+1/w2)/128；638工具给同一粗扰动态／参考的相对熵预算。因果装置及资源回收不自动完成。
- 原Dirac及Majorana对x/sqrtF实线性。两节点完整64模式的对称／差分组合中，平均块为cosh(r/R)q(m)，交叉块为-R sinh(r/R)q_u；不可改成裸中点质量。
- 多时记录用原H'和真实关联c+r演化后迹r，保内部记忆。新中点instrument的动态二阶来源图域尚未证明；706只对其原读口有效。
- 三组代码核原完整五维几何／动能测度、真实中点注能、完整64模式质量。没有完整热谱或空间连续模拟；初次结果保存为草稿。
- 382—386、425、522—524直接复用；定性Lipschitz不恢复，386／425替代，523共同实现及524局域UV探针不重证。内部H⁵、S⁴不当物理空间维数。

## 本轮合并与下一项

C01／C02／C03／C14／C17／C19／C20／C21／C22获得同一实际非线性标量块。仍保无限内部方向与全部CAR，没有成为裸单节点理论；原连续物质、动态几何及手征分支尚未共同闭合。

接[708](round708_drafts/STATUS.md)：核中点以外的幅度资料是否同时负责嵌套组合与原集体质量，以及哪些变量可真正丢弃。不要重复泛化的偏迹反例、固定图精度或Gaussian记忆。四分支、699范围及统一目标不变。
'''
write('unified_physics_condition_ledger_707.md',ledger)
mapping={'joint_local_history_source_limit':'joint_geodesic_spatial_block',
         '705':'706','706':'707','707':'708',
         '3318':'3320','3320':'3323','1427':'1430','1430':'1433','2870':'2885','2885':'2900'}
verify=replace((HERE/'verify_round706.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_regional_charge_compression as prior_model',
                      'import joint_local_history_source_limit as prior_model')
verify=verify.replace('local boundary-preserving histories and joint source jets.',
                      'original nonlinear geodesic scalar blocking, records and mass.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
verify=verify.replace("text['display_formulas']==18","text['display_formulas']==20")
verify=verify.replace('==(2,0,0)','==(3,0,0)').replace('fresh_tests=dict(run=2,','fresh_tests=dict(run=3,')
verify=verify.replace("'round706_drafts/attribution_correction_363.md'","'round705_drafts/attribution_correction_363.md'")
write('verify_round707.py',verify)
pub=replace((HERE/'publish_round706.py').read_text('utf8'),mapping|{'## 352.':'## 353.','## 257.':'## 258.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第707轮完成：** [曲目标空间分块与原共同过程]({p}research_note_707.md)'
         '原H⁵中点／相对向量给保Gauss的真实标量分块，原热来源、中点记录预算与全部质量共同输送；'
         '保内部S⁴及全部CAR，仍有记忆。'
         '三组、二十式通过，最新707／3323，1433份编号科学文件、2900份保护证据。'
         '[核验]({p}research_round_707_checks.json)、[全条件账]({p}unified_physics_condition_ledger_707.md)。'
         '不当裸单节点或连续极限，旧空间与目标不变。')
order=('**当前执行顺序（707后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[708嵌套合并与共同质量资料]({p}round708_drafts/STATUS.md)，'
       '核中点以外的幅度是否同时承担组合及原集体质量；'
       '不重复固定图精度、一般偏迹反例或Gaussian记忆，四分支及699范围保留。')
'''+pub[b:]
pub=pub.replace('保边界局部近似与真实记录来源','曲目标空间分块与原共同过程')
pub=pub.replace('局部区域、真实记录与来源共同验收','原非线性空间分块、Gauss热态与质量共同验收')
pub=pub.replace('旧空间接口继承，局部来源极限不当空间细化','旧空间接口继承，内部曲目标不当物理空间维数')
write('publish_round707.py',pub)
write('postcheck_round707.py',replace((HERE/'postcheck_round706.py').read_text('utf8'),mapping))
with (HERE/'round707_drafts/research_note_707_draft.md').open('xb') as f:f.write((HERE/'research_note_707.md').read_bytes())
review="""707 primary-agent proof/code/scope review; no independent agent.
Previous goal turn made progress: formal706 and707 entry published; no blocker audit.
All required current navigation and recent evidence inspected; no live Python work found.
617/637 serial cut vs isotropic refinement factor4 only entry, not a new scientific theorem.
640/641 quadratic blocking/history,642 full tree dictionary and no-go reused, not repeated.
Original574 H5 curvature-1/6 and complete domain; all32+32CAR and original quotient retained.
Midpoint map global in vector z; polar r>0 removes only zero-measure diagonal, no new wall.
Radial parallel frame equivariant under original G fixing origin; gauge action on r trivial.
Jacobi fields derive full10dim metric with connection and Jacobian32[sinh(2r/R)/(2r/R)]^4.
Polar measure r^4 retained; half-density radial potential verified for equal fixed node weights only.
Original tree electric radial/cross terms and all full interactions remain via exact conjugation.
Gauge projection equalsPc tensorIr, so tracing radial factor preserves physical singlet support.
Coarse space retainsS4, allCAR and loops; not bare one-node theory or finite-output construction.
Source-independent block unitary; original Gibbs partial trace keeps fullZ and C2 trace derivatives.
Mean force literature read primary III.1; no bare-bathZ subtraction, no real-time identification.
No claim logPi is operator-norm differentiable or full thermal spectrum numerically computed.
New declared midpoint instrument is not coarse grouping of original single-node records.
Midpoint singlet formula independent of gauge link; hence all electric derivative contributions vanish.
Endpoint Jacobi derivatives give exact gradient bound and complete-H injection norm; origin saturates.
Full relative entropy budget reuses638 with same original actual thermal reference.
Actual bounded finite histories keep eliminated internal r and initial correlations, not reset bath.
Dynamic second-order source domain for new instrument explicitly remains open.
Original full Dirac and Majorana matrices real-linear in hyperboloid spatial coordinates.
Complete64mode symmetric/difference rotation yields sharedcosh and nonzero contrast blocks.
Three numerical groups check full5D geometry, original record resource, all64mass modes.
Preserved first result version; strengthened final result additionally checks full metric and radial density.
Old382-386,425,522-524 retained, removedLipschitz not restored,386/425 alternatives.
No target-space dimension identified with physical space; no Einstein action or gravity generation.
Twenty display formulas; next708 tests nested composition and common mass amplitude.
"""
for name in ('research_note_707.md','joint_geodesic_spatial_block.py','joint_geodesic_spatial_block_results.json',
             'unified_physics_condition_ledger_707.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round707_drafts/final_review.txt',review)
print('707 prepared; protected evidence',protected)
