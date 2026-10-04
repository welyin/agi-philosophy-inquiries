"""Prepare677 publication without modifying frozen science."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_676.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：676非平坦边界与原泛函零支撑',
    '# 联合条件总账：677全谱调节与完整物理来源的共同极限')
ledger=ledger.replace('接[675全账](unified_physics_condition_ledger_675.md)，回填[676报告](research_note_676.md)。[结果](joint_nonflat_seed_support_results.json)、[核验](research_round_676_checks.json)。',
    '接[676全账](unified_physics_condition_ledger_676.md)，回填[677报告](research_note_677.md)。[结果](joint_rational_physical_limit_results.json)、[核验](research_round_677_checks.json)。')
ledger=ledger.replace('676指定非平坦边界表示的秩障碍与原零支撑|',
    '676指定非平坦边界表示的秩障碍与原零支撑；677全谱调节下原物理来源与完整平均的共同极限|')
ledger=ledger.replace('原物理态与辅助过程同一性仍缺|',
    '677完整未归一来源在固定有限盒与正热时收敛；原物理态与辅助过程同一性仍缺|')
ledger=ledger.replace('动态过程、机制选择及现实参数|',
    '677原有限质量来源导数与完整平均的调节极限相容；动态过程、机制选择及现实参数|')
ledger=ledger.replace('676固定Jminus边界像不能覆盖部分非平坦目标谱，须保留原泛函零支撑；',
    '676固定Jminus边界像不能覆盖部分非平坦目标谱；677改用全谱调节，已在固定有限盒完整平均中保留原泛函零支撑，不要求统一谱隙；')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 677把659已有符号共同极限接到673全固定尺寸来源及全部原平均；不是新的符号近似方法。669热矩与673两时间边界零集Haar零测度直接复用。
- 加权L1／有限质量来源导数成立于固定盒、固定正原H_b热时；不含规范／几何来源求导、物理连续极限或归一。
- 有限软投影不是精确手征投影，不能直接使用673三角解耦、675因子分離或676有限矩阵零模计数。原零支撑在极限恢复。
- 空间382—386、425、522—523逐项复用；384已消去的Lipschitz不重新列为缺口，386／425两路线不要求同时满足。

## 本轮合并与下一项

677连接C01、C14、C16、C17、C20、C26：原观测、质量、规范输送和完整辅助／配置／双Haar平均，在同一有界全谱表示下保留。没有更换原16通道、原字段或相位；当前两组数值只验证实际字典和来源，不替代完整积分的解析证明。

剩余关键条件是同一局部体表示的测度与来源身份、一般动态规范正性、非零归一、原H_F物理时间及四分支的共同连续极限。第五方向与物理时间仍区分，原Gibbs到实际仪器／空间端点的映射仍须落实。

接[678](round678_drafts/STATUS.md)：用原局部Wilson算子检验全谱块递推，核整个质量／来源泛函和体行列式。线性解算身份不能冒充反射正过程；成熟domain-wall方法可借用，具体原测度必须逐项匹配。目标与旧空间合同保持。
'''
write('unified_physics_condition_ledger_677.md',ledger)
mapping={'joint_nonflat_seed_support':'joint_rational_physical_limit',
         '675':'676','676':'677','677':'678','3257':'3259','3259':'3261',
         '1337':'1340','1340':'1343','2444':'2460','2460':'2473'}
verify=replace((HERE/'verify_round676.py').read_text('utf8'),mapping)
verify=verify.replace('Verify exact nonflat seed obstruction and original physical null support.',
    'Verify full-spectrum physical-source limit and preserved historical scope.')
verify=verify.replace('import joint_physical_boundary_average as prior_model',
    'import joint_nonflat_seed_support as prior_model')
a=verify.index("    section=result[");b=verify.index("    assert (result['tests_run']",a)
verify=verify[:a]+verify[b:]
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_677.md','round677_drafts/research_note_677_draft.md',
           'round677_drafts/final_review.txt','round678_drafts/STATUS.md',
           'round677_drafts/rational_support_probe.py','round677_drafts/rational_support_probe_results.json',
           'round677_drafts/rational_limit_probe.py','round677_drafts/rational_limit_probe_results.json',
           'round677_drafts/rational_support_entry.md','round677_drafts/entry_checks.json')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==12","assert text['display_formulas']==16")
write('verify_round677.py',verify)
publish=replace((HERE/'publish_round676.py').read_text('utf8'),mapping|{'## 322.':'## 323.','## 227.':'## 228.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第677轮完成：** [全谱有理表示与完整物理来源平均的共同极限]({p}research_note_677.md)'
         '全谱软调节连接原固定尺寸观测、质量来源与全部辅助／配置／双Haar平均；'
         '固定有限盒、正热时的加权L1极限无需统一Wilson谱隙，并保留原零支撑。'
         '两组、十六式通过，最新677／3261，1343份编号科学文件、2473份保护证据。'
         '[核验]({p}research_round_677_checks.json)、[全条件账]({p}unified_physics_condition_ledger_677.md)。'
         '局部体身份、正性、原物理时间、连续及量子GR仍开放；旧空间接口直接复用。')
order=('**当前执行顺序（677后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[678全谱局部块表示与共同来源]({p}round678_drafts/STATUS.md)，'
       '核原质量／来源及体减除，不把第五方向当物理时间；目标不改。')
'''+publish[b:]
publish=publish.replace('非平坦边界初始空间与原物理泛函的零权重支撑','全谱有理表示与完整物理来源平均的共同极限')
publish=publish.replace('非平坦边界表示与原物理来源支撑','完整平均的受控极限与旧空间接口复用')
publish=publish.replace('固定初始空间失效与物理边缘反例的区别','有限调节极限不替代实际空间与物理过程')
write('publish_round677.py',publish)
write('postcheck_round677.py',replace((HERE/'postcheck_round676.py').read_text('utf8'),mapping))
write('round677_drafts/research_note_677_draft.md',(HERE/'research_note_677.md').read_text('utf8'))
review='''677 primary-agent proof/code/scope review; no independent agent.
Prior676 completed and published;677 entry executed. Progress.
Read navigation/state/prior note/results and old382-386,425,522-523 interfaces.
Removed384 Lipschitz not resurrected.386 and425 are alternative bridges.
523 existing joint thermal/actual-direction realization preserved.
Full-spectrum tanh regulator is mature, reused659, not a new physics theory.
No fixed Jminus/rank used.0<a<1 guarantees positive fifth transfer.
Wilson h<=2d-1 derives from unitary shift and spin projectors, not space dimension.
Resolvent perturbation proof bounds sign difference without simultaneous eigenbasis.
Pointwise joint a->0,aL->infinity, not uniform gap over full configuration.
Soft P are positive contractions, not exact idempotents. No triangular identity
or exact finite nullity imported from673/675/676 at finite soft regulator.
Phi bound golden ratio. Difference bound sqrt5/4 times epsilon difference.
N skew, full source augmented Pf bounded by original field polynomial.
673 Haar-null Hzero for every spatial config;669 all original heat moments.
Fubini + DCT yields weighted L1 of every source coefficient and full average.
Mass derivatives via finite polynomial coefficients, compact parameter sets only.
No arbitrary moving-Wilson geometry/gauge source derivative, joint tau/volume limit.
Finite gauge covariance retains both fields and source inverse transpose,
unit determinant doubled Berezin map, no absolute Pf or positive projection.
676 zero physical support preserved in limit; no relative ratio to exact zero.
Numerical balanced original background has checked invertible N; derivative
inverse only numerical check, not theorem assumption. Last weight error6.7%.
L millions is spectral evaluation, not millions of constructed local layers.
No full Haar numerical calculation, finite RP, normalized state, original HF,
physical time, continuum/quantum GR or newly derived spatial dimension.
678 must audit sparse local block identity, original sources and bulk factor.
2 new numeric groups,16 displayed equations.
'''
for name in ('research_note_677.md','joint_rational_physical_limit.py','joint_rational_physical_limit_results.json','unified_physics_condition_ledger_677.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round677_drafts/final_review.txt',review)
print('677 preparation completed')

