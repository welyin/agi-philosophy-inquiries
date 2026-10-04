"""Freeze670 report and prepare audited publication, preserving earlier files."""
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


ledger=(HERE/'unified_physics_condition_ledger_669.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：669原动态标量、手征物质与规范过程',
                      '# 联合条件总账：670非平坦规范、共同观测与无质量零模')
ledger=ledger.replace('接[668全账](unified_physics_condition_ledger_668.md)，回填[669报告](research_note_669.md)。[结果](joint_dynamic_scalar_gauge_interface_results.json)、[核验](research_round_669_checks.json)。',
    '接[669全账](unified_physics_condition_ledger_669.md)，回填[670报告](research_note_670.md)。[结果](joint_nonflat_mass_measure_results.json)、[核验](research_round_670_checks.json)。')
ledger=ledger.replace('668原常数质量的正泛函及严格归一|原动态标量、一般规范背景、原全图物理时间／CAR态与连续极限|',
    '668常质量正泛函；669动态标量可积；670非平坦共同字典跨无质量零模|一般动态规范正性、指定归一、原全图物理时间／CAR态与连续极限|')
ledger=ledger.replace('668已接指定常数质量；原物理态与辅助泛函的同一完整表示仍缺',
    '670原观测／质量有无逆K的正则字典；原物理态与辅助泛函的同一完整表示仍缺')
ledger=ledger.replace('群选择及真实规范动力学仍缺|',
    '670非交换非平坦背景共同协变已核；群选择及真实规范动力学仍缺|')
ledger=ledger.replace('668固定质量已归一；669动态标量候选可积且未归一正，指定归一及一般规范测度仍缺',
    '668固定自由质量已归一；669动态标量未归一正；670质量／来源权重跨K零模。动态规范正性、指定归一与全局测度仍缺')
ledger=ledger.replace('669实际动态质量插入可积；同一规范传播、机制选择、一般背景极点及现实参数',
    '669动态插入可积；670质量与规范传播共同字典相容并有零模提升见证；动态过程、机制选择及现实参数')
ledger=ledger.replace('668质量与观测映射的同一导数|',
    '668共同导数；670非平坦来源及零模处正则求导|')
ledger=ledger.replace('667—668字典检查；669原质量与传播不同步会改变规范相关预测|',
    '667—669同步字典；670非平坦完整相位、物理二点和质量来源检查|')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 670在原615／660手征配对下，消去物理观测Φ及质量拉回中的逆Kₗ。Kₗ无质量零模不再是该有限字典的排除条件；这不移除Wilson H谱隙、固定index、动态正性和整体归一要求。

## 本轮合并与下一项

670合并C01、C14、C16、C17、C22：原投影恒等式使共同观测Φ正则且算子范数不超过黄金比；全部原质量、辅助权重、物理来源与非交换非平坦SM链路协变。原holonomy有H隙不闭、Kₗ零模、加入原质量后权重非零的具体见证，且完整局部式与物理来源式一致。

660无质量权重延伸是旧成果；本轮新增实际观测映射的无逆K身份、含质量的全部来源连接和非平坦共同检查。初次坐标发散预期失败已保存，不以放宽断言掩盖。旧空间接口不重复、384已消去条件不复活，386／425保持替代路线。

动态规范反射正性、全配置／跨index测度、指定整体归一、原完整H_F／记录身份及连续／量子GR仍开放。有限协变不等于正性，光滑谱投影不等于正半区支撑；661旧限制仍有效。605—610热支撑与电动能接口、643原Gauss有序历史全部继承。

接[671](round671_drafts/STATUS.md)：先检验反射兼容的真实空间曲率背景中同一物理泛函，再审查动态规范正时间组合与原过程映射。不得把任意不对称背景Gram问题视为整体反证，也不以有限二点正性宣称动态测度完成。认知设计后置，四分支尚未统一，目标不改。
'''
write('unified_physics_condition_ledger_670.md',ledger)
mapping={'joint_dynamic_scalar_gauge_interface':'joint_nonflat_mass_measure',
         '668':'669','669':'670','670':'671','3243':'3245','3245':'3247',
         '1316':'1319','1319':'1322','2354':'2364','2364':'2376'}
verify=replace((HERE/'verify_round669.py').read_text('utf8'),mapping)
verify=verify.replace('Verify actual scalar moment interface and joint mass/kinetic gauge requirement.',
    'Verify nonflat shared dictionary and regular physical observations at massless zero.')
verify=verify.replace('import joint_mass_auxiliary_reflection as prior_model','import joint_dynamic_scalar_gauge_interface as prior_model')
verify=verify.replace('dynamic_scalar_probe','nonflat_gauge_probe').replace('dynamic_scalar_entry','nonflat_gauge_entry')
verify=verify.replace("'round670_drafts/nonflat_gauge_probe_results.json')",
    "'round670_drafts/nonflat_gauge_probe_results.json',\n           'round670_drafts/first_attempt_mass_measure.py','round670_drafts/first_attempt_diagnostic.json')")
write('verify_round670.py',verify)
publish=replace((HERE/'publish_round669.py').read_text('utf8'),mapping|{'## 315.':'## 316.','## 220.':'## 221.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第670轮完成：** [非平坦规范场中的共同观测，以及无质量零模处的正则延伸]({p}research_note_670.md)'
         '原配对消去观测／质量字典中的逆无质量传播块；非交换非平坦链路、原质量与全部物理来源共同相容。'
         '无质量零模不再是该字典的排除条件；Wilson谱隙和动态正性边界保留。'
         '两组、十八式通过，最新670／3247，1322份编号科学文件、2376份保护证据。'
         '[核验]({p}research_round_670_checks.json)、[全条件账]({p}unified_physics_condition_ledger_670.md)。'
         '原完整过程、动态规范正性、连续及量子GR仍开放。')
order=('**当前执行顺序（670后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[671共同规范背景与物理反射泛函]({p}round671_drafts/STATUS.md)，'
       '核反射兼容的真实空间曲率背景，再接动态规范正时间组合；目标不改。')
'''+publish[b:]
publish=publish.replace('原动态标量、手征物质与共同规范过程的条件','非平坦规范场中的共同观测，以及无质量零模处的正则延伸')
publish=publish.replace('原动态标量与规范传播的联合条件','共同规范字典消去无质量零模限制')
publish=publish.replace('原热标量可积不替代共同规范动力学','正则观测字典与动态规范正性仍须区分')
write('publish_round670.py',publish)
write('postcheck_round670.py',replace((HERE/'postcheck_round669.py').read_text('utf8'),mapping))
write('round670_drafts/research_note_670_draft.md',(HERE/'research_note_670.md').read_text('utf8'))
review='''670 primary-agent algebra/code/report audit; no independent agent review.
Goal unchanged, no new application tasks, subagents, automations or image checks.
Existing Python/NumPy. Original 660,661,668,669 and670 entry reused explicitly.
Original reminders382-386,425,522-523 inherited; no revived Lipschitz condition,
no stacking alternative neighborhood bridges or redoing523 joint realization.
The first expectation of inverse-Kl coordinate divergence failed; both initial
code and diagnostic preserved. The original pairing obeys M Q+=Q+^T M.
Dv=Q+v yields exact factorization of C^T and F; inverse Kl cancels.
The new triangular map is determinant1 at Kl zeros without density assumptions.
Projector Phi is frame-free, norm bounded by golden ratio, still nonlocal in time.
This supersedes only the inverse-observation limitation in660's formula, not
the already proved massless weight extension, index changes or Wilson H gap.
Mass/source Grassmann polynomial identity survives even singular A/Nw; inverse
Wick/log formulas are used only where appropriate. Complete Pfaffian phase and
frame determinant preserved. Unit S9 bar pair Pf1 already established in660.
Actual noncommuting colour/weak/hypercharge links have genuine plaquette defect;
local gauge, scalar, auxiliary and physical sources are transported together.
Actual theta=pi/6 holonomy has H gap1 and exact massless q6 zero; original masses
give invertible physical Nambu kernel and positive witness weight, not universal
positivity. Full local covariance and Pfaffian cofactor tests included.
Numerical derivatives check a finite nonzero weight at the zero, not all models.
Sources keep Phi/L dependence. No fixed-background positive weight promoted to
full dynamic gauge RP, normalization, original Gibbs identity or continuum GR.
The external sources support only their stated original scopes, not our theorem.
Next671 physical reflection on a compatible curved spatial background; avoid
reproving661 nonlocal support or arbitrary nonsymmetric-background false no-go.
Two new groups,18 displayed numbered equations. Historical evidence immutable.
'''
for name in ('research_note_670.md','joint_nonflat_mass_measure.py','joint_nonflat_mass_measure_results.json','unified_physics_condition_ledger_670.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round670_drafts/final_review.txt',review)
print('670 report/review frozen; verification/publication prepared')
