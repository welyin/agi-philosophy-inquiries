"""Prepare732 relative source and short-time response evidence."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_731.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：732同一记录的相对来源与短时联合发展\n'+rest
ledger=ledger.replace(
 '接[730全账](unified_physics_condition_ledger_730.md)，回填[731报告](research_note_731.md)。[结果](joint_source_constraint_response_results.json)、[核验](research_round_731_checks.json)。',
 '接[731全账](unified_physics_condition_ledger_731.md)，回填[732报告](research_note_732.md)。[结果](joint_relative_source_development_results.json)、[核验](research_round_732_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**732当前增量：** 原730非选择记录的光滑态差由至多四个完整Dirac模式表达，给共同相对应力、电流与质量力；光滑双解核满足联合Ward。接731初始补偿和573固定规范波系统，得到保全部线性约束的短时相对响应。绝对基准、完整非线性反馈及实际仪器仍未完成；矩阵校准不冒充连续Einstein数值解。\n\n'+marker)
updates={
 'C04':'732实际同态相对源、原联合初始补偿及保约束的短时线性发展接通，非完整非线性半经典发展',
 'C19':'732光滑记录源差的共同反项相消不固定绝对宇宙参考；变化背景初态需另核',
 'C22':'732同一有限模式给相对应力／电流／五标量力和联合交换身份，不能仅修正能量'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 732同一记录的相对来源与联合时间发展

- 730实际非选择CAR记录给秩至多4的光滑二点差，权重可有正负；这些是来源差的模式，不是新粒子。
- 原二次作用的极化Noether身份作用于光滑双解差核，给相对应力、电流、五标量力的共同Ward。两态同一局部反项相消，不证明绝对理论已闭合。
- 原参考颜色不变，731初始右逆适用；其几何及物质速度字典一并变分，不能冻结旧K或旧速度。
- 573主部及601同阶方法直接复用。固定规范后光滑相对源驱动唯一短时线性发展，零初始约束误差齐次传播。
- 731所用正超解身份已见573，不重新作为新定理；本轮新增是实际同态来源到联合时间发展的具体连接。
- 原128维校准核四模式／完整协方差、全部12生成元和四类背景做功；不是全时空Einstein求解或连续数值证书。
- 差解不能反推两份绝对半经典解存在；634无补偿跳接限制、580／600有效项和原连续边界保持。

## 本轮合并与下一项

C04／C19／C22的相对、领先响应在同一实际记录上进一步合并。统一目标仍开放。

接[733](round733_drafts/STATUS.md)：检验有限相对模式与改变背景上的实际Hadamard参考／同一记录初始化能否共同实现，再决定有限反馈的可闭合范围。停止一般Ward、重复线性波证明和小矩阵精度优化。旧空间、604、649／699及统一目标保持。
'''
write('unified_physics_condition_ledger_732.md',ledger)
write('round732_drafts/research_note_732_draft.md',(HERE/'research_note_732.md').read_text('utf8'))
write('round733_drafts/STATUS.md','''# 第733轮入口：变化背景的实际参考与相对模式初始化

接[732](../research_note_732.md)、[全账](../unified_physics_condition_ledger_732.md)。当前已有同一实际记录的光滑相对来源、联合初值及短时线性响应，仍没有完整非线性量子—几何闭合。

1. 回查573、601、633—635、649—651、704、730—732。不得重报一般Ward、二点差光滑或同阶约束传播。
2. 先核改变几何后，同一光滑记录及其有限相对模式能否来自合法CAR／Hadamard初态。不能随便指定几份有符号波函数就称它们是实际量子态差。
3. 可检验接法：保实际记录所在的有限光滑子空间，把其协方差嵌入新背景的Hadamard参考；核正性、自对偶性、原记录后态及全部来源字典。若改变参考选择，明确记为选择而非免费制备。
4. 若初始化可行，再检验共同有限模式—物质—几何过程的PDE与约束；不得把经典带符号辅助波模型直接称为完整量子理论。
5. 原绝对背景来源、有限反项和物理制备仍须交代；相对源抵消不等于把它们设零。
6. 旧空间、原图连续、604、649／699及统一目标保持。存在真实共同增量才编号；不以调参或精度检查增加轮次。
''')
write('round732_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[dict(url='https://arxiv.org/html/1210.4031',
   use='Section4 context for background Dirac sources and state-independent renormalization.',
   caution='Its separate conserved-stress result has constant-mass restrictions. Here joint Ward is derived for smooth same-background bisolution differences from the original quadratic invariant action.'),
 dict(url='https://arxiv.org/abs/1909.01876',use='Inherited573 wave/Lorenz reduction method.',
   caution='Do not import its SU2 flat-target noncompact small-data global theorem into this compact curved-target model.'),
 dict(url='https://arxiv.org/abs/gr-qc/9211002',use='Inherited601 perturbative/order-reduction context.',
   caution='Not a full nonlinear existence theorem for the project model.')],
 new='Original actual record finite-mode relative source, complete joint Ward, existing initial right inverse and fixed-gauge short-time linear solution connected.',
 excluded='Absolute renormalized solution, nonlinear state-dependent source map, physical detector and GR generation.'
),ensure_ascii=False,indent=2)+'\n')
write('round732_drafts/scope_and_dedup_review.md','''# 732范围、推导及代码审查

主代理审查；无独立代理、无图像检验。

- 上轮731及入口均为有效进展；本轮开始没有运行中的Python。历史结果不覆盖。
- 原记录R的秩2作用给DeltaP秩至多4；后态无需Gaussian。有限模式用于二次来源，不增加物种或把有符号差当正态。
- 同背景光滑双解差可直接取局部导数；共同c-number反项／异常相消。此结论不推出绝对量子理论无反常。
- 源定义按同一原二次作用的变分固定。变化质量和连接要求联合Ward，不能套用无外场的单独应力守恒。
- 573已有约化主部与规范传播、601已有共同一阶方法、731已有联合初值右逆；新连接是当前实际记录源满足其全部接口。
- 731正超解身份在573已经出现；本轮和后续不再宣称它是新发现。
- 线性PDE的光滑源由原真实背景上的模式给出。代码未来一阶jet仅核有限局部源，不宣称已数值解573未来Einstein系统。
- 检验保全部原质量、规范作用、五标量和非零几何速度。圈／弱荷交换及几何／物质做功有真实非零项。
- 第一次运行因旧Pauli容器是list而失败，转np.asarray修复；未改变物理参数或断言范围。
- 一阶相对解不证明两份绝对解存在，也不允许634的无补偿时空跳接。高阶项／绝对源／参考准备仍保留。
- 下一项检验变化背景上有限相对数据能否嵌入真实Hadamard记录态，而不是继续小矩阵验证。
''')
files=('research_note_732.md','joint_relative_source_development.py','joint_relative_source_development_results.json',
       'unified_physics_condition_ledger_732.md')
write('round732_drafts/final_review.txt','732 primary review; no independent agent review.\n'+
 '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
 'Three checks pass. Same-record smooth relative source, complete exchange identities and conditional linear Cauchy development. No full semiclassical existence or absolute source claim. Goal remains open.\n')
v=(HERE/'verify_round731.py').read_text('utf8')
v=remap(v,{'joint_source_constraint_response':'joint_relative_source_development',
 'joint_dynamic_continuum_reference':'joint_source_constraint_response',
 'round732':'round733','round731':'round732','round730':'round731',
 '_731':'_732','_730':'_731','range(584,731)':'range(584,732)',
 '3228':'3242','3242':'3256','==20':'==18',
 'round=731':'round=732','round=730':'round=731',
 '3398':'3401','1505':'1508','gauge_counterflow_entry':'source_feedback_entry'})
write('verify_round732.py',v.replace('Reproduce731','Reproduce732'))
post=(HERE/'postcheck_round731.py').read_text('utf8')
post=remap(post,{'range(584,732)':'range(584,733)','round731':'round732','_731':'_732',
 '3242':'3256','第731':'第732','round=731':'round=732'}).replace('Check731','Check732')
write('postcheck_round732.py',post)
pub=(HERE/'publish_round731.py').read_text('utf8')
pub=remap(pub,{'732':'733','731':'732','730':'731','3398':'3401','3395':'3398',
 '1505':'1508','3242':'3256','## 377.':'## 378.','## 282.':'## 283.',
 'joint_source_constraint_response':'joint_relative_source_development',
 '原物质上的联合初始约束补偿与量子来源接口':'同一记录的相对量子来源与短时联合发展',
 '旧空间合同保持，初始约束补偿不替代完整自洽发展':'旧空间合同保持，相对线性发展不替代完整非线性反馈'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第732轮完成：** [同一记录的相对量子来源与短时联合发展]({p}research_note_732.md)原实际记录的有限光滑模式给共同应力、电流与质量力；联合Ward接原初始补偿及短时线性发展，约束沿时间保持。三组、十八式通过，最新732／3401，1508份编号科学文件、3256份保护证据。[核验]({p}research_round_732_checks.json)、[全条件账]({p}unified_physics_condition_ledger_732.md)。绝对源与完整非线性反馈仍开放。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（732后，优先于下方历史安排）：** 接[733变化背景的实际参考与相对模式]({p}round733_drafts/STATUS.md)，核有限相对资料能否共同嵌入Hadamard态和同一记录，再研究反馈。停止重复Ward、线性波及小矩阵精度优化；旧空间、604、649／699及统一目标保持。'"
write('publish_round732.py','\n'.join(lines)+'\n')
