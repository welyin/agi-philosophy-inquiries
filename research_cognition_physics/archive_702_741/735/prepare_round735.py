"""Prepare735 immutable proof scope, ledger and reproducibility scripts."""
import hashlib
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def write(name, content):
    target = HERE / name
    target.parent.mkdir(exist_ok=True)
    with target.open('x',encoding='utf8',newline='\n') as handle:
        handle.write(content)


def remap(text, mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda match:mapping[match.group()],text)


ledger = (HERE/'unified_physics_condition_ledger_734.md').read_text('utf8')
ledger = '# 联合条件总账：735一圈联合守恒处方的局部提升\n' + ledger.split('\n',1)[1]
ledger = ledger.replace('接[733全账](unified_physics_condition_ledger_733.md)，回填[734报告](research_note_734.md)。[结果](joint_retarded_reference_response_results.json)、[核验](research_round_734_checks.json)。',
    '接[734全账](unified_physics_condition_ledger_734.md)，回填[735报告](research_note_735.md)。[结果](joint_local_source_normalization_results.json)、[核验](research_round_735_checks.json)。')
marker = '## 当前共同对象及仍存在的分支'
assert marker in ledger
ledger = ledger.replace(marker,'**735当前增量：** 在明确的局部、质量无关UV规范化框架内，原完整二次费米场的一圈共同Ward处方可经有限局部jet匹配接到实际Hadamard来源；保原变化质量、非平坦连接、同一过去和记录。此为限定部门处方存在性，有限物理常数、绝对自洽解、全圈与共同连续映射仍开放。\n\n'+marker)
updates = {'C04':'735给所选局部一圈框架内的共同守恒来源；反馈导数界及非线性发展未签收',
           'C19':'735局部规范化不换过去准备；同背景实际记录来源差中局部项相消',
           'C22':'735原矩阵质量与完整联合Ward规范化经有限jet连接；只关闭声明的一圈外背景部门'}
rows = ledger.splitlines()
seen = []
for index,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value
            rows[index]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 735限定一圈部门的共同处方

- 采用标准局部、固定μ、质量无关且光滑至零质量的UV处方；它是明确物理规范化输入，不冒充认知定理。
- 原完整二次费米作用、复Y、Majorana及变化Higgs／singlet全部保留；只把背景质量视为矩阵spurion。623坐标不是新发现，也不使玻色理论变成可重整化理论。
- 一圈局部规范化需处理普通／混合反常、Abelian物质项及反场类；原对象接到共同无反常方案，未调用全圈EFT结论。
- 有限局部jet的系数匹配与唯一余项延拓给实际局部Ward提升；不假定完整非局部背景展开收敛，不把光滑函数的形式零级数当恒零。
- 同一Hadamard态可在该Wick处方取期望；原真实记录差、远过去记忆保持。原φ源及响应须含Jacobian拉回与其导数接触。
- 两组原完整矩阵检查核非交换局部对象及实际记录来源映射；不计算连续重整化应力或反常系数。
- 兼容有限常数、绝对基准解、统一来源导数界、装置来源、全Gauss及连续映射仍开放。F>0本身不提供全域统一界。

## 本轮合并与下一项

C04／C19／C22接成限定一圈、给定光滑背景的共同守恒来源合同。没有签收绝对半经典非线性解。

接[736](round736_drafts/STATUS.md)：回查已有来源正则性及高频结果，确定实际反馈所需函数空间、局部项阶数与非局部响应，避免再证明成熟Ward身份。旧空间、604、649／699及统一目标保持。
'''
write('unified_physics_condition_ledger_735.md',ledger)
write('round735_drafts/research_note_735_draft.md',(HERE/'research_note_735.md').read_text('utf8'))
write('round736_drafts/STATUS.md','''# 第736轮入口：守恒来源之后的实际反馈正则性

接[735](../research_note_735.md)、[条件账](../unified_physics_condition_ledger_735.md)。所选标准局部一圈框架内的共同来源处方已有条件性存在连接；绝对自洽背景及非线性发展仍开放。

1. 先回查632、634—635、663／667、730—735及历史高频／导数损失研究；不要重复一般光滑方向可微不等于Banach可微的口头说明。
2. 区分局部重整化项的最高导数、实际过去态的非局部响应、背景或质量坐标带来的非统一常数。F>0不等于有统一正下界，但不得把坐标边界直接判成物理奇点。
3. 对原完整来源给可检验的函数空间／有限阶展开合同，优先对接已有半经典发展定理，核原变化规范与标量背景的准确适用性。
4. 若二阶经典方程与未经降阶的一圈响应不在同一闭合估计内，区分完整方程、逐阶ħ解与有效理论降阶；不得把任选降阶当已证精确解。
5. 初始量子约束及共同准备仍须明确；不得在非解背景上直接调用围绕自洽解的稳定性定理。实际记录装置保634条件。
6. 不重算规范荷表、不扫有限矩阵精度、不新造认知原则来代替已有分析缺口。旧空间、649／699及统一目标保持。
''')
audit=dict(sources=[
    dict(url='https://arxiv.org/html/1501.07014',locations='Section8, equation8.1 and one-loop discussion following8.3',
         use='One-loop joint local anomaly normalization; specifically account for Abelian matter and antifield classes.',
         caution='Not importing the full all-orders EFT theorem or its cohomological-completeness enlargement as a physical action.'),
    dict(url='https://arxiv.org/html/1210.4031',locations='Section3.2 scaling expansion and Proposition3.14',
         use='Finite local jet expansion and extension of the singular coefficients, used with the explicitly mapped matrix-mass background.',
         caution='Section4 stress theorem is not applied to variable charged chiral masses. Matrix-mass extension and Ward normalization are separate steps in note735.')],
    independent_assumptions=['Four-dimensional ordinary spin background and original representation',
       'Only quadratic fermions quantized; smooth given backgrounds in F>0 patches',
       'Mass-independent, fixed-scale local UV normalization with standard scaling and zero-mass smoothness'],
    inference='Combine one-loop formal joint normalization with finite UV-jet matching for the original matrix mass, then evaluate the common Wick sources in the original Hadamard state.',
    excluded=['Convergence of a nonlocal effective-action background expansion','All-loop quantum gravity',
              'Numerical continuum stress coefficients','Unique finite physical constants','Nonlinear self-consistent spacetime existence'])
write('round735_drafts/literature_scope_audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
write('round735_drafts/scope_and_dedup_review.md','''# 735范围、证明与去重审查

主代理审查，无独立代理，无图像检验。目标未修改。

- 623质量坐标、628完整反常类别、730实际连续参考、732—734相对来源及退迟接触直接继承。
- 735入口和工作报告的边界保留；不覆盖其冻结文件。新增连接是一圈规范化与有限局部jet提升，不是声称前几份范围受限的文献本来就覆盖原问题。
- UV权赋给外部质量参数，是辅助局部缩放；没有把原曲目标玻色作用变成可重整化理论。
- 矩阵质量和规范连接不对易，局部递推保有序积。手征／自对偶字典及Pfaffian相位保持。
- Abelian物质反常另用全局U(1)保持的调节论证；未把有限矩阵荷守恒当连续反常计算。
- 形式提升只作用于有限局部奇异系数，余项选到Ward求导后仍唯一延拓；不推断非局部有效作用级数收敛。
- 共同有限匹配规定一份Wick来源合同；不宣称已有全部系数的显式数表。任意未指定parametrix不因此自动守恒。
- 实际态来源通过同一个Wick定义求值；更换局部规范化不等于更换准备。记录差中的局部项相消，实际响应仍有坐标及局部接触。
- 两项数值是对象映射和非零接触核验，不是连续重整化应力。额外三个历史工作probe不重复加入正式累计数。
- 兼容有限物理常数与F域上的统一估计仍开放；局部存在处方不提供半经典非线性解。
- 旧空间382—386、425、522—523不重证，不重新加入已消去的假设。原图连续映射649／699保持。
''')
main=('research_note_735.md','joint_local_source_normalization.py','joint_local_source_normalization_results.json','unified_physics_condition_ledger_735.md')
write('round735_drafts/final_review.txt','Primary-agent review only. Conditional one-loop local normalization bridge; no nonlinear semiclassical existence claimed.\n'+
      '\n'.join(name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in main)+'\n')

verify=(HERE/'verify_round734.py').read_text('utf8')
verify=remap(verify,{'735':'736','734':'735','733':'734','3270':'3284','3284':'3304','3407':'3409','1514':'1517',
                    'joint_retarded_reference_response':'joint_local_source_normalization',
                    'joint_reference_polarization_boundary':'joint_retarded_reference_response',
                    'causal_reference_entry':'counterterm_scope_entry'})
verify=verify.replace("(3,0,0)","(2,0,0)").replace('run=3,','run=2,').replace("==20","==16")
verify=verify.replace("'round735_drafts/entry_checks.json')\n    new=", "'round735_drafts/entry_checks.json',\n"+
    "           'round735_drafts/gauge_history_scope_probe.py','round735_drafts/gauge_history_scope_probe_results.json',\n"+
    "           'round735_drafts/ward_literature_scope_review.md','round735_drafts/research_note_735_working.md',\n"+
    "           'round735_drafts/check_and_publish_followup.py','round735_drafts/followup_checks.json')\n    new=")
anchor="    main=('research_note_735.md'"
pos=verify.index(anchor)
verify=verify[:pos]+"    for name,digest in core.read(HERE/'round735_drafts/followup_checks.json')['artifact_hashes'].items():\n        assert core.digest(HERE/'round735_drafts'/name)==digest,name\n"+verify[pos:]
write('verify_round735.py',verify)

publish=(HERE/'publish_round734.py').read_text('utf8')
publish=remap(publish,{'735':'736','734':'735','733':'734','3407':'3409','3404':'3407','1514':'1517','3284':'3304','380':'381','285':'286',
                      'joint_retarded_reference_response':'joint_local_source_normalization',
                      '同一过去准备的退迟来源与完整接触条件':'一圈联合守恒处方的局部提升与同一实际来源'})
start=publish.index("summary='");end=publish.index('planned={}',start)
publish=publish[:start]+'''summary='**第735轮完成（限定一圈框架）：** [联合守恒处方与实际来源]({p}research_note_735.md)在声明的局部UV规范化框架内，原完整变化质量经一圈规范化和有限jet匹配接到共同守恒来源；保实际过去与记录，响应仍需完整接触。两组、十六式通过，最新735／3409，1517份编号科学文件、3304份保护证据。[核验]({p}research_round_735_checks.json)、[条件账]({p}unified_physics_condition_ledger_735.md)。有限物理常数及非线性自洽解仍开放。'
order='**当前执行顺序（735后，优先于下方历史安排）：** 接[736实际反馈正则性]({p}round736_drafts/STATUS.md)，先复用旧高频／来源结果，区分局部项、非局部态响应和降阶合同；不再重复荷表或有限矩阵优化。旧空间、604、649／699及统一目标保持。'
'''+publish[end:]
# The unique publication marker includes its explicit scope qualification.
publish=publish.replace("'**第735轮完成：**'", "'**第735轮完成（限定一圈框架）：**'")
write('publish_round735.py',publish)
post=(HERE/'postcheck_round734.py').read_text('utf8')
post=remap(post,{'735':'736','734':'735','3284':'3304'})
post=post.replace("'**第735轮完成：**'", "'**第735轮完成（限定一圈框架）：**'")
write('postcheck_round735.py',post)
print('Prepared735 ledger, preserved scope and publication scripts.')
