"""Publish checked interface and the user's bidirectional stage objective."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checks=json.loads((HERE/'remote_witness_interface_checks.json').read_text('utf8'))
assert checks['all_reported_checks_passed']
assert checks['total_protected_including_this_review']==1017
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
snapshot=HERE/'navigation_before_bidirectional_stage_20260930'
snapshot.mkdir(exist_ok=False)
stage=('**当前阶段：双向约束与逻辑统一（用户2026-09-30调整）。** '
 '允许从成熟物理理论反推候选认知原则，与认知原则的正向推导共同约束；'
 '近期追求理论自洽相融，不再把所有条件必须先由原原则单向推出作为唯一推进门槛。'
 '统一须有共同数学对象、明确对应和相容的适用极限；逆向输入、已证关系、近似和预测机会分别标明。'
 '三维仍是近期子问题，可与动力学、记录结构联合求解。经验正确性面向可区分的预测检验，'
 '同时保持与已有观测相容；当前尚未完成逻辑统一或三维／GR推导。')
detail="""

### 当前阶段的验收与接续（2026-09-30）

用户明确允许从现有物理反求认知原则，并把当前阶段设为理论自洽相融、逻辑统一。本段覆盖旧导航中将“先证明全部新增条件来自原认知原则”作为唯一研究顺序的要求；不改写历史定理的假设和结论。严格的认知原则必然推出GR或明确不足性，保留为更强研究目标，不把它和当前统一框架的阶段完成混同。

- **双向提出条件。** 可从量子操作、因果集、相对论、引力或场论的成熟结构反推候选认知条件；按来源标记，不因是逆向输入而停止探索。与本项目已完成内容去重。
- **统一必须具体。** 列共同对象、状态、过程、组合、事件与记录；给不同理论之间的数学映射、共同接口、适用域和极限。将互不相联的定义并列，尚不算统一。
- **相容性与必然性分开。** 可以采用有效Einstein作用量、规范结构等作为候选输入来检验相容性；若这样做，不能把结果写成原FUCP独自推出该输入。只在相应假设下声明证明。
- **暂不要求所有自由参数已唯一决定。** 记录未定参数、额外结构、替代实现和区分它们的预测；经验正确性面向未来检验，并与现有观测保持相容。
- **近期主线。** 三维空间仍是子目标，但允许其与共同事件、传播和动力学联合确定；不预先要求标准模型全部完成。先利用391共同记录、287实际重叠及503回执的既有工具，避免重做一般复制或反求共同代数。

本次远端见证只是条件性接口连接，不是双向统一完成。下一项要判断：远端如何实际获得正确的来源／轮次对应和核验时机，或者从已有物理的关系观测量结构反推一种可供模型实现的合同。两类入口都允许；输入账和同一过程内的相容性仍需核验。普通标签转存、更多返回率扫描和未改变结论的矩阵放大不再增轮次。

既有自动研究已同步这一阶段目标，保留每天07:00运行及只在有意义结果、完成、失败或需要决定时通知的设置。后继证据核验入口为verify_remote_witness_interface_review.py，继承1017份保护证据；最新编号状态保持520／2552，未宣告三维、GR或逻辑统一阶段结项。
"""
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    assert '当前阶段：双向约束与逻辑统一（用户2026-09-30调整）' not in body
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
            'archive_231_/' if p.parent==RESEARCH else '')
    first,rest=body.split('\n',1)
    body=first+'\n\n'+stage+'\n'+rest
    entry=(f'**远端见证接口推论（未编号）：** [正文]({prefix}remote_witness_interface_review.md)、'
        f'[代码]({prefix}correlated_remote_witness.py)、[结果]({prefix}correlated_remote_witness_results.json)、'
        f'[核验]({prefix}remote_witness_interface_checks.json)。'
        '以新增本地颜色—记忆双翻作用把503回执与远端驻留位联系，'
        '根收到回执即可认证远端位；保留全部根读数时，本地复位只更新已知异或部门，无需远端清位。'
        '源率、重返与保存继承旧工具，故不另增521；原510 H及任意远端监测不在等价声明内。'
        '4项诊断、7式和独立终审通过；编号仍520／2552，871份编号科学文件，保护证据1017。')
    if p==HERE/'spatial_premise_closure_audit.md':
        body+='\n\n## 160. 双向研究标准与根回执—远端见证接口（未编号）\n\n'
    elif p==HERE/'three_dimensional_four_conditions_review.md':
        body+='\n\n## 65. 双向统一允许引入候选条件，空间结论仍按范围验收\n\n'
    else:
        body+='\n\n'
    body+=entry+'\n'
    if p in (RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',HERE/'spatial_premise_closure_audit.md'):
        body+=detail
    if p==HERE/'three_dimensional_four_conditions_review.md':
        body+='\n允许逆向提出空间合同，不等于这些合同已经成立；共同记录、角色协变及新驻留位均未自动给出三维。后续须在同一候选模型中检验几何、动力学与组合要求的相容性，并标记哪些维数选择仍为输入。\n'
    planned[p]=body.replace('\n',nl).encode(enc)
    label=('project_README.md' if p==ROOT/'README.md' else
           'research_'+p.name if p.parent==RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:f.write(raw)
for p,raw in baseline.items():assert p.read_bytes()==raw,str(p)
for p,data in planned.items():
    assert p.read_bytes()==baseline[p],str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{
        'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},
        f,ensure_ascii=False,indent=2)
print(json.dumps(dict(navigation_files_updated=len(paths),old_versions_saved=len(paths))))
