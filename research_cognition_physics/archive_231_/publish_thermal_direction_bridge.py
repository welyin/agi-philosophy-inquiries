"""Publish the completed conditional common-model report, preserving history."""
from pathlib import Path
import hashlib
import json
import os
import verify_thermal_direction_bridge as evidence
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
before = evidence.verify()
assert before['total_protected_including_this_review']==1052
paths = [ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
raws = {p:p.read_bytes() for p in paths}
folder = HERE/'navigation_before_thermal_direction_bridge_20260930'
folder.mkdir(exist_ok=False)
summary = (
    '**热参考与方向上界已有共同实现（未编号）：** '
    '[正式报告]({p}thermal_direction_bridge_review.md)在明确新增的连续有效几何和热来源下，'
    '把实际坐标差读、经典记录与合法qubit仪器接成同一过程。'
    '固定绝对精度给d>2，复用382反向效果上界给d≤3，得到这组条件下的d=3；'
    '482种方向设置仍有统一反向对比≥0.407523。二维方向可辨及四维大载体两个反例明确了条件边界。'
    '连续空间、完整qubit接口及绝对精度要求仍为输入，未证明FUCP唯一生成三维。'
    '6组诊断、13式及独立复核通过；[核验]({p}thermal_direction_bridge_checks.json)。'
    '编号仍520／2552，871份编号科学文件；保护证据1052份。')
next_steps = '''

### 当前接续：条件性三维共同模型已成报告，检验实际测量的局域性

新报告已将热精度下界与382上界接到同一有效场、同一Gaussian坐标差记录和qubit反馈仪器。未知测试qubit及其旧参考与场源独立，完整仪器保持完全正；平均效果与真实后态不能混用。有限482方向菜单的统一反向对比来自解析覆盖和方差界，数值只作复核。

**已知条件边界：** 二维在相同热模型中可有规模一致方向对比，绝对坐标精度却失败；四维在更大的完整载体上同时允许热精度与反向读取。因此新增绝对精度任务及完整qubit方向接口均有实质内容。它们可以作为候选原则继续研究，不必先证明为FUCP推论；但其来源不能被省略。连续空间、热源和固定UV尺度同样保留为逆向模型输入。

**下一项具体连接：** 当前球形UV截断使所谓点读口带有空间尾部，不能直接继承严格局域操作或因果组合。考虑把UV控制转移到有限分辨率的平滑探针，而保留未截断的局域自由场：复用已审Fewster–Verch测量框架，核验实际有限支撑测量是否保留热差读的红外阈值、有限回冲和qubit反向对比，并说明记录只能在共同因果未来汇集。需要给同一量子过程或明示的受控近似，不能仅在积分中乘一个截止函数就宣布局域性已经解决。

此检验服务空间读数与因果结构的相容性；不预先宣称由它生成了洛伦兹对称、原动态图连续极限或Einstein方程。若必须增加作用、状态或探针条件，按双向统一目标列明并继续检验，而不将新增条件本身视为路线失败。

后续仍需连接参考来源与维护、原始关系事件到有效几何、共同时间及几何反馈。热化不是给定初始热边缘唯一可能的解释，也不要求重复扫描同一噪声积分。候选框架的条件性结果与最终空间生成结论分开记录，当前空间阶段和完整统一未结项。

用户目标保持不变；研究结果完成后及时保存正式报告和索引。最新冻结入口verify_thermal_direction_bridge.py，保护1052份证据；编号520／2552保持。手动继续，无定时任务。
'''
manifest, planned = {}, {}
for p,raw in raws.items():
    key = ('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
    with (folder/key).open('xb') as stream: stream.write(raw)
    manifest[str(p.relative_to(ROOT))] = dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(enc).replace('\r\n','\n')
    assert '热参考与方向上界已有共同实现（未编号）' not in body
    prefix = 'research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    block = summary.format(p=prefix)
    if p in paths[:5]:
        title,tail = body.split('\n\n',1)
        body = title+'\n\n'+block+'\n\n'+tail
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 167.' not in body
        body += '\n\n## 167. 热坐标与qubit反向效果的同一实现\n\n'+block+next_steps
    else:
        assert '\n## 72.' not in body
        body += '\n\n## 72. 条件性三维的共同模型与两条不可偷换的要求\n\n'+block
        body += '\n\n热下界、实际方向效果和完整仪器已在新增连续有效模型中共同成立。原格点到连续壳并未因此证明；两条额外任务要求和来源见正式报告。下一项核验有限支撑探针与该模型的因果相容性。\n'
    if p.name=='RESEARCH_STATE.md': body += next_steps
    if p.name=='research_direction.md':
        body += '\n\n**最新接续（目标保持）：** 条件性三维共同模型已正式保存。下一项检验有限分辨率、有限支撑的真实局域探针能否替代球形UV截断，并同时保留热精度、方向对比和因果组合；不重复一般Gaussian或BU证明。详见RESEARCH_STATE.md末尾。\n'
    planned[p] = body.replace('\n',nl).encode(enc)
with (folder/'manifest.json').open('x',encoding='utf8') as stream:
    json.dump(manifest,stream,ensure_ascii=False,indent=2)
count = 0
for p,body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        count += 1
for p,raw in raws.items(): assert p.read_bytes()==raw,('concurrent change',p)
for p,body in planned.items():
    temp = p.with_name(p.name+'.thermal-direction.tmp')
    with temp.open('xb') as stream: stream.write(body)
    os.replace(temp,p)
assert evidence.verify()==before
report = dict(date='2026-09-30',scientific_base_through_round=520,
    unchanged_numbered_scientific_tests=2552,numbered_scientific_files=871,
    protected_evidence=1052,navigation_files=7,navigation_links=count,broken_links=0,
    previous_navigation_snapshots_preserved=True,independent_cross_review_completed=True,
    scientific_evidence_unchanged_after_navigation=True,active_goal_unchanged=True,
    manual_goal_start_preserved=True,stage_complete=False,all_checks_passed=True)
with (HERE/'thermal_direction_bridge_integration_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
