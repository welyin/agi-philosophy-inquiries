"""Register the completed joint-model round without overwriting earlier evidence."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauge_matter_constraints as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round531_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_531_checks.json')
assert checked['round']==531 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round531_20260930'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths}
    folder.mkdir(exist_ok=False)
    manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream:
            stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)
summary=('**第531轮完成（联合模型主线）：** [同一物质表示、反常与规范动能]({p}research_note_531.md)'
    '给出正模块迹实现三种规范耦合的精确充要条件；同一表示同时限制超荷、中心核与动能归一化。'
    '反常保留多色数竞争解，额外输入同尺度颜色／弱耦合相等才在该族内选三色；这不是空间三维。'
    '历史三代一圈基准可满足加权条件，未做当前实验拟合或完整谱模型实现。'
    '9组检查与独立终审通过，最新531／2631，905份编号科学文件，'
    f"{checked['cumulative_unique_protected_evidence_files']}份保护证据；"
    '[核验]({p}research_round_531_checks.json)、[联合条件更新]({p}unified_physics_condition_ledger_531.md)。'
    '本段及下述接续优先于保留的旧轮次下一步；锚点草稿未计入本轮，统一目标仍未完成。')
next_steps='''

### 第531轮后：审计共同权重是否能覆盖整个模型

531将共同候选、物质表示、反常、正迹动能、全局表示核及同尺度匹配放入同一条件账。新增精确接口为x_Y=3x_w−(N²−1)x_c/(2N)、x_w>Nx_c/4（x_c>0），这是在给定表示及共同正迹规则下的充要条件。颜色数、空间维数及概率迹分别审计；没有把成熟一般色数族或无权耦合失配计为新发现。

下一项检查同一正权w_q,w_l与有限表示、Yukawa／Dirac、实结构及组合的兼容性，并将其影响同时追踪到Higgs和几何部门。规范作用的迹、概率归一化的迹和认知组合不是仅凭名称就能等同。若完整结构要求等权，则必须处理已有匹配缺口；若允许不等权，则记录自由度来源、其它部门限制及可检验后果。直接对接成熟谱模型，不另造相同运行曲线，也不回到无限优化锚点来源。

维数、共同g、群族、表示类型、作用和Majorana机制仍有输入；真实主体、参考、微观到场论、熵、宇宙学和量子引力未闭合。条件账首版冻结保留，531回填版为当前接续。正式新研究仍及时交付编号报告、复算结果及范围；下一编号532。应用目标已由用户调整为联合物理目标；不再修改目标，不创建任务或定时任务。
'''
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第531轮完成（联合模型主线）：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 177.' not in body
        body+='\n\n## 177. 联合规范—物质条件与空间维数的边界\n\n'+block+next_steps
    else:
        assert '\n## 82.' not in body
        body+='\n\n## 82. 颜色数的条件选择不是空间三维\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为530／2622','最新科学轮次与检查数为531／2631')
        body=body.replace('完成231—530轮。','完成231—531轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第531轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|531|[规范—物质联合条件](research_note_531.md)|[代码](joint_gauge_matter_constraints.py)、'
            '[结果](joint_gauge_matter_constraints_results.json)、[核验](research_round_531_checks.json)|\n')
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body:
        continue
    temp=path.with_name(path.name+'.round531.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream:
            stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=531,cumulative_tests=2631,numbered_scientific_files=905,
    unique_protected_evidence_files=checked['cumulative_unique_protected_evidence_files'],
    navigation_files=7,navigation_links=links,broken_links=0,navigation_snapshots_preserved=True,
    recovered_partial_navigation_publish=resuming,complete_report_published=True,
    existing_results_preserved=True,next_round=532,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
