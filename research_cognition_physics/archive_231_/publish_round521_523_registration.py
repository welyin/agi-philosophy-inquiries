"""Restore consecutive research-note navigation for three completed studies."""
from pathlib import Path
import hashlib
import json
import os
import verify_round521_523_integration as evidence
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
before=evidence.verify()
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',
       HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
raws={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round521_523_registration_20260930'
folder.mkdir(exist_ok=False)
summary=(
    '**第521—523轮已完整归档：** '
    '[521：正势下空白探针接入]({p}research_note_521.md)、'
    '[522：热参考精度与维数边界]({p}research_note_522.md)、'
    '[523：热参考与qubit方向的共同实现]({p}research_note_523.md)。'
    '三篇均含完整假设、推导／反例、代码、结果和范围，不是链接占位。'
    '此前误作未编号接口的18组检查现一次性纳入编号统计，未另做或重复登记实验。'
    '最新523／2570，880份编号科学文件，1058份唯一保护证据；'
    '[整合核验]({p}round521_523_integration_checks.json)。'
    '后继完整研究继续使用research_note_524.md起的编号；下文旧未编号表述和520计数保留为登记前历史，以本段为准。')
next_steps='''

### 第523轮后：保持连续编号，及时保存完整报告

按用户2026-09-30明确要求，521—523已从现存完整接口报告登记为正式编号研究。原报告、代码、JSON和历史冻结检查全部保留；18组检查从未编号转入编号统计，不计作额外新实验。521和522独立接续520之后的既有接口；523使用522热参考结果。完整正文、来源和范围均未缩写为摘要。

以后每个完成的完整研究单元及时交付research_note_编号.md、相应代码／结果及检查记录，并更新索引。简单勘误、纯文献核对或同一结果的重述仍作为原轮补充，不能为了编号重复实验。目标内容保持不变。

**下一轮为524：** 沿523已记录的具体连接，检验有限分辨率、有限支撑的实际探针能否替代球形UV截断，在同一过程保留热差读阈值、有限回冲、qubit方向对比和因果组合。优先复用已审成熟测量定理，明确新增模型条件；不把平滑一个积分核直接当成局域仪器的实现。条件性三维已有共同模型，但原动态图到有效几何、来源与控制、统一几何反馈仍开放，空间阶段尚未结项。

最新核验入口verify_round521_523_integration.py；最新编号523／2570，保护1058份证据。手动继续，不设置定时任务。
'''
manifest,planned={},{}
for p,raw in raws.items():
    key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
    with (folder/key).open('xb') as stream:stream.write(raw)
    manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    assert '第521—523轮已完整归档' not in body
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if p in paths[:5]:
        title,tail=body.split('\n\n',1)
        body=title+'\n\n'+block+'\n\n'+tail
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 168.' not in body
        body+='\n\n## 168. 第521—523轮连续编号与正式报告登记\n\n'+block+next_steps
    else:
        assert '\n## 73.' not in body
        body+='\n\n## 73. 521—523完整研究报告登记\n\n'+block
    if p==HERE/'README.md':
        body+='\n\n## 第521—523轮正式研究索引\n\n|轮次|完整报告|代码、结果与核验|\n|---|---|---|\n'
        rows=[(521,'正势下空白探针接入','positive_reference_acquisition','positive_reference_acquisition.py'),
              (522,'热参考精度与维数边界','thermal_reference_dimension','thermal_reference_dimension_model.py'),
              (523,'热参考与方向仪器的共同实现','thermal_direction_bridge','thermal_direction_bridge_model.py')]
        for n,title,stem,code in rows:
            body+=f'|{n}|[{title}](research_note_{n}.md)|[代码]({code})、[结果]({stem}_results.json)、[核验](research_round_{n}_checks.json)|\n'
    if p.name=='RESEARCH_STATE.md':body+=next_steps
    if p.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为520／2552','最新科学轮次与检查数为523／2570')
        body=body.replace('完成231—520轮。','完成231—523轮。')
        body+='\n\n**报告格式与下一轮：** 用户要求每个完整研究及时写入research_note_编号.md；521—523已补齐，原未编号材料保留，后继从524继续。研究目标不变，具体下一问题见RESEARCH_STATE.md末尾。\n'
    planned[p]=body.replace('\n',nl).encode(enc)
with (folder/'manifest.json').open('x',encoding='utf8') as stream:
    json.dump(manifest,stream,ensure_ascii=False,indent=2)
links=0
for p,body in planned.items():
    for link in core.link_parser()(body.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
for p,raw in raws.items():assert p.read_bytes()==raw,('concurrent change',p)
for p,body in planned.items():
    temp=p.with_name(p.name+'.round521-523.tmp')
    with temp.open('xb') as stream:stream.write(body)
    os.replace(temp,p)
assert evidence.verify()==before
report=dict(date='2026-09-30',rounds=[521,522,523],stage_saved_tests=2570,
    numbered_scientific_files=880,unique_protected_evidence_files=1058,
    existing_diagnostics_registered_once=18,new_experiments_for_registration=0,
    complete_numbered_notes=3,navigation_files=7,navigation_links=links,broken_links=0,
    old_reports_code_results_unchanged=True,navigation_snapshots_preserved=True,
    next_round=524,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with (HERE/'round521_523_navigation_checks.json').open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
