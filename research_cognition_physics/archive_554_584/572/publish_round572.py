"""Publish round 572; preserve all navigation snapshots and frozen science."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_gauss_einstein_initial_data as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round572_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_572_checks.json')
assert checked['round']==572 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round572_20261001'
resuming=folder.exists()
if resuming:
    manifest=core.read(folder/'manifest.json')
    raws={p:(folder/manifest[str(p.relative_to(ROOT))]['snapshot']).read_bytes() for p in paths}
    for p,raw in raws.items():
        assert hashlib.sha256(raw).hexdigest()==manifest[str(p.relative_to(ROOT))]['sha256']
else:
    raws={p:p.read_bytes() for p in paths}; folder.mkdir(exist_ok=False); manifest={}
    for p,raw in raws.items():
        key=('root_' if p.parent==ROOT else 'research_' if p.parent==RESEARCH else 'archive_')+p.name
        with (folder/key).open('xb') as stream: stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)

summary='**第572轮完成：** [共同规范物质与全局引力初值]({p}research_note_572.md)证明原571来源不能原样接入平坦共形周期CMC分支；同一Higgs与singlet的有限反流保Gauss后，可同时满足全部经典Einstein初始约束，并有固定资料下的唯一正共形因子。八组、十六式及独立终审通过，最新572／2915，1028份编号科学文件，1487份保护证据；[核验]({p}research_round_572_checks.json)、[条件账]({p}unified_physics_condition_ledger_572.md)。来源已改变，维数、二导数作用和有效匹配仍输入；低曲率、量子连续与统一仍开放。'
next_steps='\n\n### 第572轮后：先合并共同实现，再攻参考与演化接口\n\n共同Higgs、singlet和规范源现在能在明确反流后，同时满足Gauss及全部周期Einstein初始约束。原来源在平坦共形CMC分支的总动量失败证书保留；不能把有限改变的状态称为原样接入。固定资料下的共形因子唯一，不代表全模型或参数唯一。\n\n继续落实“先合并可复用条件，再攻高连接价值卡点”：已经成立的动量映射、正性和同尺度规范相容直接继承；分清对象映射、逻辑蕴涵、同模型共同存在。若只是成熟工具的直接应用，作审计即可，不新增研究轮次；不会等待做完所有容易小题再处理已知矛盾。\n\n573入口是同一反流后来源的局部发展与可用物质参考。先回查548—549及实际读取记录，核Einstein—Yang–Mills—sigma模型的成熟适定性；同一新解的参考雅可比和有效范围才是实际待核接口。反流改变速度，不能沿用旧局部解的满秩结论；初值成立不能替代有限图量子全流的连续极限。\n\n新增M₀²=2、CMC、TT及三维种子等输入明确保留，不把指定二导数作用的相容解称为引力生成。完整势常数保留，550—553裸曲率限制不撤销。停止重复网格、单项预算与参数精扫，整体目标保持。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第572轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 218.' not in body
        body+='\n\n## 218. 共同规范物质与全局Einstein初值\n\n'+block+next_steps
    else:
        assert '\n## 123.' not in body
        body+='\n\n## 123. 全局约束相容仍以维数和作用为输入\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[573共同初值、局部发展与物质参考候选](archive_231_/round573_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为571／2907','最新科学轮次与检查数为572／2915')
        body=body.replace('完成231—571轮。','完成231—572轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第572轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|572|[共同规范物质与全局引力初值](research_note_572.md)|[代码](joint_gauss_einstein_initial_data.py)、'
            '[结果](joint_gauss_einstein_initial_data_results.json)、[核验](research_round_572_checks.json)|\n')
        body+='\n**573候选：** [共同初值、局部发展与物质参考](round573_drafts/STATUS.md)，尚未计轮次；先复用成熟发展接口，再核同一新解的参考秩与有效范围。\n'
    planned[path]=body.replace('\n',newline).encode(encoding)
links=0
for path,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(),(path,link)
        links+=1
for path,raw in raws.items():
    assert path.read_bytes() in (raw,planned[path]),('concurrent edit',path)
for path,body in planned.items():
    if path.read_bytes()==body: continue
    temp=path.with_name(path.name+'.round572.tmp')
    if temp.exists(): assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-10-01',latest_round=572,cumulative_tests=2915,numbered_scientific_files=1028,
    unique_protected_evidence_files=1487,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=573,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
