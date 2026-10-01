"""Publish round 540 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_higgs_flow_obstruction as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round540_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_540_checks.json')
assert checked['round']==540 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round540_20260930'
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
        with (folder/key).open('xb') as stream:
            stream.write(raw)
        manifest[str(p.relative_to(ROOT))]=dict(snapshot=key,sha256=hashlib.sha256(raw).hexdigest())
    with (folder/'manifest.json').open('x',encoding='utf8') as stream:
        json.dump(manifest,stream,ensure_ascii=False,indent=2)
summary='**第540轮完成：** [共同Higgs边界的全族运行障碍]({p}research_note_540.md)用完成平方、矩阵行列式增长界和精确有理数证书，证明所列两活动ν小破缺、一圈＋树匹配分支的全部top／阈值选择都有λ>0.17009，不能达到历史MS目标0.12604。保持比较结构的负修订预算至少881/30000。7组与独立终审通过，最新540／2702，932份编号科学文件，1193份保护证据；[核验]({p}research_round_540_checks.json)、[条件更新]({p}unified_physics_condition_ledger_540.md)。只排除明确近似边值问题，有限匹配及完整统一仍开放。'
next_steps='\n\n### 第540轮后：以真实有限匹配检验修订预算\n\n540对冻结规范、top-only带电、至多两活动ν列、高能D0=det(Y†Y)≤10^-6及λ树级连续匹配的分支，证明全族低能λ>0.17009。完成平方和精确100步有理下包络替代top／阈值扫描；保持同一比较结构而达到历史MSλ=.12604，负初值、跃变与连续源总预算须至少881/30000。该必要条件有明确正性bootstrap，不是任意新理论的充分修补方案。\n\n下一正式编号541对接成熟type-I seesaw有限一圈匹配，核同一复相关强耦合及小破缺分支的λ阈值符号和量级，再与预算比较；需要同时交代Higgs二次项、场规范化和其它部门。不得自由手调负跃变后称其已由模型产生，不继续扫描已排除无修正分支。\n\n历史MS中央值不是当前精密拟合；全量物理模型、高圈、其它物质或不同边界尚未排除。共同几何、实际记录、内部资源、维数及Lorentz来源仍开放。每项实质结果及时编号报告并回填总账，不无故更改目标，不设定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第540轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 186.' not in body
        body+='\n\n## 186. 共同参数族的全尺度限制\n\n'+block+next_steps
    else:
        assert '\n## 91.' not in body
        body+='\n\n## 91. 全族反例与整体生成命题的范围\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为539／2695','最新科学轮次与检查数为540／2702')
        body=body.replace('完成231—539轮。','完成231—540轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第540轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|540|[共同Higgs边界的全族运行障碍](research_note_540.md)|[代码](joint_higgs_flow_obstruction.py)、'
            '[结果](joint_higgs_flow_obstruction_results.json)、[核验](research_round_540_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round540.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=540,cumulative_tests=2702,numbered_scientific_files=932,
    unique_protected_evidence_files=1193,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=541,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
