"""Publish round 542 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import protected_pair_dim6_feedback as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round542_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_542_checks.json')
assert checked['round']==542 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round542_20260930'
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
summary='**第542轮完成：** [质量符号与维六反馈]({p}research_note_542.md)证明同一保护配对的负质量目标和负系数迹，使树级维六单插入反馈沿低能方向非负；所列SEFT子系统仍有λ>0.16259。迹保持负号而矩阵负半定性不保持，已给反例。7组与独立终审通过，最新542／2716，938份编号科学文件，1207份保护证据；[核验]({p}research_round_542_checks.json)、[条件更新]({p}unified_physics_condition_ledger_542.md)。完整EFT观测、一般破缺及统一仍开放。'
next_steps='\n\n### 第542轮后：检验能够改变共同边界的真实结构修订\n\n542接入树级Higgs—轻子流算符的直接一圈反馈。负低能m²目标和负系数迹分别由乘法方程保持，向下λ新增源非负，不能供应540负修订预算；结论对声明的闭合SEFT子系统成立，不升级为完整SMEFT重求和或实际观测提取。一般非简并破缺及高能谱质量边界仍未关闭。\n\n下一正式编号543对接成熟中性标量／Higgs机制（Chamseddine—Connes 1208.1030），先检验其在现有有限代数、表示及内涨落中的来源，再联立Majorana、标量稳定性与共同参数。可能需要改变代数／允许操作或添加自由度，须明确列账；不能只手调负阈值。先查相关历史和成熟来源，不把复述文献算新轮次，不继续对已受全族障碍覆盖的精确保护对微调。\n\n全局C编号继续遵循首版：质量和参数=C17—C18、运行=C20、观测=C27。几何、维数、Lorentz、实际记录、资源与全尺度统一仍开放。目标保持不变，每项实质结果及时形成编号报告。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第542轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 188.' not in body
        body+='\n\n## 188. 质量与运行的联合方向限制\n\n'+block+next_steps
    else:
        assert '\n## 93.' not in body
        body+='\n\n## 93. 单插入有效运行的适用边界\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为541／2709','最新科学轮次与检查数为542／2716')
        body=body.replace('完成231—541轮。','完成231—542轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第542轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|542|[质量符号与维六反馈](research_note_542.md)|[代码](protected_pair_dim6_feedback.py)、'
            '[结果](protected_pair_dim6_feedback_results.json)、[核验](research_round_542_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round542.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=542,cumulative_tests=2716,numbered_scientific_files=938,
    unique_protected_evidence_files=1207,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=543,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
