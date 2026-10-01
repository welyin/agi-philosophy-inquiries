"""Publish round 539 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import rank_two_breaking_bridge as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round539_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_539_checks.json')
assert checked['round']==539 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round539_20260930'
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
summary='**第539轮完成：** [非零轻质量与共同Higgs边界]({p}research_note_539.md)给两活动列的精确二／四次矩与Higgs分段最小值，并用完整Y／Majorana矩阵运行连接两个历史质量差目标。小破缺可实现非零质量，但两例对Higgs下界的改变极小；4λC₅约定已独立核清。8组与独立终审通过，最新539／2695，929份编号科学文件，1185份保护证据；[核验]({p}research_round_539_checks.json)、[条件更新]({p}unified_physics_condition_ledger_539.md)。目标质量是输入，有限匹配、完整味及统一仍开放。'
next_steps='\n\n### 第539轮后：检验共同边界的Higgs低能相容性\n\n539在同一三ν_R内容中引入第二独立微小Yukawa列，证明Q=S²−2|det M_R|ab；简并高能时S≥M(a+b)，并完整给出λ分段最小值。保留rank-one而仅改Majorana不能在树级给两个非零轻质量。破缺后RG生成Majorana对角项，正式计算已保留完整矩阵。\n\n两个历史NuFIT6质量差目标通过矩阵RG与正确4λC₅运行接入；m1=0、v0=246树级转换、带电top-only、共同尺度树匹配均为明确输入。尚非真实PMNS／极点拟合，小阈值劈裂不自动界定相消后轻质量相对误差。\n\n下一正式编号540优先检查同一谱边界是否容许Higgs与top低尺度目标。尽可能以全族比较界替代重复扫描：若一圈＋树匹配分支不可能，明确需要改变的边界或匹配量；随后对接成熟有限一圈匹配。禁止把这个近似分支的反例宣称为全模型、认知原则或统一计划的反证。继续及时编号报告、条件账和核验，不改目标、不设定时任务。\n'
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第539轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 185.' not in body
        body+='\n\n## 185. 非零轻谱与Higgs边界的联立\n\n'+block+next_steps
    else:
        assert '\n## 90.' not in body
        body+='\n\n## 90. 轻质量目标与生成命题的范围\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为538／2687','最新科学轮次与检查数为539／2695')
        body=body.replace('完成231—538轮。','完成231—539轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第539轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|539|[非零轻质量与共同Higgs边界](research_note_539.md)|[代码](rank_two_breaking_bridge.py)、'
            '[结果](rank_two_breaking_bridge_results.json)、[核验](research_round_539_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round539.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=539,cumulative_tests=2695,numbered_scientific_files=929,
    unique_protected_evidence_files=1185,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=540,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
