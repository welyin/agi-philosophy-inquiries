"""Publish round 532 navigation, retaining all earlier text and snapshots."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import shared_trace_unimodularity as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round532_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_532_checks.json')
assert checked['round']==532 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round532_20260930'
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
summary=('**第532轮完成：** [共同权重、群选择与反常]({p}research_note_532.md)'
    '证明不等模块权重可兼容给定有限表示与D／J，但额外要求同Z同时承担unimodularity时，'
    '非零超荷与实际反常抵消强迫等权，不能沿用531不等权匹配。普通群选择＋加权作用分支仍保留；'
    '同一权重的Higgs／Majorana及身份谱矩已同步核算。7组与独立终审通过，最新532／2638，'
    '908份编号科学文件，1137份保护证据；[核验]({p}research_round_532_checks.json)、'
    '[联合条件更新]({p}unified_physics_condition_ledger_532.md)。'
    '只关闭该附加同迹条件下的扩展，不否定全部加权模型或认知统一；三维、GR及整体目标仍开放。')
next_steps='''

### 第532轮后：同时匹配规范、Higgs与几何

532把531动能修补接到群选择和实际手征异常：原代数表示A的Tr ZA=0不同于提升后恒零的Tr ZB；在所列表示及非零超荷下，前者联立异常要求w_l=w_q。此条件U是新增模型输入，不能从“共同认知结构”未经证明地推出。普通unimodularity＋加权谱作用仍通过本轮有限代数检查。

继续研究保留分支：把同一正权的规范动能、Higgs二／四次矩、Majorana项与几何身份迹放入一套归一化和匹配条件，检验共同参数域。直接对接成熟谱作用公式；保留Euclidean到Lorentz、高阶曲率、参考与实际测量、经验检验等开放接口。谱矩在特殊补偿参数下可不变，不能把泛型依赖说成普遍必变。

不重做已知反常族、不重新扫描相同一圈曲线、不为得到预期答案擅加阈值物质。新增可自由参数、边界条件和反例均回填联合条件账；下一正式编号533。完整研究须及时交付research_note_编号.md及复算、结果、核验。目标继续，不改应用目标、不设定时任务、不做图像检验。
'''
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第532轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 178.' not in body
        body+='\n\n## 178. 共同模型中两种迹约束的相容边界\n\n'+block+next_steps
    else:
        assert '\n## 83.' not in body
        body+='\n\n## 83. 联合规范接口仍未选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为531／2631','最新科学轮次与检查数为532／2638')
        body=body.replace('完成231—531轮。','完成231—532轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第532轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|532|[同迹约束的联合边界](research_note_532.md)|[代码](shared_trace_unimodularity.py)、'
            '[结果](shared_trace_unimodularity_results.json)、[核验](research_round_532_checks.json)|\n')
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
    temp=path.with_name(path.name+'.round532.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=532,cumulative_tests=2638,numbered_scientific_files=908,
    unique_protected_evidence_files=1137,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=533,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
