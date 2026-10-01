"""Publish round 551 with existing navigation/history and concurrent-write checks."""
from pathlib import Path
import hashlib
import json
import os
import verify_interaction_rounds as core
import joint_vacuum_hierarchy_matching as model

HERE=Path(__file__).resolve().parent
RESEARCH,ROOT=HERE.parent,HERE.parent.parent
TARGET=HERE/'round551_navigation_checks.json'
assert not TARGET.exists(),'already published'
checked=core.read(HERE/'research_round_551_checks.json')
assert checked['round']==551 and checked['all_reported_checks_passed']
assert model.run()==core.read(model.TARGET)
for group in ('new_file_hashes','preserved_draft_hashes'):
    for name,sha in checked[group].items():
        assert core.digest(HERE/name)==sha,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
       RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
       HERE/'three_dimensional_four_conditions_review.md']
folder=HERE/'navigation_before_round551_20260930'
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

summary='**第551轮完成：** [真空匹配与物质／引力层级]({p}research_note_551.md)分类同一h＋s常场正F解并证明径向线性稳定；固定规范比值与代数目时，仅调真空常数不能在固定弱曲率窗口内得到任意大F/h²。再明示一个共同质量匹配参数，可给平直大层级正面族；新增输入尚未解释。8组、15式及独立终审通过，最新551／2780，965份编号科学文件，1282份保护证据；[核验]({p}research_round_551_checks.json)、[条件更新]({p}unified_physics_condition_ledger_551.md)。常场真空本身不提供满秩参考，量子匹配、完整稳定性与统一仍开放。'
next_steps='''

### 第551轮后：先合并可用条件，再核运行接口

548—549已把同一物质参考与局部几何反作用连接；550暴露共同裸系数的动态曲率限制；551把显式真空与共同质量匹配后的常场可行性、尺度和径向稳定联立。按用户顺序优先继承已证连接，不以单部门无限精细优化替代条件压缩。

下一候选552检查同一个质量匹配比例与544已有一圈运行的接口。整体齐次缩放及共同初值分裂已有依据，不能重报新轮次；需要新的联合约束或共同可行证书才正式编号。曲率耦合、真空项、真实量子读取及完整作用近似仍需分别核，不能用物质子系统的成功签收整体。

统一目标与冻结历史保持；阶段尚未完成。
'''
planned={}
for path,raw in raws.items():
    encoding='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(encoding).replace('\r\n','\n')
    assert '**第551轮完成：**' not in body
    prefix='research_cognition_physics/archive_231_/' if path.parent==ROOT else 'archive_231_/' if path.parent==RESEARCH else ''
    block=summary.format(p=prefix)
    if path in paths[:5]:
        title,tail=body.split('\n\n',1); body=title+'\n\n'+block+'\n\n'+tail
    elif path.name=='spatial_premise_closure_audit.md':
        assert '\n## 197.' not in body
        body+='\n\n## 197. 真空、物质和引力系数的共同匹配\n\n'+block+next_steps
    else:
        assert '\n## 102.' not in body
        body+='\n\n## 102. 尺度共同可行仍不选择空间维数\n\n'+block
    if path.name in ('RESEARCH_STATE.md','research_direction.md'):
        body+=next_steps
        body+='\n[552共同质量与运行接口草稿](archive_231_/round552_drafts/STATUS.md)已登记，未计完成轮次。\n'
    if path.name=='research_direction.md':
        body=body.replace('最新科学轮次与检查数为550／2772','最新科学轮次与检查数为551／2780')
        body=body.replace('完成231—550轮。','完成231—551轮。')
    if path==HERE/'README.md':
        body+=('\n\n## 第551轮研究索引\n\n|轮次|完整报告|复算与核验|\n|---|---|---|\n'
            '|551|[真空匹配与物质／引力层级](research_note_551.md)|[代码](joint_vacuum_hierarchy_matching.py)、'
            '[结果](joint_vacuum_hierarchy_matching_results.json)、[核验](research_round_551_checks.json)|\n')
        body+='\n**552候选：** [共同质量与运行接口](round552_drafts/STATUS.md)，尚未完成；已有齐次缩放推论不另计轮次。\n'
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
    temp=path.with_name(path.name+'.round551.tmp')
    if temp.exists():
        assert temp.read_bytes()==body
    else:
        with temp.open('xb') as stream: stream.write(body)
    assert path.read_bytes()==raws[path]
    os.replace(temp,path)
report=dict(date='2026-09-30',latest_round=551,cumulative_tests=2780,numbered_scientific_files=965,
    unique_protected_evidence_files=1282,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,recovered_partial_navigation_publish=resuming,
    complete_report_published=True,existing_results_preserved=True,next_round=552,
    stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:
    stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
