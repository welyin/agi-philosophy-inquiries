"""Publish the 1308 goal as one reversible navigation layer."""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import relocation_original_bytes

def sha(raw):return hashlib.sha256(raw).hexdigest()

def main():
    out=HERE/'layer.json'
    if out.exists():raise SystemExit('Navigation already recorded.')
    mp=BASE/'_migration/closure_1288_1307/manifest.json'
    closure=json.loads(mp.read_text(encoding='utf-8'))
    boundary={e['destination']:e for e in closure['entries']}
    prefix=('\n## 当前阶段：群、代际与参数关系的条件性生成\n\n'
        '2026-10-18，用户授权新目标，应用目标已启用。研究进入[archive_1308_](archive_1308_/README.md)，'
        '[正式目标与验收](archive_1308_/_shared/goal.md)、[范围审阅](archive_1308_/_shared/goal_scope_review.md)和'
        '[复用／首项筛查](archive_1308_/_shared/start_plan.md)已保存。最新正式科学成果仍是1307；启动与查重计0。\n\n'
        '目标是让规范群／表示、代际数和参数关系成为更一般、具有独立认知或操作含义的条件下的推论，'
        '并严格判清剩余自由。暂不直接公理化标准群、恰三代或现实参数，不默认最小性。'
        '三条线须有明确结论并检验共同分支非空；单轮或单一候选反例不作为整体结项。\n\n'
        '当前先审开放颜色族与跨部门交换的新增选择力，不重做旧色数反常、n≥3相位或更多整数样本。'
        '动态场论只在影响选择／有限预测时补充。两旧阶段已正式收尾，迁移与复算见'
        '[归档记录](_migration/closure_1288_1307/README.md)。以下条目保留原时点状态，当前方向以本节和新目标为准。\n')
    updates=[]
    for filename in ('README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'):
        path=BASE/filename;raw=path.read_bytes();name=path.relative_to(ROOT).as_posix()
        assert sha(raw)==boundary[name]['current_sha256'],name
        enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
        text=raw.decode(enc);at=text.index('\n')+1
        current=(text[:at]+prefix+text[at:]).encode(enc)
        entry=dict(original=name,destination=name,original_sha256=sha(raw),current_sha256=sha(current),
            original_bytes=len(raw),current_bytes=len(current),edits=[dict(start=at,end=at,before='',after=prefix)])
        assert relocation_original_bytes(current,entry)==raw
        updates.append((path,raw,current,entry))
    support=dict(closure['unchanged_support_sha256'])
    for p in (mp,ROOT/'scripts/run_research_structure_selection_closed.py'):
        support[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    for name,expected in support.items():assert sha((ROOT/name).read_bytes())==expected,name
    stage='research_cognition_physics/archive_1308_'
    entry_points=[stage+'/README.md',stage+'/_shared/goal.md']
    assert all((ROOT/name).is_file() for name in entry_points)
    layer=dict(schema='conditional_selection_active_navigation_v1',date='2026-10-18',stage_path=stage,
        stage_first_round=1308,entry_points=entry_points,latest_accepted_round=1307,in_progress_round=None,
        application_goal_status_last_observed='active',startup_state='admission screening; no new scientific acceptance',
        entries=[e for _,_,_,e in updates],unchanged_support_sha256=support)
    assert all(p.read_bytes()==raw for p,raw,_,_ in updates)
    for p,_,current,_ in updates:p.write_bytes(current)
    with out.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(layer,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':'ACTIVE_NAVIGATION_CREATED','entries':4,'support':len(support),
                      'layer_sha256':sha(out.read_bytes()),'new_scientific_rounds':0}))

if __name__=='__main__':main()
