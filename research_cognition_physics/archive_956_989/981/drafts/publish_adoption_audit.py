"""Add one scope-audit pointer to living indexes without advancing scientific rounds."""
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
paths += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md',
                           '跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
old={p:p.read_bytes() for p in paths}
updates={}
heading='## 980后采用复核：先合并父模型，停止局部扩建'
anchor='## 980：有限热材料与有效遗忘'
for p,raw in old.items():
    text=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert heading not in text and text.count(anchor)==1
    prefix='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'
        +f'[981工作报告]({prefix}981/drafts/research_note_981_working.md)核对共同父模型与有限有效域：'
        +'复用已有材料、表示、来源及热结果，停止局部器件扩建。下一项先冻结父描述的物种、作用阶数、匹配、参考和共同任务，再补决定采用的接口。'
        +'本次是旧证据与投入审计，新增科学试验组0；正式980／累计3765，981科学轮尚未结项。'
        +f'[审计核验]({prefix}981/drafts/common_adoption_audit_checks.json)。应用目标及957假说保持，整体目标未完成。\n\n')
    text=text.replace(anchor,block+anchor,1)
    if b'\r\n' in raw:
        text=text.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+text.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in old.items())
for p,raw in updates.items():
    p.write_bytes(raw)
print('Updated eight living navigation entries; formal round remains 980.')
