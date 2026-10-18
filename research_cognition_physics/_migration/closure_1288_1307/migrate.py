"""Record reversible link/navigation edits; native PowerShell performs renames."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'research_cognition_physics'
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes

STAGES = [dict(source='research_cognition_physics/archive_1288_',
               destination='research_cognition_physics/archive_1288_1289', first_round=1288,last_round=1289),
          dict(source='research_cognition_physics/archive_1290_',
               destination='research_cognition_physics/archive_1290_1307', first_round=1290,last_round=1307)]
NAV = {BASE/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md')}
PATTERN = re.compile(r'archive_(1288|1290)_(?![0-9])')
REPLACEMENT = {'1288':'archive_1288_1289','1290':'archive_1290_1307'}

def sha(raw): return hashlib.sha256(raw).hexdigest()
def dump(path,data):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
def edited(raw,edits):
    if not edits: return raw
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    text=raw.decode(enc)
    for e in reversed(edits):
        assert text[e['start']:e['end']]==e['before']
        text=text[:e['start']]+e['after']+text[e['end']:]
    return text.encode(enc)

def prepare():
    before=json.loads((HERE/'before.json').read_text(encoding='utf-8'))
    files={ROOT/e['path'] for e in before['files']}
    for e in before['files']:
        raw=(ROOT/e['path']).read_bytes()
        assert sha(raw)==e['sha256'] and len(raw)==e['bytes'],e['path']
    for s in STAGES:
        assert (ROOT/s['source']).is_dir() and not (ROOT/s['destination']).exists()
    external=set(NAV)
    for p in BASE.rglob('*.md'):
        if p.is_relative_to(HERE) or p in files: continue
        if PATTERN.search(p.read_text(encoding='utf-8-sig')): external.add(p)
    entries=[]
    for p in sorted(files|external):
        raw=p.read_bytes();name=p.relative_to(ROOT).as_posix();dest=name
        for s in STAGES:
            if name.startswith(s['source']+'/'):dest=s['destination']+name[len(s['source']):]
        edits=[]
        if p.suffix=='.md':
            text=raw.decode('utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8')
            for m in PATTERN.finditer(text):
                edits.append(dict(start=m.start(),end=m.end(),before=m.group(),after=REPLACEMENT[m.group(1)]))
            at=text.index('\n')+1
            addition=''
            if p in NAV:
                addition=('\n## 1288—1307两阶段正式归档\n\n'
                    '2026-10-18。两阶段已按原验收范围收尾：'
                    '[1288—1289](archive_1288_1289/README.md)对应论文12，'
                    '[1290—1307](archive_1290_1307/README.md)对应论文13。'
                    '20轮报告、代码、结果及验收全部保留；新增科学轮次与公理均为0。'
                    '[迁移证据与历史复算](_migration/closure_1288_1307/README.md)。\n\n'
                    '后继从1308开始，研究让规范群、代际数及参数关系成为具有独立认知／操作含义的条件下的推论。'
                    '用户明确暂不直接公理化这些输出；不默认自然界选最少物种或最少代际。'
                    '本节以下状态为各时点历史，后继目标以新阶段正式记录为准。\n')
            for s in STAGES:
                if name==s['source']+'/README.md':
                    title=text[:at].rstrip('\r\n')
                    if not title.startswith(f"# {s['first_round']}—{s['last_round']}"):
                        edits.append(dict(start=0,end=len(title),before=title,
                            after=f"# {s['first_round']}—{s['last_round']}：已结项研究档案"))
                    addition=('\n## 正式归档与后继\n\n'
                        f"2026-10-18。目录收尾为 `{Path(s['destination']).name}`。"
                        '本阶段原验收已完成，以下待审／继续表述按其原时点保留；不改变冻结的科学结论。'
                        '代码与JSON逐字节保留，Markdown仅更新路径与归档导航。'
                        '[迁移证据与当前复算入口](../_migration/closure_1288_1307/README.md)。\n\n'
                        '后继由1308起研究群、代际及参数关系的条件性选择；'
                        '本次归档本身不计新研究轮次，也不采用新公理。\n')
            if addition:edits.append(dict(start=at,end=at,before='',after=addition))
        edits.sort(key=lambda e:(e['start'],e['end']))
        current=edited(raw,edits)
        entry=dict(original=name,destination=dest,original_sha256=sha(raw),current_sha256=sha(current),
                   original_bytes=len(raw),current_bytes=len(current),edits=edits)
        assert relocation_original_bytes(current,entry)==raw
        entries.append(entry)
    layer=BASE/'_migration/structure_selection_active_1288/layer.json'
    support=dict(json.loads(layer.read_text(encoding='utf-8'))['unchanged_support_sha256'])
    for p in (layer,ROOT/'scripts/run_research_structure_selection_active.py'):
        support[p.relative_to(ROOT).as_posix()]=sha(p.read_bytes())
    for name,expected in support.items():assert sha((ROOT/name).read_bytes())==expected,name
    plan=dict(schema='research_closure_1288_1307_v1',date='2026-10-18',stages=STAGES,entries=entries,
              unchanged_support_sha256=support,moved_files=len(files),external_markdown_files=len(external),
              directory_copy_used=False,zip_used=False,scientific_results_rewritten=False)
    dump(HERE/'plan.json',plan)
    print(json.dumps({'status':'PLAN_READY','moved_files':len(files),'external_files':len(external),
                      'edited_markdown':sum(bool(e['edits']) for e in entries)}))

def apply():
    plan=json.loads((HERE/'plan.json').read_text(encoding='utf-8'))
    for s in STAGES:assert not (ROOT/s['source']).exists() and (ROOT/s['destination']).is_dir()
    for e in plan['entries']:assert sha((ROOT/e['destination']).read_bytes()) in (e['original_sha256'],e['current_sha256']),e['destination']
    for e in plan['entries']:
        p=ROOT/e['destination'];raw=p.read_bytes()
        if sha(raw)!=e['current_sha256']:p.write_bytes(edited(raw,e['edits']))
        relocation_original_bytes(p.read_bytes(),e)
    dump(HERE/'manifest.json',plan)
    assert (HERE/'plan.json').read_bytes()==(HERE/'manifest.json').read_bytes()
    (HERE/'plan.json').unlink()
    print(json.dumps({'status':'MIGRATED','files':len(plan['entries']),'manifest_sha256':sha((HERE/'manifest.json').read_bytes())}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['prepare','apply'])
    prepare() if parser.parse_args().mode=='prepare' else apply()
