"""Close one research stage; move files once and record reversible Markdown edits."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess

from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'research_cognition_physics'
SOURCE = BASE/'archive_1234_'
TARGET = BASE/'archive_1234_1280'
MIG = BASE/'_migration/closure_1234_1280'
MANIFEST = MIG/'manifest.json'
RUNNER = 'scripts/run_research_gr_response_closed.py'
NAV = {'README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'}
LINK = re.compile(r'(?P<lead>\[[^\]\n]*\]\()(?P<target>[^)\n]+)(?P<tail>\))')

def sha(raw): return hashlib.sha256(raw).hexdigest()
def rel(p): return p.relative_to(ROOT).as_posix()
def mapped(p): return TARGET/p.relative_to(SOURCE) if p.is_relative_to(SOURCE) else p
def rename_text(s): return re.sub(r'archive_1234_(?!\d)', 'archive_1234_1280', s)

def notice():
    return ('## 1234—1280阶段归档\n\n'
      '已归档为[archive_1234_1280](archive_1234_1280/README.md)，保持一个完整阶段。'
      '47轮共同服务于[论文10](10_认知操作公理与有效广义相对论_阶段论文.md)，'
      '阶段内按五段研究进程导航；报告、代码、结果及原验收全部保留。'
      '本目标已结项，未启动新目标或1281轮。'
      '[目录核验与复算](_migration/closure_1234_1280/README.md)。\n\n')

def overview():
    return ('## 阅读导航\n\n'
      '本目录是一个已结项阶段：**共同几何响应与有效广义相对论重建**。'
      '下列五段是同一研究链的阅读分组，不是五个独立结项阶段。'
      '它们共享父模型、来源账和最终论文，因此不拆分或复制共享材料。\n\n'
      '|轮次|研究进程|报告入口|\n|---|---|---|\n'
      '|1234—1244|热角色边界、CT与同源几何接口|[1234](research_note_1234.md) → [1244](research_note_1244.md)|\n'
      '|1245—1252|动态几何、记录、量子桥与有限预测|[1245](research_note_1245.md) → [1252](research_note_1252.md)|\n'
      '|1253—1270|原生量子材料、局部来源与因果读口|[1253](research_note_1253.md) → [1270](research_note_1270.md)|\n'
      '|1271—1278|治理机制、热交换、共识边界与参考选择|[1271](research_note_1271.md) → [1278](research_note_1278.md)|\n'
      '|1279—1280及收尾|公理采用后的GR重建与共同有效验收|[1279](research_note_1279.md) → [1280](research_note_1280.md) → [结项](_shared/jacobson_stage_closure_20261017.md)|\n\n'
      '**先看结论：**[阶段论文10](../10_认知操作公理与有效广义相对论_阶段论文.md) · '
      '[逐轮文件索引](文件索引.md) · [结项审计](_shared/jacobson_stage_closure_20261017.md)。\n\n'
      '**目录格式：**`research_note_1234.md`等报告在本目录根部；各轮代码、结果、证明和验收在同编号子目录；'
      '`_shared`保存跨轮共同稿与账本，`_admission`保存候选准入和非编号接续。共享内容只维护一份。\n\n'
      '**复算：**使用[归档只读入口](../_migration/closure_1234_1280/README.md)。'
      '科学Python/JSON字节保持不变；旧路径及旧文档哈希通过记录的映射在内存中恢复。'
      '不需要ZIP、旧目录副本或路径替身。下方记录保留当时状态，当前结论以结项审计为准。\n\n')

def rewrite(p, raw):
    if p.suffix.lower() != '.md': return raw, []
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
    before = raw.decode(enc)
    after = LINK.sub(lambda m:m['lead']+rename_text(m['target'])+m['tail'],before)
    if p.parent == BASE and p.name in NAV:
        head,tail=after.split('\n',1)
        after=head+'\n\n'+notice()+tail.lstrip('\r\n')
    if p.parent in (SOURCE,TARGET) and p.name=='README.md':
        _,tail=after.split('\n',1)
        after='# 1234—1280：共同几何响应与有效广义相对论（研究档案）\n\n'+overview()+tail.lstrip('\r\n')
    if p.parent in (SOURCE,TARGET) and p.name=='文件索引.md':
        _,tail=after.split('\n',1)
        after=('# 1234—1280轮文件索引\n\n'
          '[主题导航](README.md) · [归档核验与只读复算](../_migration/closure_1234_1280/README.md)。'
          '以下原命令、路径与研究状态保留历史含义；当前文件归属于本完整编号目录。\n\n'+tail.lstrip('\r\n'))
    edits=text_edits(before,after)
    for e in edits:e['kind']='archive_link_or_closure_navigation'
    return after.encode(enc),edits

def guards():
    for p in (SOURCE,TARGET,MIG):
        absolute=p.resolve()
        if not absolute.is_relative_to(ROOT.resolve()) or absolute in (ROOT.resolve(),BASE.resolve()):
            raise ValueError('Unsafe resolved target: '+str(absolute))
    if SOURCE.parent.resolve()!=BASE.resolve() or TARGET.parent.resolve()!=BASE.resolve():
        raise ValueError('Only the explicitly named sibling archive may move')

def plan():
    guards()
    if MANIFEST.exists() or TARGET.exists() or not SOURCE.is_dir():
        raise ValueError('Unexpected source, target, or previous migration')
    receipt=json.loads((SOURCE/'_shared/jacobson_goal_completion_20261017.json').read_text(encoding='utf-8'))
    if receipt['goal']['status']!='complete':raise ValueError('Research goal not completed')
    expected=set(range(1234,1281))
    actual={int(p.stem.split('_')[-1]) for p in SOURCE.glob('research_note_*.md')}
    if actual!=expected:raise ValueError('Round-note inventory mismatch')
    if any(not (SOURCE/str(n)).is_dir() for n in expected):raise ValueError('Missing round directory')
    moved=[]
    for p in [SOURCE,*SOURCE.rglob('*')]:
        if p.is_symlink() or getattr(p,'is_junction',lambda:False)():raise ValueError('Reparse point: '+str(p))
        if p.is_file():moved.append(p)
    if any(p.suffix.lower() in {'.zip','.pyc','.pyo'} for p in moved):
        raise ValueError('Unexpected ZIP/bytecode; do not silently remove')
    search=subprocess.run(['rg','-l','--hidden','-g','*.md','-g','!**/.git/**',
        '-g','!**/.research_runtime/**','-g','!**/.tmp/**','-g','!**/closure_1234_1280/**',
        'archive_1234_', '.'],cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    if search.returncode not in (0,1):raise RuntimeError(search.stderr)
    candidates={ROOT/s.removeprefix('.\\').removeprefix('./') for s in search.stdout.splitlines()}
    candidates.update(BASE/n for n in NAV)
    external=[p for p in candidates if p not in moved and rewrite(p,p.read_bytes())[0]!=p.read_bytes()]
    entries=[]
    for p in sorted(set(moved+external)):
        raw=p.read_bytes(); current,edits=rewrite(p,raw)
        e=dict(original=rel(p),destination=rel(mapped(p)),original_sha256=sha(raw),current_sha256=sha(current),
               original_bytes=len(raw),current_bytes=len(current),edits=edits)
        assert relocation_original_bytes(current,e)==raw
        entries.append(e)
    support=[p for p in (ROOT/'scripts').glob('*.py') if p.name not in {Path(__file__).name,Path(RUNNER).name}]
    support += [p for p in (BASE/'_migration').rglob('*') if p.is_file() and not p.is_relative_to(MIG)
                and p.suffix.lower() in {'.json','.py'}]
    return dict(schema='research_closure_1234_1280_v1',date='2026-10-17',
        stages=[dict(source=rel(SOURCE),destination=rel(TARGET),first_round=1234,last_round=1280)],
        entries=entries,unchanged_support_sha256={rel(p):sha(p.read_bytes()) for p in sorted(set(support))},
        source_file_counts={'archive_1234_':len(moved)},source_suffix_counts=dict(Counter(p.suffix for p in moved)),
        external_markdown_files=len(external),scientific_python_json_unchanged=True,new_scientific_rounds=0,
        goal_status='complete',next_round=1281,new_goal_started=False,directory_copy_used=False,
        old_directory_shell=False,zip_created=False,replay_entry=RUNNER)

def summary(d):
    return {k:v for k,v in d.items() if k not in {'entries','unchanged_support_sha256'}} | dict(
        manifest_entries=len(d['entries']),reversible_markdown_files=sum(bool(e['edits']) for e in d['entries']),
        frozen_support_files=len(d['unchanged_support_sha256']))

def apply():
    d=plan()
    for e in d['entries']:assert sha((ROOT/e['original']).read_bytes())==e['original_sha256']
    MIG.mkdir(parents=True,exist_ok=True)
    with MANIFEST.open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
    renamed=False;patched=[]
    try:
        SOURCE.rename(TARGET);renamed=True
        for e in d['entries']:
            if not e['edits']:continue
            p=ROOT/e['destination'];raw=p.read_bytes()
            assert sha(raw)==e['original_sha256']
            new,edits=rewrite(p,raw)
            assert edits==e['edits'] and sha(new)==e['current_sha256']
            p.write_bytes(new);patched.append((p,raw,new))
    except Exception:
        for p,raw,new in reversed(patched):
            if p.read_bytes()!=new:raise RuntimeError('Concurrent edit prevents rollback: '+str(p))
            p.write_bytes(raw)
        if renamed:TARGET.rename(SOURCE)
        MANIFEST.unlink()
        raise
    print(json.dumps(summary(d)|dict(status='MOVED'),ensure_ascii=False))

def verify():
    guards();d=json.loads(MANIFEST.read_text(encoding='utf-8'))
    for e in d['entries']:
        raw=(ROOT/e['destination']).read_bytes();prior=relocation_original_bytes(raw,e)
        if Path(e['destination']).suffix.lower() in {'.py','.json'} and raw!=prior:raise ValueError('Scientific bytes changed')
    for p,h in d['unchanged_support_sha256'].items():assert sha((ROOT/p).read_bytes())==h,p
    expected={e['destination'] for e in d['entries'] if e['original'].startswith(rel(SOURCE)+'/')}
    actual={rel(p) for p in TARGET.rglob('*') if p.is_file()}
    assert not SOURCE.exists() and actual==expected,'Physical inventory mismatch'
    print(json.dumps(summary(d)|dict(status='PASS',old_directory_absent=True),ensure_ascii=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group()
    g.add_argument('--apply',action='store_true');g.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    if args.apply:apply()
    elif args.verify:verify()
    else:print(json.dumps(summary(plan()),ensure_ascii=False))
