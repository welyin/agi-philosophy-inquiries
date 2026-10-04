"""Save735 working evidence without increasing formal research counts."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
RESEARCH = ARCHIVE.parent
ROOT = RESEARCH.parent
sys.path.insert(0, str(ARCHIVE))
import verify_interaction_rounds as core
import gauge_history_scope_probe as probe

TARGET = HERE / 'followup_checks.json'
assert not TARGET.exists(), 'Published followup is immutable; use its receipt for read-only checks.'
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', ARCHIVE/'README.md']
before = {path: path.read_bytes() for path in paths}

result = json.loads(json.dumps(probe.run()))
assert result == core.read(probe.TARGET)
assert not result['completed_new_round'] and result['latest_completed_round'] == 734
for name, digest in result['dependencies'].items():
    assert core.digest(ARCHIVE/name) == digest, name

history = dict(core.read(ARCHIVE/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for number in range(584, 735):
    receipt = core.read(ARCHIVE/f'research_round_{number}_checks.json')
    for key in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, digest in receipt[key].items():
            assert name not in history or history[name] == digest
            history[name] = digest
history.update(core.read(ARCHIVE/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history) == 3284
for name, digest in history.items():
    assert core.digest(ARCHIVE/name) == digest, name
entry = core.read(HERE/'entry_checks.json')
for name, digest in entry['artifact_hashes'].items():
    assert core.digest(HERE/name) == digest, name

texts = [core.text_checks(HERE/'research_note_735_working.md'),
         core.text_checks(HERE/'ward_literature_scope_review.md')]
assert [item['display_formulas'] for item in texts] == [5, 1]
for name in ('gauge_history_scope_probe.py', Path(__file__).name):
    ast.parse((HERE/name).read_text('utf8'))
links = 0
for name in ('research_note_735_working.md', 'ward_literature_scope_review.md'):
    for link in core.link_parser()((HERE/name).read_text('utf8')):
        destination = (HERE/link).resolve()
        assert destination.exists() or destination == TARGET.resolve(), (name, link)
        links += 1

heading = '**735工作报告已保存，正式仍734／3407：**'
summary = (heading + ' [实际规范历史与来源规范化]( {p}round735_drafts/research_note_735_working.md)'
           '原完整含时规范字典、同一过去参考与记录通过三项补充校准；保留时间连接是必要条件。'
           '态／记录协变不能独自认证绝对来源；相关文献的非零物质背景与Yukawa范围已核。'
           '下一项仍是原联合Wick余项及同源局部修正，不重算荷表或重复质量坐标，未新增完成轮次。')
summary = summary.replace(']( ', '](')
planned = {}
for path, raw in before.items():
    encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline = '\r\n' if b'\r\n' in raw else '\n'
    text = raw.decode(encoding).replace('\r\n', '\n')
    assert heading not in text
    prefix = 'research_cognition_physics/archive_231_/' if path.parent == ROOT else 'archive_231_/' if path.parent == RESEARCH else ''
    head, rest = text.split('\n\n', 1)
    planned[path] = (head+'\n\n'+summary.format(p=prefix)+'\n\n'+rest).replace('\n', newline).encode(encoding)

navigation_links = 0
for path, raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (path.parent/link).resolve().exists(), (path, link)
        navigation_links += 1
assert all(path.read_bytes() == raw for path, raw in before.items()), 'Concurrent navigation edit'
folder = ARCHIVE/'navigation_before_round735_followup_20261004'
folder.mkdir(exist_ok=False)
manifest = {}
for index, (path, raw) in enumerate(before.items()):
    name = f'{index}_{path.name}'
    with (folder/name).open('xb') as handle:
        handle.write(raw)
    manifest[path.relative_to(ROOT).as_posix()] = dict(snapshot=name, sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x', encoding='utf8') as handle:
    json.dump(manifest, handle, ensure_ascii=False, indent=2)
for path, raw in planned.items():
    assert path.read_bytes() == before[path]
    temporary = path.with_name(path.name+'.round735followup.tmp')
    with temporary.open('xb') as handle:
        handle.write(raw)
    os.replace(temporary, path)
for name, digest in history.items():
    assert core.digest(ARCHIVE/name) == digest, name
assert all(path.read_bytes() == raw for path, raw in planned.items())

artifacts = ('gauge_history_scope_probe.py', 'gauge_history_scope_probe_results.json',
             'ward_literature_scope_review.md', 'research_note_735_working.md', Path(__file__).name)
receipt = dict(date='2026-10-04', round_in_progress=735, latest_completed_round=734,
               formal_test_count_unchanged=3407, completed_new_round=False,
               supplementary_probe_checks=3, saved_results_reproduced=True,
               previous_protected_artifacts=3284, previous_hashes_unchanged=True,
               entry_artifacts_unchanged=True, text_checks=texts, local_report_links=links,
               navigation_files=5, navigation_links=navigation_links, broken_links=0,
               navigation_snapshots_preserved=True, active_goal_unchanged=True,
               no_scheduled_task_created=True, no_visual_check=True,
               artifact_hashes={name: core.digest(HERE/name) for name in artifacts},
               navigation_hashes={path.relative_to(ROOT).as_posix(): core.digest(path) for path in paths},
               all_checks_passed=True)
with TARGET.open('x', encoding='utf8', newline='\n') as handle:
    handle.write(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({key: receipt[key] for key in ('round_in_progress', 'latest_completed_round',
                 'navigation_links', 'previous_protected_artifacts', 'all_checks_passed')}))
