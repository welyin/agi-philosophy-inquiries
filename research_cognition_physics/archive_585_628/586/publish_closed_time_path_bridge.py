"""Publish a method review, without modifying frozen rounds or their counts."""
import ast
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
RESEARCH, ROOT = HERE.parent, HERE.parent.parent
SNAPSHOT = HERE/'navigation_before_closed_time_path_bridge_20261001'
TARGET = HERE/'closed_time_path_bridge_review_checks.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert not SNAPSHOT.exists() and not TARGET.exists()
    saved = json.loads((HERE/'closed_time_path_bridge_entry_results.json').read_text('utf8'))
    assert saved['checks_passed'] and saved['cumulative_numbered_checks_unchanged'] == 2989
    for name, expected in saved['source_hashes'].items():
        assert sha(HERE/name) == expected, name
    old = json.loads((HERE/'research_round_586_checks.json').read_text('utf8'))
    for key in ('new_file_hashes', 'preserved_draft_hashes'):
        for name, expected in old[key].items():
            assert sha(HERE/name) == expected, name
    for filename in ('closed_time_path_bridge_entry.py', 'publish_closed_time_path_bridge.py'):
        ast.parse((HERE/filename).read_text('utf8'))
    review = HERE/'closed_time_path_bridge_review_586.md'
    report_text = review.read_text('utf8')
    assert report_text.count('$$') % 2 == 0
    assert '586／2989' in report_text and '不新增编号轮次' in report_text
    paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
             RESEARCH/'RESEARCH_STATE.md', HERE/'README.md']
    before = {p: p.read_bytes() for p in paths}
    planned = {}
    marker = '**路径积分方法对接（2026-10-01，不增轮次）：**'
    for p, raw in before.items():
        prefix = 'research_cognition_physics/archive_231_/' if p.parent == ROOT else 'archive_231_/' if p.parent == RESEARCH else ''
        enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
        newline = '\r\n' if b'\r\n' in raw else '\n'
        body = raw.decode(enc).replace('\r\n', '\n')
        assert marker not in body
        head, rest = body.split('\n\n', 1)
        addition = (marker+f' [闭合时间路径与共同量子过程]({prefix}closed_time_path_bridge_review_586.md)'
                    '将586原记录仪器接为双分支核：同报告概率下的非对角差重现既有源能源差。'
                    '采用CTP／影响泛函／过程张量对接记录、反作用和尺度，'
                    f'[入口复算]({prefix}closed_time_path_bridge_entry_results.json)通过。'
                    '科学基线仍586／2989；587能源流候选保留，动态几何、主体自主实现及连续极限仍开放。')
        planned[p] = (head+'\n\n'+addition+'\n\n'+rest).replace('\n', newline).encode(enc)
    links = 0
    for p, raw in {**planned, review: review.read_bytes()}.items():
        for link in core.link_parser()(raw.decode('utf-8-sig')):
            assert (p.parent/link).resolve().exists(), (p, link)
            links += 1
    assert all(p.read_bytes() == raw for p, raw in before.items()), 'concurrent navigation change'
    SNAPSHOT.mkdir(exist_ok=False)
    manifest = {}
    for i, (p, raw) in enumerate(before.items()):
        filename = f'{i}_{p.name}'
        with (SNAPSHOT/filename).open('xb') as f:
            f.write(raw)
        manifest[p.relative_to(ROOT).as_posix()] = dict(snapshot=filename, sha256=hashlib.sha256(raw).hexdigest())
    with (SNAPSHOT/'manifest.json').open('x', encoding='utf8') as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    for p, raw in planned.items():
        assert p.read_bytes() == before[p], 'concurrent navigation change'
        temp = p.with_name(p.name+'.ctp-bridge.tmp')
        with temp.open('xb') as f:
            f.write(raw)
        os.replace(temp, p)
    new = ('closed_time_path_bridge_review_586.md', 'closed_time_path_bridge_entry.py',
           'closed_time_path_bridge_entry_results.json', 'publish_closed_time_path_bridge.py')
    result = dict(kind='non-numbered method review', date='2026-10-01', latest_scientific_round=586,
                  cumulative_numbered_checks=2989, navigation_files=len(paths), checked_local_links=links,
                  broken_links=0, source_hashes=saved['source_hashes'],
                  artifact_hashes={name: sha(HERE/name) for name in new},
                  navigation_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
                  historical_files_unchanged=True, round587_drafts_preserved=True,
                  goal_unchanged=True, independent_agent_review=False, all_checks_passed=True)
    with TARGET.open('x', encoding='utf8', newline='\n') as f:
        f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('latest_scientific_round', 'navigation_files',
                     'checked_local_links', 'broken_links', 'all_checks_passed')}))


if __name__ == '__main__':
    main()
