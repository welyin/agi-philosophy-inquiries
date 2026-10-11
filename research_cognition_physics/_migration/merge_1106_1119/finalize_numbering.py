"""Apply the user's numbering exception to this still-open maintenance operation.

No old scientific or historical migration bytes change. Recompute the new
outer layer against its exact original preimages, not against an intermediate
numbered layout. Refuse reruns and concurrent edits.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

HERE = Path(__file__).resolve().parent
BASE = ROOT / 'research_cognition_physics'
NAME = '条件生成链与物理选择的剩余自由_阶段论文.md'
OLD_NAME = '06_' + NAME


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    manifest_path = HERE / 'manifest.json'
    manifest_before = manifest_path.read_bytes()
    manifest = json.loads(manifest_before)
    assert len(manifest['paper_order']) == 6, 'Numbering exception already applied'
    assert (BASE / OLD_NAME).is_file() and not (BASE / NAME).exists()
    planned = []
    for entry in manifest['entries']:
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        assert sha(raw) == entry['current_sha256']
        if path.suffix != '.md':
            continue
        original = relocation_original_bytes(raw, entry)
        encoding = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf-8'
        before = raw.decode(encoding)
        after = before.replace(OLD_NAME, NAME)
        if path.parent == BASE and path.name in {
                'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'}:
            row = '|06|[条件生成链与剩余自由（1009—1043）](' + NAME + ')|'
            assert row in after
            after = after.replace(row + ('\r\n' if '\r\n' in after else '\n'), '')
            old_sentence = '第06篇是场论与引力方向的条件综合，不表示前五篇已满足其全部前提。'
            assert old_sentence in after
            after = after.replace(old_sentence,
                '另保留[条件生成链与物理选择的剩余自由](' + NAME + ')，'
                '按用户要求暂不参与阶段论文编号，文件位置和内容保留。')
        destination = BASE / NAME if path == BASE / OLD_NAME else path
        if raw == after.encode(encoding) and path == destination:
            continue
        assert destination.resolve().is_relative_to(ROOT)
        if destination != path:
            assert not destination.exists()
        new = after.encode(encoding)
        original_encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        entry.update(destination=destination.relative_to(ROOT).as_posix(),
                     current_sha256=sha(new), current_bytes=len(new),
                     edits=text_edits(original.decode(original_encoding), after))
        assert relocation_original_bytes(new, entry) == original
        planned.append((path, destination, raw, new))
    manifest['paper_order'] = [p for p in manifest['paper_order'] if p['number'] <= 5]
    manifest['unnumbered_papers'] = [dict(name=NAME,
        reason='User requested no sequence number; filename suffix, location and content retained.')]
    manifest['numbering_amendment'] = dict(date='2026-10-11',
        preceding_manifest_sha256=sha(manifest_before),
        instruction='条件生成链文档不用给序号；其余五篇编号不变。',
        changed_markdown=sum(raw != new for _, _, raw, new in planned),
        renamed_files=sum(p != q for p, q, _, _ in planned))
    destinations = [e['destination'] for e in manifest['entries']]
    assert len(destinations) == len(set(destinations))
    assert manifest_path.read_bytes() == manifest_before
    for path, destination, raw, new in planned:
        assert path.read_bytes() == raw
        if path != destination:
            path.rename(destination)
        if raw != new:
            destination.write_bytes(new)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', amendment=manifest['numbering_amendment'],
                         numbered_papers=5, unnumbered=NAME), ensure_ascii=False))


if __name__ == '__main__':
    main()
