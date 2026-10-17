"""Read-only replay of 952 through its hash-verified pre-formula document view.

This is an archive compatibility helper, not a scientific implementation.
It neither changes old evidence nor skips the old script's comparisons.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import Layout, ResearchRuntime, no_archive_or_research_writes
from normalize_markdown_math import inverse_text

TARGETS = (
    'research_cognition_physics/archive_301_341/research_note_303.md',
    'research_cognition_physics/archive_585_628/research_note_588.md',
)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


class TimingReplayLayout(Layout):
    def __init__(self):
        super().__init__()
        manifest = ROOT / 'research_cognition_physics/_migration/markdown_math_20261010/manifest.json'
        entries = {e['path']: e for e in json.loads(manifest.read_text('utf-8'))['entries']}
        result = json.loads((ROOT / 'research_cognition_physics/archive_935_955/952/portal_lapse_source_results.json').read_text('utf-8'))
        expected = {k.replace('\\', '/'): v for k, v in result['source_hashes'].items()}
        self.verified_old = {}
        self.documents = {}
        for rel in TARGETS:
            raw = (ROOT / rel).read_bytes()
            if sha(raw) == expected[rel]:
                old = raw
            else:
                e = entries[rel]
                assert sha(raw) == e['current_sha256'], rel
                old = inverse_text(raw.decode(e['encoding']), e['edits']).encode(e['encoding'])
                assert sha(old) == e['original_sha256'], rel
            assert sha(old) == expected[rel], rel
            self.verified_old[rel] = old
            self.documents[rel] = {'current_sha256': sha(raw), 'historical_sha256': sha(old)}

    def read(self, path, historical=True, view='legacy'):
        try:
            rel = self.actual(path, view).resolve().relative_to(ROOT).as_posix()
        except ValueError:
            rel = None
        if historical and rel in self.verified_old:
            self.read_paths.add(rel)
            return self.verified_old[rel]
        return super().read(path, historical=historical, view=view)


def main():
    assert not sys.argv[1:], 'Read-only helper takes no arguments'
    sys.dont_write_bytecode = True
    store = TimingReplayLayout()
    runtime = ResearchRuntime(store, 'late')
    sys.path.append(str(ROOT / '.research_runtime'))
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed():
        runtime.run_script('archive_935_955/952/portal_lapse_source.py', [])
    for rel, hashes in store.documents.items():
        assert sha((ROOT / rel).read_bytes()) == hashes['current_sha256']
    print(json.dumps({'status': 'PASS_ARCHIVED_952', 'documents': store.documents,
                     'research_files_read_only': True, 'new_scientific_groups': 0}, ensure_ascii=False))


if __name__ == '__main__':
    main()
