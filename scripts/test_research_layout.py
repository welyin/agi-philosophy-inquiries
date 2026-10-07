"""Regression checks for relocation, evidence integrity, and read-only replay."""
import importlib
import importlib.util
import builtins
import json
import tempfile
import unittest
from pathlib import Path

from research_layout import (Layout, ResearchRuntime, ROOT, original_bytes,
    relocation_original_bytes, sha, no_archive_or_research_writes)


class CurrentLayoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.layout = Layout()

    def test_both_script_names_select_the_same_physical_file(self):
        old = self.layout.lookup_script('joint_wick_source_anomaly.py')
        new = self.layout.lookup_script('archive_764_/775/joint_wick_source_anomaly.py')
        self.assertEqual(old, new)
        self.assertTrue(self.layout.entries[old][0].is_file())

    def test_scientific_source_is_read_from_its_current_directory(self):
        key = self.layout.lookup_script('joint_wick_source_anomaly.py')
        actual, entry = self.layout.entries[key]
        self.assertEqual(self.layout.read(ROOT / key), actual.read_bytes())
        self.assertEqual(sha(actual.read_bytes()), entry['original_sha256'])

    def test_markdown_history_is_verified_from_current_bytes(self):
        key = 'research_cognition_physics/archive_231_/research_note_775.md'
        actual, entry = self.layout.entries[key]
        current = actual.read_bytes()
        original = original_bytes(self.layout.read_pre_split(actual), entry)
        self.assertNotEqual(original, current)
        self.assertEqual(sha(original), entry['original_sha256'])
        self.assertEqual(self.layout.read(ROOT / key, historical=False), current)
        self.assertEqual(actual.read_bytes(), current)

    def test_current_file_tampering_is_rejected(self):
        key = self.layout.lookup_script('joint_wick_source_anomaly.py')
        actual, entry = self.layout.entries[key]
        with self.assertRaises(ValueError):
            original_bytes(actual.read_bytes() + b'changed', entry)

    def test_incorrect_original_hash_is_not_accepted(self):
        key = 'research_cognition_physics/archive_231_/research_note_775.md'
        actual, entry = self.layout.entries[key]
        with self.assertRaises(ValueError):
            original_bytes(self.layout.read_pre_split(actual), dict(entry, original_sha256='0'*64))

    def test_cross_round_read_resolves_without_old_directory(self):
        runtime = ResearchRuntime(self.layout)
        key = self.layout.lookup_script('joint_wick_source_anomaly.py')
        here = runtime.path(self.layout.entries[key][0]).resolve().parent
        note = here / 'research_note_735.md'
        self.assertTrue(note.exists())
        self.assertEqual(sha(note.read_bytes()), self.layout.entries[self.layout.key(note)][1]['original_sha256'])
        self.assertFalse((ROOT / 'research_cognition_physics/archive_231_').exists())

    def test_archive_and_navigation_backup_reads_are_rejected(self):
        for path in ('fake.zip', 'navigation_before_round773_20261004/0_README.md'):
            with self.assertRaises(PermissionError):
                no_archive_or_research_writes('open', (str(ROOT / path), 'r', 0))

    def test_replay_cannot_write_evidence(self):
        runtime = ResearchRuntime(self.layout)
        key = self.layout.lookup_script('joint_wick_source_anomaly.py')
        with self.assertRaises(PermissionError):
            runtime.path(ROOT / key).write_text('must not be written')
        with self.assertRaises(PermissionError):
            no_archive_or_research_writes('os.remove', (str(ROOT / key), -1))


class SplitLayoutTests(unittest.TestCase):
    """Small real files exercise both history layers even before the real move."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'project'
        self.base = self.root / 'research_cognition_physics'
        self.migration = self.base / '_migration'
        self.migration.mkdir(parents=True)
        self.prior = 'research_cognition_physics/archive_764_/'
        self.entries = []
        self.original = b'\xef\xbb\xbf# 775\r\n[x](old.md)\r\n'
        self.intermediate = b'\xef\xbb\xbf# 775\r\n[x](../archive_700/old.md)\r\n'
        self.latest = b'\xef\xbb\xbf# 775\r\n[x](../archive_700_750/old.md)\r\n'
        old_text = self.original.decode('utf-8-sig')
        self.add_file('research_note_775.md',
            'research_cognition_physics/archive_764_776/research_note_775.md',
            self.intermediate, self.latest, '../archive_700/old.md', '../archive_700_750/old.md')
        legacy_entry = dict(original='research_note_775.md',
            destination='archive_764_/research_note_775.md',
            original_sha256=sha(self.original), current_sha256=sha(self.intermediate),
            rewrite_markdown_links=True,
            link_edits=[dict(start=old_text.index('old.md'), end=old_text.index('old.md') + 6,
                            before='old.md', after='../archive_700/old.md')])
        (self.migration / 'manifest.json').write_text(json.dumps(dict(entries=[legacy_entry])), 'utf8')
        early = self.migration / 'layout_001_230_20261004'
        early.mkdir()
        (early / 'manifest.json').write_text(json.dumps(dict(entries=[])), 'utf8')
        for n in (776, 777):
            self.add_file(f'research_note_{n}.md',
                f'research_cognition_physics/archive_{"764_776" if n == 776 else "777_800"}/research_note_{n}.md',
                f'# {n}\n'.encode())
        self.add_file('776/data.json', 'research_cognition_physics/archive_764_776/776/data.json', b'{"n": 7}\n')
        self.add_file('777/peer.py', 'research_cognition_physics/archive_777_800/777/peer.py',
            b'from pathlib import Path\nHERE = Path(__file__).resolve().parent\ndef value(): return HERE.parent.name\n')
        driver = '''from pathlib import Path
import json, importlib.util, hashlib
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
def run():
    with open(HERE/'data.json', encoding='utf8') as stream: data=json.load(stream)
    spec=importlib.util.spec_from_file_location('_split_test_peer',STAGE/'777/peer.py')
    peer=importlib.util.module_from_spec(spec);spec.loader.exec_module(peer)
    return dict(n=data['n'], peer=peer.value(), stage=STAGE.name,
        here=str(HERE.relative_to(ROOT)).replace(chr(92),'/'),
        note_hash=hashlib.sha256((STAGE/'research_note_775.md').read_bytes()).hexdigest(),
        glob=sorted(p.name for p in STAGE.glob('research_note_*.md')),
        recursive=sorted(p.name for p in STAGE.parent.rglob('research_note_*.md')),
        direct_recursive=sorted(p.name for p in HERE.rglob('*.json')))
def attempt_write():
    with open(HERE/'data.json','w') as stream: stream.write('changed')
'''
        self.add_file('776/driver.py', 'research_cognition_physics/archive_764_776/776/driver.py', driver.encode())
        self.navigation_before = '# 进展\r\n第1008轮\r\n'.encode('utf-8-sig')
        self.navigation_after = '# 进展\r\n阶段整理\r\n'.encode('utf-8-sig')
        self.add_file('navigation.md', 'research_cognition_physics/README.md',
            self.navigation_before, self.navigation_after, '第1008轮', '阶段整理')
        self.entries[-1]['original'] = 'research_cognition_physics/README.md'
        self.entries[-1]['edits'][0]['kind'] = 'navigation_or_reference'
        loose = self.base / 'archive_001_003/research_note_1.md'
        loose.parent.mkdir()
        loose.write_text('# 1\n', 'utf8')
        folder = self.migration / 'layout_764_1008_20261008'
        folder.mkdir()
        (folder / 'manifest.json').write_text(json.dumps(dict(version=1, entries=self.entries)), 'utf8')
        self.layout = Layout(self.root)

    def add_file(self, name, destination, original, current=None, before=None, after=None):
        current = original if current is None else current
        path = self.root / destination
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(current)
        edits = []
        if before is not None:
            start = original.decode('utf-8-sig').index(before)
            edits.append(dict(start=start, end=start + len(before), before=before, after=after, kind='link'))
        self.entries.append(dict(original=self.prior + name, destination=destination,
            original_sha256=sha(original), current_sha256=sha(current),
            original_bytes=len(original), current_bytes=len(current), edits=edits))

    def test_two_markdown_eras_are_distinct_and_bom_crlf_are_preserved(self):
        old = self.root / 'research_cognition_physics/archive_231_/research_note_775.md'
        prior = self.root / (self.prior + 'research_note_775.md')
        actual = self.root / self.entries[0]['destination']
        self.assertEqual(self.layout.read(old), self.original)
        self.assertEqual(self.layout.read(old, view='current'), self.original)
        self.assertEqual(self.layout.read(prior, view='current'), self.intermediate)
        self.assertEqual(self.layout.read(actual, view='current'), self.intermediate)
        self.assertEqual(self.layout.read(old, historical=False), self.latest)
        self.assertEqual(actual.read_bytes(), self.latest)
        self.assertFalse((self.base / 'archive_764_').exists())

    def test_legacy_counts_do_not_include_new_relocation_entries(self):
        self.assertEqual(self.layout.verify(), dict(mapped_files=1, checked_files=1,
            unchanged_files=0, verified_link_only_changes=1,
            optional_local_files_skipped=0, zip_files_read=0))
        current = self.layout.verify_relocation()
        self.assertEqual(current['checked_files'], 7)
        self.assertEqual(current['rewritten_markdown_files'], 2)
        self.assertEqual(current['recorded_text_edits'], 2)

    def test_root_navigation_semantic_edit_is_a_separate_reversible_layer(self):
        navigation = self.base / 'README.md'
        self.assertEqual(self.layout.read(navigation, view='current'), self.navigation_before)
        self.assertEqual(self.layout.read(navigation, historical=False, view='current'), self.navigation_after)
        current_path = self.layout.path_class('current')
        self.assertEqual(str(current_path(self.base)), str(self.base))
        self.assertEqual(str(current_path(self.root)), str(self.root))

    def test_current_script_names_importlib_open_and_globs_replay(self):
        old_key = self.layout.lookup_script(self.prior + '776/driver.py')
        new_key = self.layout.lookup_script('archive_764_776/776/driver.py')
        self.assertEqual(old_key, new_key)
        runtime = ResearchRuntime(self.layout)
        before_open = builtins.open
        before_spec = importlib.util.spec_from_file_location
        with runtime.installed():
            spec = runtime.spec('_split_test_driver', new_key)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            result = module.run()
            with self.assertRaises(PermissionError):
                module.attempt_write()
        self.assertIs(builtins.open, before_open)
        self.assertIs(importlib.util.spec_from_file_location, before_spec)
        self.assertEqual(result, dict(n=7, peer='archive_764_', stage='archive_764_',
            here='research_cognition_physics/archive_764_/776',
            note_hash=sha(self.intermediate),
            glob=['research_note_775.md', 'research_note_776.md', 'research_note_777.md'],
            recursive=['research_note_1.md', 'research_note_775.md', 'research_note_776.md', 'research_note_777.md'],
            direct_recursive=['data.json']))

    def test_relocation_rejects_tampered_files_and_bad_edit_records(self):
        entry = self.entries[0]
        with self.assertRaises(ValueError):
            relocation_original_bytes(self.latest + b'changed', entry)
        bad_edit = dict(entry['edits'][0], before='wrong')
        with self.assertRaises(ValueError):
            relocation_original_bytes(self.latest, dict(entry, edits=[bad_edit]))
        with self.assertRaises(ValueError):
            relocation_original_bytes(self.latest, dict(entry, original_sha256='0' * 64))

    def test_non_markdown_edits_are_not_authorized_by_the_manifest(self):
        entry = dict(self.entries[0], destination='science.py')
        with self.assertRaises(ValueError):
            relocation_original_bytes(self.latest, entry)


if __name__ == '__main__':
    unittest.main()
