"""Regression checks for relocation, evidence integrity, and read-only replay."""
import importlib
import unittest
from pathlib import Path

from research_layout import Layout, ResearchRuntime, ROOT, original_bytes, sha, no_archive_or_research_writes


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
        original = original_bytes(current, entry)
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
            original_bytes(actual.read_bytes(), dict(entry, original_sha256='0'*64))

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


if __name__ == '__main__':
    unittest.main()
