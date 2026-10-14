"""Read-only historical replay beneath the 1219 thermal-geometry navigation."""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import run_research_thermodynamics_closed as previous
from research_layout import ROOT, ResearchLoader, no_archive_or_research_writes, relocation_original_bytes

THIS = Path(__file__).resolve()
STAGE = 'research_cognition_physics/archive_1219_'
LAYER = ROOT / 'research_cognition_physics/_migration/thermal_geometry_active_1219/layer.json'
NAVIGATION = {'research_cognition_physics/' + name for name in
              ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}


class ThermalGeometryLayout(previous.ThermalLayout):
    def __init__(self):
        raw = LAYER.read_bytes()
        self.geometry_layer, self.geometry_sha256 = json.loads(raw), previous.digest(raw)
        if self.geometry_layer.get('schema') != 'thermal_geometry_active_navigation_v1':
            raise ValueError('Unexpected thermal-geometry navigation layer')
        self.geometry_entries = {}
        for entry in self.geometry_layer['entries']:
            old, new = previous.previous.validate_entry(entry)
            if old not in NAVIGATION or new != old or old in self.geometry_entries:
                raise ValueError('Only four distinct root navigation files may change')
            self.geometry_entries[old] = entry
        if set(self.geometry_entries) != NAVIGATION:
            raise ValueError('Exactly four root navigation files are required')
        super().__init__()
        for name, entry in self.geometry_entries.items():
            boundary = self.thermal_entries[name]
            if (entry['original_sha256'], entry['original_bytes']) != (boundary['current_sha256'], boundary['current_bytes']):
                raise ValueError('Active navigation does not join closed thermal boundary: ' + name)

    def before_thermal(self, path):
        name = self.thermal_name(path)
        if previous.moved(name) and name not in self.thermal_entries:
            raise FileNotFoundError('File absent from thermal closure inventory: ' + name)
        raw = self.thermal_actual(path).read_bytes()
        active = self.geometry_entries.get(name)
        if active:
            raw = relocation_original_bytes(raw, active)
        entry = self.thermal_entries.get(name)
        prior = relocation_original_bytes(raw, entry) if entry else raw
        if entry and hasattr(self, 'read_paths'):
            self.read_paths.add(entry['destination'])
            if raw != prior:
                self.restored_markdown.add(entry['destination'])
        return prior

    def script_key(self, path):
        name = previous.relative(path)
        if previous.inside(name, STAGE) and name.endswith('.py') and Path(str(path)).is_file():
            return name
        return super().script_key(path)

    def script_actual(self, key):
        return ROOT / key if previous.inside(key, STAGE) else super().script_actual(key)

    def verify_active_navigation(self):
        if previous.digest(LAYER.read_bytes()) != self.geometry_sha256:
            raise ValueError('Active layer changed during replay')
        for name, entry in self.geometry_entries.items():
            relocation_original_bytes((ROOT / name).read_bytes(), entry)
        for name, expected in self.geometry_layer['unchanged_support_sha256'].items():
            if previous.digest((ROOT / previous.safe_relative(name)).read_bytes()) != expected:
                raise ValueError('Frozen support changed: ' + name)
        return dict(status='PASS',files=4,manifest_sha256=self.geometry_sha256,
                    closed_boundary_verified=True,original_navigation_recoverable=True,
                    frozen_archives_modified=False)

    def verify_layout(self):
        return dict(active_navigation=self.verify_active_navigation(),
                    frozen_layers=super().verify_layout(),all_checks_passed=True)


class ThermalGeometryLoader(ResearchLoader):
    def exec_module(self, module):
        if not previous.inside(self.key, STAGE):
            return previous.ThermalLoader(self.runtime, self.key).exec_module(module)
        actual = self.runtime.store.script_actual(self.key)
        module.__file__ = str(actual)
        old_path = list(sys.path)
        sys.path.insert(0, str(actual.parent))
        try:
            exec(compile(actual.read_bytes(), str(actual), 'exec'), module.__dict__)
        finally:
            sys.path[:] = old_path


class ThermalGeometryRuntime(previous.ThermalRuntime):
    def spec(self, name, key):
        return importlib.util.spec_from_loader(name, ThermalGeometryLoader(self, key),
                                              origin=str(self.store.script_actual(key)))


@contextmanager
def closed_entry_route():
    old_run = subprocess.run
    def mapped(command, *args, **kwargs):
        if isinstance(command, (list, tuple)) and len(command) >= 6:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(sys.executable)
            if (same_python and list(command[1:4]) == ['-B', '-X', 'utf8'] and
                    previous.relative(command[4]) == 'scripts/run_research_thermodynamics_closed.py' and
                    command[5] in ('--verify-closure','--verify-navigation','--verify-layout','--script')):
                command = [*command[:4], str(THIS), *command[5:]]
        return old_run(command, *args, **kwargs)
    subprocess.run = mapped
    try:
        yield
    finally:
        subprocess.run = old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-layout', action='store_true')
    choice.add_argument('--verify-closure', action='store_true')
    choice.add_argument('--verify-navigation', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ThermalGeometryLayout()
    runtime = ThermalGeometryRuntime(store)
    sys.path.append(str(ROOT / '.research_runtime'))
    routes = previous.previous.combined.merged.previous.previous
    routes.THIS = THIS
    previous.THIS = THIS
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), routes.execution_routes(previous.ExecutionLayout(store), runtime), previous.gr_verification_route(), closed_entry_route():
        if args.script:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)
        else:
            result = (store.verify_layout() if args.verify_layout else
                      store.verify_closure() if args.verify_closure else store.verify_active_navigation())
            print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
