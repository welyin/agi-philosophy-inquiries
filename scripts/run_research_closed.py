"""Read-only replay of the closed 1009--1043 stage, without a legacy directory.

Reuses the frozen research_layout runtime.  Only this stage adds a path layer,
runpy dispatch and forwarding of the known Python-script subprocess form.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import types

from research_layout import (Layout, ResearchRuntime, ResearchLoader, ROOT,
                             relocation_original_bytes, no_archive_or_research_writes)

OLD = 'research_cognition_physics/archive_1009_'
NEW = 'research_cognition_physics/archive_1009_1043'
MANIFEST = ROOT/'research_cognition_physics/_migration/closure_1009_1043_20261008/manifest.json'
THIS = Path(__file__).resolve()


class ClosedLayout(Layout):
    def __init__(self):
        super().__init__()
        data = json.loads(MANIFEST.read_text('utf8'))
        assert data['schema'] == 'byte_preserving_stage_closure_v1'
        self.closure = data
        for name, digest in data['unchanged_runtime_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
        for entry in data['entries']:
            logical, actual = entry['original'], ROOT/entry['destination']
            assert logical.startswith(OLD+'/') and entry['destination'].startswith(NEW+'/')
            assert actual.resolve().is_relative_to(ROOT)
            assert logical not in self.relocations
            self.relocations[logical] = (actual, entry)
            self.relocation_reverse[os.path.normcase(os.path.abspath(actual))] = logical
            self.current_files.add(logical)
            self.current_directories.update(p.as_posix() for p in Path(logical).parents)

    def closed_key(self, path):
        try:
            relative = Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
        except ValueError:
            return None
        if relative == NEW or relative.startswith(NEW+'/'):
            return OLD+relative[len(NEW):]
        if relative == OLD or relative.startswith(OLD+'/'):
            return relative
        return None

    def key(self, path, view='legacy'):
        closed = self.closed_key(path)
        if closed:
            return closed
        if view == 'current':
            # This stage was frozen AFTER the 764--1008 split. Its references
            # to other physical stages use their current bytes and names,
            # not the older pre-split historical view.
            try:
                return Path(os.path.abspath(str(path))).relative_to(ROOT).as_posix()
            except ValueError:
                return None
        return super().key(path, view)

    def actual(self, path, view='legacy'):
        key = self.closed_key(path)
        if key:
            return ROOT/(NEW+key[len(OLD):])
        if view == 'current':
            return Path(str(path))
        return super().actual(path, view)

    def read(self, path, historical=True, view='legacy'):
        if view == 'current' and not self.closed_key(path):
            return Path(str(path)).read_bytes()
        return super().read(path, historical=historical, view=view)

    def read_pre_split(self, path):
        key = self.closed_key(path)
        if key and key not in self.relocations:
            return self.actual(path).read_bytes()
        return super().read_pre_split(path)

    def glob_paths(self, directory, pattern, view):
        key = self.closed_key(directory)
        if not key:
            yield from super().glob_paths(directory, pattern, view)
            return
        physical = self.actual(directory)
        if physical.is_dir():
            yield from (ROOT/self.closed_key(p) for p in sorted(physical.glob(pattern)))

    def verify_closure(self):
        assert not (ROOT/OLD).exists(), 'An old shell directory must not be retained'
        changed_markdown = 0
        for entry in self.closure['entries']:
            actual = ROOT/entry['destination']
            raw = actual.read_bytes()
            original = relocation_original_bytes(raw, entry)
            if actual.suffix in ('.py', '.json'):
                assert raw == original
            changed_markdown += raw != original
        return dict(files=len(self.closure['entries']),
                    reversible_markdown_files=changed_markdown,
                    scientific_python_json_bytes_unchanged=True,
                    frozen_runtime_files_unchanged=True,
                    zip_used=False, old_directory_retained=False)


@contextmanager
def execution_routes(store, runtime):
    old_run_path, old_run = runpy.run_path, subprocess.run

    def mapped_run_path(path_name, init_globals=None, run_name=None):
        key = store.script_key(path_name)
        if not key:
            return old_run_path(path_name, init_globals=init_globals, run_name=run_name)
        name = run_name or '<run_path>'
        module = types.ModuleType(name)
        if init_globals:
            module.__dict__.update(init_globals)
        module.__name__ = name
        module.__spec__ = runtime.spec(name, key)
        previous = sys.modules.get(name)
        sys.modules[name] = module
        try:
            ResearchLoader(runtime, key).exec_module(module)
            return module.__dict__.copy()
        finally:
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    def mapped_run(command, *args, **kwargs):
        # This deliberately supports only the actual Python invocation form in
        # the frozen stage. It does not reinterpret arbitrary shell commands.
        if isinstance(command, (list, tuple)) and len(command) >= 5:
            same_python = os.path.normcase(os.path.abspath(str(command[0]))) == os.path.normcase(os.path.abspath(sys.executable))
            if same_python and list(command[1:4]) == ['-B', '-X', 'utf8']:
                key = store.script_key(command[4])
                if key:
                    command = [command[0], '-B', '-X', 'utf8', str(THIS),
                               '--script', key, '--', *command[5:]]
                    if kwargs.get('cwd') is not None:
                        kwargs['cwd'] = str(store.actual(kwargs['cwd']))
        return old_run(command, *args, **kwargs)

    runpy.run_path, subprocess.run = mapped_run_path, mapped_run
    try:
        yield
    finally:
        runpy.run_path, subprocess.run = old_run_path, old_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--script')
    choice.add_argument('--verify-closure', action='store_true')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    sys.dont_write_bytecode = True
    store = ClosedLayout()
    runtime = ResearchRuntime(store)
    sys.path.append(str(ROOT/'.research_runtime'))
    sys.addaudithook(no_archive_or_research_writes)
    with runtime.installed(), execution_routes(store, runtime):
        if args.verify_closure:
            print(json.dumps(store.verify_closure(), ensure_ascii=False))
        else:
            arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
            runtime.run_script(args.script, arguments)


if __name__ == '__main__':
    main()
