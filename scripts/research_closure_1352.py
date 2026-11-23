"""Read-only replay of the single physical archive for rounds 1323--1352."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import builtins
import importlib.machinery
import importlib.util
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from urllib.parse import unquote

from research_layout import ROOT, no_archive_or_research_writes
from research_navigation_1353 import before_stage

BASE = ROOT / 'research_cognition_physics'
OLD = BASE / 'archive_1323_'
ARCHIVE = BASE / 'archive_1323_1352'
MIGRATION = BASE / '_migration/closure_1323_1352'
MANIFEST = MIGRATION / 'manifest.json'
PRODUCTS = MIGRATION / 'products.json'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def relative(path):
    if not isinstance(path, (str, bytes, os.PathLike)):
        return None
    return Path(os.path.abspath(os.fsdecode(path)))


def physical(path):
    if path is not None and path.is_relative_to(OLD):
        return ARCHIVE / path.relative_to(OLD)
    return path


def bounded(base, name):
    rel = PurePosixPath(name)
    require(not rel.is_absolute() and '..' not in rel.parts and '\\' not in name
            and ':' not in name and rel.as_posix() == name, 'Invalid path: ' + name)
    path = base / name
    require(path.resolve().is_relative_to(base), 'Path leaves archive: ' + name)
    return path


def restore(raw, entry):
    require(digest(raw) == entry['current_sha256'], 'Closure document changed: ' + entry['logical'])
    # Reverse the mechanical directory-link spelling first, then the text edits.
    text = raw.decode('utf8').replace('archive_1323_1352/', 'archive_1323_/')
    pieces, at = [], 0
    for start, end, original in entry['reverse_edits']:
        require(at <= start <= end <= len(text), 'Invalid edit: ' + entry['logical'])
        pieces.extend((text[at:start], original))
        at = end
    raw = (''.join(pieces) + text[at:]).encode('utf8')
    require(digest(raw) == entry['original_sha256'], 'Failed historical digest: ' + entry['logical'])
    return raw


class StageClosure:
    def __init__(self):
        self.data = json.loads(MANIFEST.read_bytes())
        require(self.data['schema'] == 'closed_research_1323_1352_v1', 'Unexpected closure schema')
        require(not OLD.exists() and ARCHIVE.is_dir(), 'Expected a single closed archive')
        self.cache = {}
        for entry in self.data['changed_markdown']:
            path = physical(bounded(ROOT, entry['logical']))
            require(path.suffix == '.md' and path not in self.cache, 'Invalid duplicate Markdown edit')
            self.cache[path] = restore(before_stage(path), entry)

    def before(self, path):
        actual = physical(path)
        return self.cache[actual] if actual in self.cache else actual.read_bytes()

    @contextmanager
    def historical_paths(self):
        # The old verification code sees its own names and original Markdown.
        # Relative imports of the original repository's modules stay untouched.
        path_open, native_open, path_stat = Path.open, builtins.open, Path.stat
        def stream(file, mode='r', buffering=-1, encoding=None, errors=None, newline=None, closefd=True, opener=None):
            path = relative(file)
            actual = physical(path)
            if actual != path or actual in self.cache:
                if any(c in mode for c in 'wax+'):
                    raise PermissionError('Historical archive is read-only: ' + str(file))
                if actual in self.cache:
                    require(mode in ('r', 'rt', 'rb') and opener is None and closefd, 'Invalid historical read mode')
                    raw = io.BytesIO(self.cache[actual])
                    return raw if 'b' in mode else io.TextIOWrapper(raw, encoding=encoding or 'utf8', errors=errors, newline=newline)
                file = actual
            return native_open(file, mode, buffering, encoding, errors, newline, closefd, opener)

        def opened(path, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
            candidate = path if path.is_absolute() else relative(path)
            if candidate.is_relative_to(OLD) or physical(candidate) in self.cache:
                return stream(path, mode, buffering, encoding, errors, newline)
            return path_open(path, mode, buffering, encoding, errors, newline)

        def stat(path, **kwargs):
            candidate = path if path.is_absolute() else relative(path)
            return path_stat(physical(candidate) if candidate.is_relative_to(OLD) else path, **kwargs)

        spec_from_file = importlib.util.spec_from_file_location
        class HistoricalSource(importlib.machinery.SourceFileLoader):
            def get_code(self, fullname):
                return compile(Path(self.path).read_bytes(), self.path, 'exec')

        def source_spec(name, location=None, *args, **kwargs):
            path = relative(location)
            if path is not None and path.is_relative_to(OLD) and path.suffix == '.py' and not kwargs.get('loader'):
                kwargs['loader'] = HistoricalSource(name, str(path))
            return spec_from_file(name, location, *args, **kwargs)

        try:
            importlib.util.spec_from_file_location = source_spec
            Path.open, Path.stat, builtins.open = opened, stat, stream
            yield
        finally:
            importlib.util.spec_from_file_location = spec_from_file
            Path.open, Path.stat, builtins.open = path_open, path_stat, native_open

    def verify(self):
        inventory = self.data['prior_archive_sha256']
        require(all(digest(self.before(OLD / name)) == expected for name, expected in inventory.items()),
                'Inventory of original archive files changed')
        for name, expected in self.data['unchanged_manifests_sha256'].items():
            require(digest(bounded(ROOT, name).read_bytes()) == expected, 'Frozen manifest changed: ' + name)
        accepted, checks = 0, []
        for n in range(1323, 1353):
            acceptance = json.loads((ARCHIVE / str(n) / 'acceptance.json').read_bytes())
            require(acceptance['round'] == n and acceptance['status'].startswith('ACCEPTED'), 'Invalid acceptance')
            for name, expected in acceptance['sha256'].items():
                logical = bounded(ROOT, name)
                require(logical.is_relative_to(OLD) and digest(self.before(logical)) == expected, 'Evidence changed: ' + name)
                accepted += 1
            for record in acceptance['executed_checks']:
                spec = record['script'].split()
                path = bounded(ROOT, spec[0]) if spec[0].startswith('research_cognition_physics/') else bounded(OLD / str(n), spec[0])
                require(path.is_relative_to(OLD / str(n)) and path.suffix == '.py' and record['exit_code'] == 0,
                        'Invalid accepted scientific run: ' + str(path))
                extra = spec[1:]
                compare = n == 1351 and extra == ['certificate()', 'vs', 'results.json']
                require(extra in ([], ['--check']) or compare, 'Unexpected scientific replay arguments')
                checks.append((str(physical(path)), [] if compare else extra, compare))
        products = json.loads(PRODUCTS.read_bytes())
        require(products['closure_sha256'] == digest(MANIFEST.read_bytes()), 'Products and closure differ')
        for name, expected in products['files_sha256'].items():
            require(digest(before_stage(bounded(ROOT, name))) == expected, 'Closure product changed: ' + name)
        files = set(self.cache) | {bounded(ROOT,n) for n in products['files_sha256'] if n.endswith('.md')}
        links = 0
        for path in files:
            for match in re.finditer(r'!?\[[^\]\n]*\]\((<[^>]+>|[^)\n]+)\)', path.read_text(encoding='utf-8-sig')):
                name = match[1].strip()
                name = name[1:-1] if name.startswith('<') else re.sub(r'\s+"[^"\n]*"$', '', name)
                name = unquote(re.sub(r':\d+$', '', name.partition('#')[0]))
                if not name or re.match(r'^(https?|mailto|app|codex|plugin):', name, re.I):
                    continue
                if re.match(r'^/[A-Za-z]:[/\\]', name):
                    name = name[1:]
                require((path.parent / name).exists(), f'Broken link: {path} -> {name}')
                links += 1
        return dict(status='PASS', original_files=len(inventory), modified_markdown_with_verified_preimages=len(self.cache),
                    accepted_rounds=30, accepted_evidence_bindings=accepted, recorded_scientific_runs=len(checks),
                    frozen_manifests=len(self.data['unchanged_manifests_sha256']), checked_links=links,
                    products=len(products['files_sha256']), closure_sha256=products['closure_sha256']), checks


def execute_science(path, arguments, compare):
    script = Path(path).resolve()
    sys.dont_write_bytecode = True
    sys.addaudithook(no_archive_or_research_writes)
    logical = OLD / script.relative_to(ARCHIVE)
    # Some early receipts named just the script, even where the frozen CLI
    # requires its explicit read-only switch or publishes on a bare call.
    requires_check = {
        '1333/check.py', '1333/independent_check.py', '1334/check.py', '1334/independent_check.py',
        '1335/check.py', '1335/optimal_bound_check.py', '1336/check.py',
        '1336/window/independent_check.py', '1337/swap_independent_check.py',
    }
    if not arguments and script.relative_to(ARCHIVE).as_posix() in requires_check:
        arguments = ['--check']
    sys.argv = [str(logical), *arguments]
    sys.path.insert(0, str(script.parent))
    # Early certificates store logical path strings as well as per-file digests.
    # Keep the physical import path but execute source under its recorded name.
    namespace = dict(__name__='certificate_replay' if compare else '__main__', __file__=str(logical))
    closure = StageClosure()
    with closure.historical_paths():
        exec(compile(logical.read_bytes(), str(logical), 'exec'), namespace)
        if compare:
            computed = namespace['certificate']()
            require(json.loads(json.dumps(computed, allow_nan=False)) == json.loads(script.with_name('results.json').read_bytes()), 'Saved certificate changed')


def replay_science(spec):
    path, arguments, compare = spec
    code = ('import sys;sys.path.insert(0,"scripts");from research_closure_1352 import execute_science;'
            'execute_science(sys.argv[1],sys.argv[3:],sys.argv[2]=="compare")')
    proc = subprocess.run([sys.executable,'-B','-X','utf8','-c',code,path,'compare' if compare else 'run',*arguments],
                          cwd=ROOT, capture_output=True, text=True, encoding='utf8')
    require(proc.returncode == 0, f'{path}: {proc.stderr or proc.stdout}')
    return Path(path).relative_to(ARCHIVE).as_posix()


def run(main):
    closure = StageClosure()
    flags = sys.argv[1:]
    if flags in (['--verify-stage-closure'], ['--replay-stage-science']):
        result, checks = closure.verify()
        if flags[0] == '--replay-stage-science':
            with ThreadPoolExecutor(max_workers=2) as pool:
                result['scientific_replays_passed'] = list(pool.map(replay_science, checks))
        print(json.dumps(result, ensure_ascii=False))
        return
    if flags == ['--verify-current-progress']:
        sys.argv[1:] = ['--verify-reclassification']
    if flags[:1] == ['--script'] and len(flags) > 1:
        path = relative(flags[1])
        if path.is_relative_to(ARCHIVE):
            sys.argv[2] = str(OLD / path.relative_to(ARCHIVE))
        if physical(path).is_relative_to(ARCHIVE):
            sys.path.insert(0, str(physical(path).parent))
    with closure.historical_paths():
        main()

