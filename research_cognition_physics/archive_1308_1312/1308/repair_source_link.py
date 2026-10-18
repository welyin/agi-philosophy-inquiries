"""One-time preacceptance URL correction; preserves the old result, changes no scientific values."""
import importlib.util
import json
from hashlib import sha256
from pathlib import Path

root = Path(__file__).resolve().parent
old_bytes = (root / 'results.json').read_bytes()
old = json.loads(old_bytes)
spec = importlib.util.spec_from_file_location('round1308_author', root / 'check.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
new = module.compute()
key = 'unchanged_source_sha256'
assert {k: v for k, v in old.items() if k != key} == {k: v for k, v in new.items() if k != key}
changed = [p for p in old[key] if old[key][p] != new[key][p]]
assert changed == ['research_cognition_physics/archive_1308_/1308/proof.md'], changed
draft = root / 'drafts'
draft.mkdir(exist_ok=True)
with (draft / 'results_before_source_link_fix.json').open('xb') as handle:
    handle.write(old_bytes)
new_bytes = (json.dumps(new, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
(root / 'results.json').write_bytes(new_bytes)
receipt = {'reason': 'Roy 1907.10238 source URL changed from unavailable HTML to verified PDF',
           'scientific_values_unchanged': True, 'changed_binding': changed,
           'old_results_sha256': sha256(old_bytes).hexdigest(),
           'new_results_sha256': sha256(new_bytes).hexdigest()}
with (root / 'source_link_fix.json').open('x', encoding='utf-8') as handle:
    json.dump(receipt, handle, indent=2)
    handle.write('\n')
print(json.dumps(receipt))
