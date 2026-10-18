"""Preserve and correct the preacceptance author version after independent review."""
import importlib.util
import json
from hashlib import sha256
from pathlib import Path

base = Path(__file__).resolve().parent
proof = base / 'proof.md'
old_proof = proof.read_bytes()
old_results = (base / 'results.json').read_bytes()
drafts = base / 'drafts'
drafts.mkdir(exist_ok=True)
for name, data in [('proof_before_review_corrections.md', old_proof),
                   ('results_before_review_corrections.json', old_results)]:
    with (drafts / name).open('xb') as handle:
        handle.write(data)
text = old_proof.decode('utf-8')
format_count = text.count('\x0crac')
text = text.replace('\x0crac', r'\frac')
before = '|E|$(N_c,n_g)=(5,3)$ 在局部 R/P/W 结构比较中|$b_w=-7/6$，但精确R禁止单位交换；此例不额外授予所有 $\\mathrm{Spin}^{\\mathbb Z_4}$ 背景的全局资格|'
after = '|E的经典诊断（不是量子R反模型）|$(N_c,n_g)=(5,3)$|$b_w=-7/6$且经典作用可保r；但普通弱instanton已有量子r反常，不能据此证明E对精确量子R独立必要|'
assert before in text
text = text.replace(before, after)
before = '删 E 的例子刻意只针对定理声明的局部及领先阶合同；若另要求无额外扇区的全 $\\mathrm{Spin}^{\\mathbb Z_4}$ 背景，$N_c=5,n_g=3$ 会被必要纯计数排除，不能当那一更强理论的反模型。本定理不依赖这些删弱例的全局升级。'
after = ('删E行只保留经典作用／算符层诊断：$N_c=5,n_g=3$的一个普通$SU(2)$ instanton有'
         '$n_g(N_c+1)=18$个左手零模，测度在$r$下乘$i^{18}=-1$。所以即使不考虑更广的'
         '$\\mathrm{Spin}^{\\mathbb Z_4}$背景，也不能把它当作保精确量子R的删弱反模型。'
         '本轮未证明E对完整量子R独立必要。主必要链只使用R的作用不变后果；'
         '最终三色各代的量子资格另按1060原背景与可加性承接，该修正不改变式(2)。')
assert before in text
text = text.replace(before, after)
assert not any(ord(c) < 32 and c not in '\n\r\t' for c in text)
proof.write_text(text, encoding='utf-8', newline='\n')
spec = importlib.util.spec_from_file_location('round1308_revised', base / 'check.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
new = mod.compute()
old = json.loads(old_results)
keys = [k for k in old if old[k] != new[k]]
assert set(keys) == {'deletion_witnesses', 'unchanged_source_sha256'}, keys
assert new['selected_Nc_ng'] == old['selected_Nc_ng']
payload = (json.dumps(new, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
(base / 'results.json').write_bytes(payload)
record = {'review_changes': ['removed misleading quantum-r claim from without-E diagnostic',
                             'repaired form-feed in displayed frac'],
          'format_repair_count': format_count,
          'main_theorem_and_numeric_values_unchanged': True,
          'changed_result_keys': keys,
          'old_proof_sha256': sha256(old_proof).hexdigest(),
          'new_proof_sha256': sha256(proof.read_bytes()).hexdigest(),
          'old_results_sha256': sha256(old_results).hexdigest(),
          'new_results_sha256': sha256(payload).hexdigest()}
with (base / 'preacceptance_review_corrections.json').open('x', encoding='utf-8') as handle:
    json.dump(record, handle, indent=2)
    handle.write('\n')
print(json.dumps(record))
