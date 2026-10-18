"""Register the user-authorized conditional goal without changing historical bytes."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from research_layout import relocation_original_bytes

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1290_'
LAYER = BASE / '_migration/structure_selection_active_1288/layer.json'
RECEIPT = STAGE / '_shared/goal_restart_navigation_after1306.json'

def main():
    if RECEIPT.exists():
        raise SystemExit('Restart already recorded; do not repeat.')
    old_layer_raw = LAYER.read_bytes()
    layer = json.loads(old_layer_raw)
    for name, expected in layer['unchanged_support_sha256'].items():
        assert sha((ROOT / name).read_bytes()) == expected, name
    goal = STAGE / '_shared/goal_after1306.md'
    historical_goal = STAGE / '_shared/goal.md'
    assert sha(historical_goal.read_bytes()) == '4bd9fa800ac8385b0850e13526f4f5356b62a287464bb13b1cacc1a22370072e'
    prefix = '''
## 当前目标：标准模型的条件重建与共同相容性

2026-10-18。用户删除旧目标后要求重启，并明确“本次只需达成条件性成立”。[新正式目标](archive_1290_/_shared/goal_after1306.md)已经启用。结项要求为前件／依赖表、条件结构重建链、同一父对象与公开有效范围证书、论文和独立验收；不再要求唯一选择群或三代、排尽替代或消除全部额外输入。

研究接续于1290—1306共17轮，历史结果和[原目标](archive_1290_/_shared/goal.md)保留。[1307准入](archive_1290_/_shared/round1307_admission.md)开始检验新增中性N的Dirac质量耦合能否迫使其与带荷物质采用同一领先应力源权，并核对完整父作用到1305背景过程的真实映射。1307尚未验收，新目标未结项，本次目标登记计0。

采用标准群、表示、三代和具体参数时明确作为条件；已推出的系数关系与这些输入分开。SM加右手中微子的Dirac分支与严格最小SM分别标明；背景层过程不冒作全部相互作用量子预测。下方旧“当前长程目标”和状态均属历史，以本节为准。
'''
    output_hashes = {}
    for entry in layer['entries']:
        path = ROOT / entry['destination']
        raw = path.read_bytes()
        original = relocation_original_bytes(raw, entry)
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        assert len(entry['edits']) == 1
        edit = entry['edits'][0]
        assert edit['start'] == edit['end'] and edit['before'] == ''
        marker = '\n## 当前长程目标：从认知组合到标准模型结构的共同生成'
        assert marker in edit['after']
        historical = edit['after'][edit['after'].index(marker):]
        edit['after'] = prefix + historical
        old_text = original.decode(encoding)
        updated = (old_text[:edit['start']] + edit['after'] + old_text[edit['end']:]).encode(encoding)
        entry['current_sha256'] = sha(updated)
        entry['current_bytes'] = len(updated)
        assert relocation_original_bytes(updated, entry) == original
        path.write_bytes(updated)
        output_hashes[entry['destination']] = sha(updated)
    stage_path = STAGE / 'README.md'
    raw = stage_path.read_bytes()
    text = raw.decode('utf-8')
    marker = '\n- [正式目标、五条主线与整体结项标准]'
    assert marker in text
    history = text[text.index(marker):]
    stage_text = '# 1290起：标准模型结构研究\n' + prefix.replace('archive_1290_/', '') + '\n## 原长程计划与历史接续\n' + history
    stage_path.write_bytes(stage_text.encode('utf-8'))
    output_hashes[str(stage_path.relative_to(ROOT)).replace('\\', '/')] = sha(stage_path.read_bytes())
    layer['latest_goal_path'] = str(goal.relative_to(ROOT)).replace('\\', '/')
    layer['active_goal_path'] = layer['latest_goal_path']
    layer['application_goal_status_last_observed'] = 'active'
    layer['application_goal_status_note'] = 'User deleted previous goal; create_goal succeeded for conditional SM reconstruction on 2026-10-18.'
    layer['in_progress_round'] = 1307
    layer['current_state'] = 'Conditional SM reconstruction goal active. 1290-1306 remain accepted (17 rounds); 1307 admitted for neutral Dirac source classification and parent/background mapping. No new acceptance or axioms yet.'
    new_layer = (json.dumps(layer, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
    LAYER.write_bytes(new_layer)
    receipt = dict(schema='conditional_sm_goal_restart_navigation_v1', date='2026-10-18', application_goal_status='active',
                   latest_accepted_round=1306, in_progress_round=1307, scientific_round_increment=0, new_axioms=0,
                   historical_goal_sha256=sha(historical_goal.read_bytes()), new_goal_sha256=sha(goal.read_bytes()),
                   old_layer_sha256=sha(old_layer_raw), new_layer_sha256=sha(new_layer),
                   originals_recoverable=len(layer['entries']), support_bindings_verified=len(layer['unchanged_support_sha256']),
                   files_sha256=output_hashes)
    with RECEIPT.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps({'status':'PASS', 'goal':'active', 'latest_accepted_round':1306, 'round_in_progress':1307,
                      'originals_recoverable':len(layer['entries']), 'unchanged_support':len(layer['unchanged_support_sha256'])}))

if __name__ == '__main__':
    main()
