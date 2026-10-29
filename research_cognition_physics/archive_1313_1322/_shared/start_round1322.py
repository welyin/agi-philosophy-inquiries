"""Record admission once; retain the single reversible active navigation layer."""
import copy
import json
from publish_progress_after1316 import ROOT, STAGE, LAYER, sha, encoded, atomic
from research_layout import relocation_original_bytes
from run_research_interaction_generation_active import InteractionGenerationLayout


def main():
    receipt = STAGE / '_shared/round1322_start_receipt.json'
    assert not receipt.exists()
    InteractionGenerationLayout().verify_interaction_navigation()
    raw = LAYER.read_bytes()
    layer = json.loads(raw)
    assert layer['latest_accepted_round'] == 1321
    assert layer['next_round'] == 1322 and layer['in_progress_round'] is None
    admission = STAGE / '_admission/round1322_admission.md'
    assert admission.is_file()
    original_layer = copy.deepcopy(layer)
    updates = []
    for entry in layer['entries']:
        p = ROOT / entry['destination']
        before = p.read_bytes()
        original = relocation_original_bytes(before, entry)
        prefix = entry['edits'][0]['after']
        start = prefix.index('下一编号1322先审')
        end = prefix.index('\n\n上一阶段', start)
        prefix = prefix[:start] + ('[1322已准入](archive_1313_/_admission/round1322_admission.md)，'
            '正在检验同一左门户、非零H参考与强交接任务中的右弱角色替代。'
            '保持实际有限父过程和剩余保护，分别登记未破缺表示、高次质量作用及右规范源变化；'
            '不把关系投影当作真实弱流消失。科学验收仍止于1321，A/B/C/J目标继续。') + prefix[end:]
        enc = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf8'
        content = original.decode(enc)
        at = entry['edits'][0]['start']
        after = (content[:at] + prefix + content[at:]).encode(enc)
        entry['edits'][0]['after'] = prefix
        entry.update(current_sha256=sha(after), current_bytes=len(after))
        assert relocation_original_bytes(after, entry) == original
        updates.append((p, before, after))
    layer['in_progress_round'] = 1322
    layer['startup_state'] = 'round1322 admitted; scientific acceptance remains at1321'
    updates.append((LAYER, raw, encoded(layer)))
    record = dict(status='STARTED_NOT_ACCEPTED', round=1322, new_axioms=0,
                  goal_status='active', admission_sha256=sha(admission.read_bytes()),
                  previous_layer_sha256=sha(raw), current_layer_sha256=sha(encoded(layer)),
                  previous_layer=original_layer)
    for p, before, _ in updates:
        assert p.read_bytes() == before, str(p)
    written = []
    try:
        for p, before, after in updates:
            atomic(p, after)
            written.append((p, before, after))
        with receipt.open('xb') as stream:
            stream.write(encoded(record))
    except BaseException:
        for p, before, after in reversed(written):
            if p.read_bytes() == after:
                atomic(p, before)
        raise
    print(json.dumps({k: v for k, v in record.items() if k != 'previous_layer'}))


if __name__ == '__main__':
    main()
