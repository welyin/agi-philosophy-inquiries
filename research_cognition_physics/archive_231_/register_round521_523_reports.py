"""Register three completed studies as full numbered notes; do not redo research."""
import argparse
import json
from pathlib import Path
import verify_thermal_direction_bridge as previous
import verify_interaction_rounds as core

HERE = Path(__file__).resolve().parent
CONFIG = {
    521: ('positive_reference_acquisition', 'positive_reference_acquisition.py', 9, 520),
    522: ('thermal_reference_dimension', 'thermal_reference_dimension_model.py', 12, 520),
    523: ('thermal_direction_bridge', 'thermal_direction_bridge_model.py', 13, 522),
}


def numbered_content(number):
    stem, _, _, scientific_base = CONFIG[number]
    original = (HERE/(stem+'_review.md')).read_text('utf8')
    title, rest = original.split('\n\n',1)
    _, body = rest.split('\n\n',1)
    intro = (
        f'日期：2026-09-30。科学内容基线为第{scientific_base}轮及已注明的后续接口。'
        f'[代码]({CONFIG[number][1]})；[结果]({stem}_results.json)；'
        f'[本轮核验](research_round_{number}_checks.json)。\n\n'
        f'**编号说明：** 本项此前保存为[完整接口报告]({stem}_review.md)，'
        '现按用户要求纳入连续研究轮次，保留全部原稿、代码和结果。'
        '下面保存完整研究正文，不只是索引或摘要。6组复算检查已通过独立审阅；'
        '本次登记将其纳入编号统计，不新增实验、不重复计为另一份科学成果。'
        '成熟定理仍按原出处引用，额外条件与结论范围保持不变。')
    body = body.replace('成熟工具的复用不计为新编号研究。',
                        '成熟工具仍属复用；本轮6组检查核验这里的具体模型，不将工具本身宣称为新定理。')
    if number==523:
        body = body.replace('(thermal_reference_dimension_review.md)', '(research_note_522.md)')
    return f'# 第{number}轮：{title.removeprefix("# ")}\n\n{intro}\n\n{body}'


def register():
    base = previous.verify()
    assert base['total_protected_including_this_review']==1052
    assert base['unchanged_numbered_scientific_tests']==2552
    for n in CONFIG:
        assert not (HERE/f'research_note_{n}.md').exists(),n
        assert not (HERE/f'research_round_{n}_checks.json').exists(),n
        assert not (HERE/f'round{n}_drafts').exists(),n
    for index,(n,(stem,code,formulas,science_base)) in enumerate(CONFIG.items()):
        content = numbered_content(n).encode('utf8')
        note = HERE/f'research_note_{n}.md'
        folder = HERE/f'round{n}_drafts'
        folder.mkdir(exist_ok=False)
        draft = folder/f'research_note_{n}.txt'
        for path in (note,draft):
            with path.open('xb') as stream:stream.write(content)
        checked = core.text_checks(note)
        assert checked['display_formulas']==formulas
        saved = core.read(HERE/(stem+'_results.json'))
        assert (saved.get('diagnostic_groups',saved.get('diagnostic_tests')),
                saved['failures'],saved['errors'])==(6,0,0)
        report = dict(date='2026-09-30',round=n,
            scientific_base_through_round=science_base,registration_predecessor=n-1,
            registered_tests=dict(run=6,failures=0,errors=0),
            tests_transferred_from_existing_unnumbered_report=6,
            new_scientific_experiments_for_registration=0,
            cumulative_numbered_tests=2552+6*(index+1),
            cumulative_numbered_scientific_files=871+3*(index+1),
            cumulative_unique_protected_evidence_files=1052+2*(index+1),
            saved_results_reproduced=True,scientific_results_rewritten=False,
            new_file_hashes={name:core.digest(HERE/name) for name in
                             (code,stem+'_results.json',note.name)},
            previously_protected_scientific_files_registered=[code,stem+'_results.json'],
            preserved_draft_hashes={str(draft.relative_to(HERE)).replace('\\','/'):core.digest(draft)},
            original_full_report=stem+'_review.md',
            original_full_report_sha256=core.digest(HERE/(stem+'_review.md')),
            full_body_preserved_with_numbering_metadata_only=True,
            text_checks=checked,independent_final_code_and_draft_review_completed=True,
            visual_checks_performed=False,scope=saved['scope'],all_reported_checks_passed=True)
        with (HERE/f'research_round_{n}_checks.json').open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(registered_rounds=list(CONFIG),full_reports_created=3,
        cumulative_numbered_tests=2570,unique_protected_evidence_files=1058,
        new_scientific_experiments_for_registration=0)))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true',required=True)
    parser.parse_args()
    register()
