"""Publish the reviewed D-P text; preserve all frozen round cores."""
from pathlib import Path
import hashlib
import json

STAGE=Path(__file__).resolve().parents[1]
BASE=STAGE.parent
PAPER=BASE/'认知动力学公理与狭义相对论基础动力学_阶段论文.md'

def replace_all(path, pairs):
    original=path.read_bytes()
    text=original.decode('utf-8-sig').replace('\r\n','\n')
    for before,after in pairs:
        assert text.count(before)==1,(str(path),before,text.count(before))
        text=text.replace(before,after)
    return original,text.encode('utf-8')

def main():
    receipt=STAGE/'_shared/adoption_record.json'
    assert not receipt.exists(),'Adoption already published; do not overwrite receipt'
    assert hashlib.sha256(PAPER.read_bytes()).hexdigest()=='bdc48b366ee8e074bc051dbfd279915a82a2e502b2025ba058c53c4d33d9970f'
    ledger_path=STAGE/'_shared/dependency_ledger.json'
    ledger=json.loads(ledger_path.read_text(encoding='utf-8'))
    assert ledger['completed_rounds']==[1116,1117,1118,1119]
    assert not ledger['in_progress_rounds'] and not ledger['new_adopted_axioms']
    for n in ledger['completed_rounds']:
        a=json.loads((STAGE/str(n)/'acceptance.json').read_text(encoding='utf-8'))
        assert a['status']=='PASS'
    changes={}
    changes[PAPER]=replace_all(PAPER,[
      ('阶段论文｜1116起研究｜2026年10月11日','阶段论文｜研究轮次1116—1119｜2026年10月11日'),
      ('**审阅稿。** 1116—1118已验收，1119共同实现证书正在核验；D-P尚未正式采用。本稿提出明确的采用文本和结项范围，最终状态以本阶段采用记录及逐项目标审计为准。',
       '**正式阶段论文。** 1116—1119四轮已验收，本文第3节的D-P正式加入认知公理体系；有限已知世界线版本、材料部门及未覆盖量词保持明确。数学和范围独立审阅通过，采用与依赖见[采用记录](archive_1116_/_shared/adoption_and_dependencies.md)，逐项验收见[阶段结项](archive_1116_/_shared/stage_completion.md)。历史报告中“候选／未采用”保留其原研究时点。'),
      ('本文提出完整过程平移接续公理D-P','本文采用完整过程平移接续公理D-P'),
      ('本文的拟结项对象是','本文的结项对象是'),
      ('D-P（本稿拟采用）','D-P（本文正式采用）'),
      ('D-P采用时成为认知公理体系的组成部分。','D-P现已成为认知公理体系的组成部分。'),
      ('本节证书由[1119报告](archive_1116_/research_note_1119.md)给出，当前待独立验收。',
       '本节证书由[1119报告](archive_1116_/research_note_1119.md)给出，已通过[独立审阅](archive_1116_/1119/independent_review.md)及主线验收。'),
      ('本阶段拟签收的结论为：','本阶段正式签收的结论为：'),
      ('当前1116—1118共17组数值代数检查通过，1119待验收。','1116—1119共22组数值代数检查通过，各轮独立审阅及主线验收完成。'),
      ('完整采用文本、来源分类和逐项验收将收于','完整采用文本、来源分类和逐项验收见'),
      ('本稿不自动启动场论、引力或下一研究目标。','本阶段不自动启动场论、引力或下一研究目标。'),
    ])
    adoption=STAGE/'_shared/adoption_and_dependencies.md'
    changes[adoption]=replace_all(adoption,[
      ('# D-P采用审议与狭义相对论基础动力学依赖收敛','# 认知动力学公理D-P：正式采用与依赖收敛'),
      ('2026-10-11。**拟采用文本，联合审阅尚未结束；当前D-P未采用。** 本页记录用户已授权的公理扩充如何落实，不再请求重复许可。正式状态待1119及论文审计后更新。',
       '2026-10-11。**D-P正式采用。** 按用户已授权的公理扩充，经1116—1119全部验收、论文独立数学审阅和目标范围审计，采用下述明确版本。新增一组动力学公理、四款操作条款；科学校准增量0。采用凭据见[adoption_record.json](adoption_record.json)。'),
      ('## 1. 拟采用版本','## 1. 正式采用版本'),
      ('## 5. 拟结项范围','## 5. 正式结项范围'),
      ('本页待1119独审及论文审计后正式采用；审议历史见[候选](dynamics_axiom_candidate.md)、[最终作用域预审](final_axiom_scope_audit.md)。',
       '1119独审和主线验收、[论文数学审阅](stage_paper_math_review.md)及[目标范围审阅](goal_completion_scope_review.md)均通过；本页据此正式采用。审议历史见[候选](dynamics_axiom_candidate.md)、[最终作用域预审](final_axiom_scope_audit.md)，其中“未采用”保留原研究时点，不覆盖本页当前状态。'),
    ])
    completion=STAGE/'_shared/stage_completion.md'
    changes[completion]=replace_all(completion,[
      ('2026-10-11。**四轮研究已验收；公理采用与论文最终签审进行中，目标仍active。** 本页先将原目标逐项对照完整证据，不因四轮复算通过就提前宣布结项。',
       '2026-10-11。**四轮、D-P采用和正式论文均已验收，主线签收本目标明确的有限有效域基础动力学结论。** 数学审阅、目标范围审阅及原目标逐项核验通过。目标工具回执在最终发布时另存；主线未以四轮数值通过代替解析或公理审查。'),
      ('本阶段的拟结项结论是：','本阶段的正式结论是：'),
      ('若上表中D-P和论文尚未完成采用签审，本页不把其状态写成已完成。共同证书的单个具体参数是相容见证；公理化定理量词由解析承担。',
       '上表D-P已按[正式采用记录](adoption_and_dependencies.md)加入体系，论文数学终核及目标范围审阅通过。共同证书的单个具体参数是相容见证；公理化定理量词由解析承担。'),
      ('本目标一旦正式结项就停止续轮','本目标正式结项后停止续轮'),
    ])
    for p,(old,new) in changes.items():
        assert p.read_bytes()==old,('Concurrent edit',str(p))
    for p,(_,new) in changes.items(): p.write_bytes(new)
    ledger.update(new_adopted_axioms=['D-P'],adopted_clause_count=4,
        PT_status='D-P finite known-worldline version adopted; PT2 sector and PT3 provenance classified',
        key_open_connections=[],stage_status='verified_ready_for_goal_completion',
        theorem_scope='stable positive finite rest types; fixed type continuous action with original ports; certified full effective boundaries',
        remaining_freedom=['absolute mass and material spectra','couplings and permitted interactions','stability and bound-state existence','coverage beyond declared material and readout domain'],
        final_paper=PAPER.name)
    ledger_path.write_text(json.dumps(ledger,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    evidence={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest()
      for p in [PAPER,adoption,STAGE/'_shared/stage_paper_math_review.md',STAGE/'_shared/goal_completion_scope_review.md',STAGE/'1119/acceptance.json']}
    record=dict(status='ADOPTED',date='2026-10-11',axiom_group='D-P',clauses=4,
        user_authorization='Explicit new goal: derive PT from current axioms where possible, otherwise openly axiomatize justified extensions',
        version='finite known-worldline types; fixed identifications; original ports retained',
        formal_text=str(PAPER.relative_to(BASE))+'#3-D-P',
        axiom_group_is_not_proven_minimal_or_independent=True,
        massive_sector_not_universal_axiom=True,SR_S_optional=True,G_adopted=False,
        reviewed_draft_sha256='bdc48b366ee8e074bc051dbfd279915a82a2e502b2025ba058c53c4d33d9970f',
        publication_changes='adoption, completion status, totals and source links only; mathematical clauses unchanged',
        evidence_sha256=evidence,scientific_calibration_increment=0,
        goal_tool_status='active_pending_final_receipt',external_human_review=False)
    receipt.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='ADOPTED',group='D-P',clauses=4,paper_sha256=evidence[str(PAPER.relative_to(BASE))]),ensure_ascii=False))

if __name__=='__main__':main()
