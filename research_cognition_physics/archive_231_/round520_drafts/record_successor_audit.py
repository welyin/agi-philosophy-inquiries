"""Save successor deduplication without opening a new scientific round."""
from pathlib import Path
import hashlib
import json
import re

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
checked=json.loads((HERE/'round520_integration_checks.json').read_text('utf8'))
assert checked['all_reported_checks_passed']
assert checked['total_protected_evidence_hashes']==1008
additions={
HERE/'spatial_premise_closure_audit.md': '''

### 158.1 520后立即接续的候选去重

独立审核及父级复核确认三种接法的范围。第一，两个估计器实际收到同一完整记录、使用同一H与标签字典时，从不同固定满秩猜测也会趋同；将520的d换成初猜最小本征值的倒数，再用三角不等式即可。这不对所有满秩猜测给统一常数，也不证明仅凭各自局部记录自动形成共识。第二，按已有记录选择的同一后继定位仪器直接由CPTP收缩继承误差，不重跑实验；若未知旧R重新参与，数学联合收缩不等于计算器已知R的状态。

第三，把条件估计图态的平均路径长度直接叫作位置仍失败：[468](research_note_468.md)适用固定逐标签度数、全部标签的任意图人口；[489](research_note_489.md)又给完整六主体标尺的统一小尺度障碍。因此学习到当前态不会让该均距接口获得例外。467的单树限制不能滥用于均值，469允许选定端点的任意维数近似也继续保留。

这三项均是既有结果的合成，不开521，不增加科学文件或检查数。真正待补的是同一实际过程为何选出某端点族、尺度与定位报告，并使其在新端点加入、路径组合和尺度改变时仍有效。只要求声明的定位菜单，不额外要求位置决定全部私人未来。520有实质进展，目标保持活动；下一项不再以更多滤波曲线或两个相同算法趋同充数。
''',
RESEARCH/'research_direction.md': '''

**520后候选去重：** 已审共同完整记录下不同满秩猜测趋同、接同一未来菜单的误差继承，均为520及成熟滤波工具直接推论，不开521。468／489仍限制将完整六主体的条件图均距直接认作局部空间标尺。下一项须减少端点族、尺度或定位报告的选择输入，详见[审计158.1](archive_231_/spatial_premise_closure_audit.md)。
''',
RESEARCH/'RESEARCH_STATE.md': '''

**520后立即接续：** [审计158.1](archive_231_/spatial_premise_closure_audit.md)完成三项候选去重；共同记录估计趋同与未来通道收缩不重复编号，条件图均距仍受468／489范围限制。当前520／2552，1008份证据，不新增521；有来源的宏观位置及三维仍是下一实质缺口。
'''}
before={p:p.read_bytes() for p in additions}
snapshot=HERE/'round520_drafts/navigation_before_successor_audit'
snapshot.mkdir(exist_ok=False)
planned={}
for p,tail in additions.items():
    raw=before[p]
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    label=('research_' if p.parent==RESEARCH else 'archive_')+p.name
    with (snapshot/label).open('xb') as f:f.write(raw)
    planned[p]=(raw.decode(enc).replace('\r\n','\n')+tail).replace('\n',nl).encode(enc)
for p,raw in before.items():assert p.read_bytes()==raw,str(p)
for p,raw in planned.items():
    assert p.read_bytes()==before[p],str(p)
    p.write_bytes(raw)
manifest={str(p.relative_to(ROOT)):{'before':hashlib.sha256(before[p]).hexdigest(),
    'after':hashlib.sha256(raw).hexdigest()} for p,raw in planned.items()}
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump(manifest,f,ensure_ascii=False,indent=2)
# This is a navigation check, not an extra physics experiment.
import sys
sys.path.insert(0,str(HERE))
import verify_interaction_rounds as core
docs=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
      RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
      HERE/'three_dimensional_four_conditions_review.md',HERE/'research_note_520.md']
links=0
for doc in docs:
    for link in core.link_parser()(doc.read_text('utf-8-sig')):
        assert (doc.parent/link).resolve().exists(),(doc,link)
        links+=1
for name,sha in json.loads((HERE/'research_round_520_checks.json').read_text('utf8'))['new_file_hashes'].items():
    assert core.digest(HERE/name)==sha,name
result=dict(round=520,numbered_tests=2552,protected_evidence=1008,
    new_numbered_round=False,new_scientific_checks=0,navigation_documents_checked=len(docs),
    local_links_checked=links,broken_links=0,independent_successor_review_completed=True,
    all_reported_checks_passed=True)
with (HERE/'round520_successor_navigation_checks.json').open('x',encoding='utf8') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result))
