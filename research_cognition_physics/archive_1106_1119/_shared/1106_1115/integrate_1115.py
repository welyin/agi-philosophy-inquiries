"""Publish the reviewed domain-realization certificate, preserving frozen history."""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes
from run_research_motion_active import MotionLayout

BASE = ROOT / 'research_cognition_physics'
STAGE = BASE / 'archive_1106_'
FOLDER = STAGE / '1115'
sha = lambda b: hashlib.sha256(b).hexdigest()
originals, updates = {}, {}

def read(p):
    p = Path(p)
    originals[p] = p.read_bytes()
    return originals[p].decode('utf-8-sig')

def put(p, text):
    p = Path(p)
    if p not in originals:
        originals[p] = p.read_bytes() if p.exists() else None
    if originals[p] and b'\r\n' in originals[p]:
        text = text.replace('\r\n', '\n').replace('\n', '\r\n')
    updates[p] = text.encode('utf-8')

def putjson(p, obj):
    put(p, json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

assert not (FOLDER / 'acceptance.json').exists(), 'Already published'
result = json.loads(read(FOLDER / 'results.json'))
assert result['numerical_check_groups'] == 6 and result['scientific_calibration'] == 0
review = read(FOLDER / 'independent_review.md')
assert 'PASS_DECLARED_DOMAIN_COMPLETE_PARENT_PROCESS_REALIZATION' in review
core = {}
for name in ['proof.md', 'check.py', 'results.json', 'NEXT.md', '../research_note_1115.md']:
    p = FOLDER / name
    read(p)
    core[name] = sha(originals[p])
    assert core[name] in review, ('Unreviewed core', name)
for name, expected in result['source_sha256'].items():
    assert sha((STAGE / name).read_bytes()) == expected, ('Changed source', name)

put(FOLDER / 'mainline_review.md', '''# 1115主线验收：声明有效域共同实现证书通过

2026-10-11。主线实读完整证明、编号报告、保存结果与独立审阅，核对最终五核心及九份冻结来源哈希，接受**精确认知父族在声明有限有效域中的共同实现证书**。

父族保留完整FUCP；材料层实现预先定义的逻辑仪器、原生钟／等待、完整物质散射及其允许组合。两层通过同一真实端口、原生H／四动量、完整输出和事件字典连接。原生钟、反冲、环境和程序保留全部量子对象身份，未认证的母体操作明确列为实现范围以外，未以受限交换子子群替代完整C。

式(3)覆盖所有合法未知输入与任意参考；式(11)的菜单等同要求指定实际接收因子上的局部输出，不能依赖分布式全局解码。参考运输、守恒和全部纳入影响分别承接1114；它们不是从通道近似误差小自动获得。内部钟记录来自已有耦合，不授予中途外接无扰读取或免费复位。

独立审阅者实际运行默认入口，6组通过且保存结果一致。主线检查了其执行证据，未为相同字节重复大型科学计算。科学校准增量0，累计3860；没有外部人类同行审稿或经验检验。

这是轮次验收。阶段是否完成另按原目标逐条审核；本收据时目标仍active。1113全材料固定窗口覆盖仍为可选加强项，不自动作为本证书前置。

[证明](proof.md) · [结果](results.json) · [独立审阅](independent_review.md) · [报告](../research_note_1115.md) · [目标范围审计](../_shared/after_1114_goal_scope_audit.md)。
''')
putjson(FOLDER / 'acceptance.json', dict(
    schema='round1115_acceptance_v1', date='2026-10-11', status='PASS',
    author='/root/influence_bridge_review', independent_reviewer='/root/typed_model_adversarial',
    mainline_reviewer='/root', external_human_review=False, core_sha256=core,
    review_sha256=sha(originals[FOLDER / 'independent_review.md']),
    independent_default_run=dict(groups=6, saved_results_match=True, exit_code=0),
    mainline_verification='complete proof, input/output roles, event-menu reflection, source/core hashes and independent execution evidence',
    scope='精确认知父族在声明有限有效域中的完整过程共同实现',
    new_adopted_axioms=0, scientific_calibration_increment=0, scientific_count_total=3860,
    goal_status='active', goal_complete=False))

summary = ('**已完成1106—1112、1114—1115，共九份报告。** '
 '[1115](archive_1106_/research_note_1115.md)将精确量子父族、原CO位置、实际钟、物质反冲、全部资源和带事件任务接入同一完整过程字典；'
 '承接1112的原端口逻辑编译及1114的整网参考与全部已纳入影响边界。阶段共59组构造检查，新增科学校准0，累计3860。'
 '结论是声明有限有效域中的共同实现，未把原生材料受限菜单改称完整C。'
 '按[原目标范围审计](archive_1106_/_shared/after_1114_goal_scope_audit.md)进行逐项阶段验收；'
 '[1113](archive_1106_/1113/boundary_port_admission.md)仍为可选加强覆盖候选，不计为完成报告，也不自动作为本阶段前置。'
 '当前目标保持active，完成状态尚待原目标验收。历史复算使用 scripts/run_research_motion_active.py。')
lp = BASE / '_migration/motion_active_1106/layer.json'
layer = json.loads(read(lp))
marker = '## 历史导航（截至1105，保留原研究时点）'
for entry in layer['entries']:
    p = ROOT / entry['destination']
    if p.name not in {'README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md'}:
        continue
    before = read(p)
    baseline = relocation_original_bytes(originals[p], entry)
    changed, count = re.subn(r'\*\*已完成1106—1112及1114，共八份报告。\*\*[^\r\n]+', summary, before, count=1)
    assert count == 1
    assert before.split(marker, 1)[1] == changed.split(marker, 1)[1]
    put(p, changed)
    entry.update(current_sha256=sha(updates[p]), current_bytes=len(updates[p]),
                 edits=text_edits(baseline.decode('utf-8-sig'), updates[p].decode('utf-8')))
    assert relocation_original_bytes(updates[p], entry) == baseline
layer.update(scope='Live navigation and SR paper; nine accepted reports 1106-1112 and 1114-1115; stage audit pending',
             new_formal_reports=9, latest_completed_round=1115)
putjson(lp, layer)

p = STAGE / 'README.md'
before = read(p)
changed, count = re.subn(r'已完成1106—1112及1114八份报告[^\r\n]+',
    '已完成1106—1112、1114—1115九份报告，共59组构造检查，科学校准增量0，累计3860。'
    '[1115](research_note_1115.md)已通过独立审阅和主线验收，给声明有效域内的完整父过程共同实现。'
    '阶段原目标逐项验收正在进行，目标保持active；1113仍为可选覆盖候选，不计完成。', before, count=1)
assert count == 1
changed, count = re.subn(r'下一项按\[1114下一步\][^\r\n]+',
    '下一项按[1114后目标范围审计](_shared/after_1114_goal_scope_audit.md)与[1115证书](1115/proof.md)，'
    '核对阶段全部要求、同一对象和实际过程覆盖；不把原生全母体固定期限覆盖自动升级为门槛。'
    '旧轮的“下一步”保留原研究时点，本段给当前执行优先级。', changed, count=1)
assert count == 1
put(p, changed)

p = STAGE / '文件索引.md'
before = read(p)
needle = '已完成1106—1112及1114共八份报告，53组构造检查；1113候选不计入已完成轮次。最新完成1114，完整共同模型验收继续。'
assert needle in before
row = '|1115|[完整父对象与有限材料的共同实现](research_note_1115.md)：完成，独审及主线PASS|[证明](1115/proof.md) · [结果](1115/results.json) · [验收](1115/mainline_review.md) · [下一步](1115/NEXT.md)|'
pos = before.index(needle)
changed = before[:pos].rstrip('\r\n') + '\n' + row + '\n\n' + '已完成1106—1112、1114—1115共九份报告，59组构造检查；1113候选不计入完成轮次。阶段验收继续，目标active。' + before[pos + len(needle):]
changed += '\n- [1114后原目标范围审计](_shared/after_1114_goal_scope_audit.md)\n'
put(p, changed)

p = STAGE / '_shared/dependency_ledger.json'
ledger = json.loads(read(p))
ledger.update(completed_rounds=[1106,1107,1108,1109,1110,1111,1112,1114,1115],
              latest_completed_round=1115, new_completed_rounds=9, numerical_check_groups=59)
ledger['progress'].update(meaningful_new_result=1115, no_progress_streak=0, goal_complete=False,
    previous_goal_turn='progress: 1115 complete parent/domain/port/resource realization independently verified and published')
for c in ledger['open_connections']:
    if c['id'] == 'JOINT':
        c.update(status='declared_domain_complete_process_realization_verified',
                 completed_delta='1115给完整父对象、事前实验域、全资源过程等距及同一事件任务菜单；未签全材料精确C',
                 next='按原目标逐项阶段验收；1113更强全材料定时覆盖非自动前置')
ledger['next_admission'] = dict(path='research_cognition_physics/archive_1106_/_shared/after_1114_goal_scope_audit.md',
    result='1115限定共同实现证书已验收，阶段全部目标尚待独立审计', next='逐项目标证据、依赖消冗余与冲突审计',
    completed_round=False, additional_device_expansion_priority=False,
    complete_C_of_all_physical_mother_objects_proved=False)
ledger['realization_certificate'] = dict(round=1115,
    status='verified_declared_domain_complete_process_relation',
    exact_parent_FUCP=True, native_mother_full_C_claimed=False,
    all_resources_retained=True, receiver_factorization_required=True,
    new_axiom=False, whole_physical_SR_claimed=False)
putjson(p, ledger)

for p, raw in originals.items():
    assert (p.read_bytes() if p.exists() else None) == raw, ('Concurrent change', p)
for p, raw in updates.items():
    p.write_bytes(raw)
compat = MotionLayout().verify_motion_layer()
out = FOLDER / 'integration_checks.json'
assert not out.exists()
out.write_text(json.dumps(dict(schema='round1115_publication_v1',date='2026-10-11',
    core_sha256=core, frozen_source_hashes_match=9, compatibility=compat,
    historical_navigation_tails_unchanged=True, completed_reports=9, check_groups=59,
    candidate1113_not_counted=True, goal_complete=False), ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(json.dumps(dict(status='PASS', completed=1115, reports=9, groups=59, compatibility=compat), ensure_ascii=False))
