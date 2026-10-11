"""One-time publication of accepted rounds 1112 and 1114; no science rerun.

Maintain the existing reversible navigation layer and preserve historical
tails and the current SR paper. Refuse changed concurrent inputs.
"""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes
from run_research_motion_active import MotionLayout

BASE = ROOT/'research_cognition_physics'
STAGE = BASE/'archive_1106_'
sha = lambda raw: hashlib.sha256(raw).hexdigest()
updates = {}
originals = {}


def read(path):
    path=Path(path)
    raw=path.read_bytes()
    originals[path]=raw
    return raw.decode('utf-8-sig')


def set_text(path, text):
    path=Path(path)
    if path not in originals:
        originals[path]=path.read_bytes() if path.exists() else None
    if originals[path] and b'\r\n' in originals[path]:
        text=text.replace('\r\n','\n').replace('\n','\r\n')
    updates[path]=text.encode('utf-8')


def set_json(path, obj):
    set_text(path, json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


checks = {}
for number, groups, author, reviewer in [
    (1112,8,'/root/influence_bridge_review','/root/sr_next_goal_audit'),
    (1114,6,'/root','/root/typed_model_adversarial')]:
    folder=STAGE/str(number)
    assert not (folder/'acceptance.json').exists(), f'{number} already published'
    result=json.loads(read(folder/'results.json'))
    assert result['numerical_check_groups']==groups and result['scientific_calibration']==0
    review=read(folder/'independent_review.md')
    assert 'PASS' in review
    core={p:sha((folder/p).read_bytes()) for p in
          ['proof.md','check.py','results.json','NEXT.md',f'../research_note_{number}.md']}
    for value in core.values():
        assert value in review, f'Unreviewed {number} core: {value}'
    scope=('持续运行装置上的完整逻辑装卸、任意参考、完整记录与原CO端口固定窗口编译'
           if number==1112 else '明确有限材料准备族的整网参考闭合、全部开放影响未来支持及SR-F')
    set_text(folder/'mainline_review.md',f'''# {number}主线验收：条件通过

2026-10-11。主线实读完整证明、保存结果及[独立审阅](independent_review.md)，接受其限定结论：{scope}。

独立审阅者实际运行默认入口，{groups}组通过且保存结果一致。主线核验最终五核心哈希与审阅对应；未为文案调整重复相同科学计算。本轮没有外部人类同行审稿或经验数据。

有效材料权限和普遍认知公理分列。完整输出、原端口与实际窗口不能被约化逻辑通道替代；未认证的非标量材料全仪器、实际跨部门输入及整个SR仍未签收。当前目标保持active。新增科学校准0，累计3860。

[证明](proof.md) · [结果](results.json) · [编号报告](../research_note_{number}.md) · [下一步](NEXT.md)。
''')
    receipt=dict(schema=f'round{number}_acceptance_v1',date='2026-10-11',status='PASS',
      author=author, independent_reviewer=reviewer, mainline_reviewer='/root',
      external_human_review=False, core_sha256=core,
      review_sha256=sha(originals[folder/'independent_review.md']),
      independent_default_run=dict(groups=groups,saved_results_match=True,exit_code=0),
      mainline_verification='proof, scope and actual final-file hashes; independent execution evidence inspected',
      scope=scope,new_adopted_axioms=0,scientific_calibration_increment=0,
      scientific_count_total=3860,goal_status='active',goal_complete=False)
    set_json(folder/'acceptance.json',receipt)
    checks[str(number)]=dict(groups=groups,core_sha256=core)

summary=('**已完成1106—1112及1114，共八份报告。** '
 '[1112](archive_1106_/research_note_1112.md)给同一持续运行装置的完整逻辑装卸、记录及原CO端口固定窗口编译；'
 '[1114](archive_1106_/research_note_1114.md)给明确有限材料族的整网参考运输、全部已纳入影响的共同因果界，以及跨所有允许有限准备的SR-F。'
 '两轮独立审阅通过；阶段共53组构造检查，新增科学校准0，累计3860。'
 '[1113](archive_1106_/1113/boundary_port_admission.md)仍是未立项的接口候选，不计为已完成报告。'
 '下一步核验共同模型必要的完整材料覆盖与原端口接口：已有守恒逼近不能自动保证增加参考资源后原材料的端口和期限不变；先审量词与同一入出编码，停止扩建器件。'
 '精确操作层与有限材料层继续按[有效实现字典](archive_1106_/_shared/effective_realization_dictionary.md)核验。'
 '整个实际SR目标保持active，历史复算使用 scripts/run_research_motion_active.py。')
layer_path=BASE/'_migration/motion_active_1106/layer.json'
layer=json.loads(read(layer_path))
marker='## 历史导航（截至1105，保留原研究时点）'
for entry in layer['entries']:
    path=ROOT/entry['destination']
    if path.name not in {'README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'}:
        continue
    before=read(path)
    original=relocation_original_bytes(originals[path],entry)
    changed,count=re.subn(r'\*\*1106—1111六份报告已完成并验收。\*\*[^\r\n]+',summary,before,count=1)
    assert count==1
    assert before.split(marker,1)[1]==changed.split(marker,1)[1]
    set_text(path,changed)
    entry['current_sha256']=sha(updates[path])
    entry['current_bytes']=len(updates[path])
    entry['edits']=text_edits(original.decode('utf-8-sig'),updates[path].decode('utf-8'))
    assert relocation_original_bytes(updates[path],entry)==original
layer.update(scope='Live navigation and current SR paper, eight completed reports 1106-1112 and 1114; 1113 candidate only',
             new_formal_reports=8,latest_completed_round=1114)
set_json(layer_path,layer)

p=STAGE/'README.md'
text=read(p)
text=re.sub(r'已完成1106—1111六份报告[^\r\n]+',
 '已完成1106—1112及1114八份报告，共53组构造检查，新增科学校准0，累计3860。'
 '[1112](research_note_1112.md)补上完整逻辑装卸及原CO端口；[1114](research_note_1114.md)补上整网实际参考族与全部已纳入影响的共同边界。'
 '均经独立审阅；1113只有候选准入，未计完成。阶段目标继续active。',text,count=1)
text=re.sub(r'下一项按\[共同模型收束准入\][^\r\n]+',
 '下一项按[1114下一步](1114/NEXT.md)及[1113原端口准入](1113/boundary_port_admission.md)，'
 '先审完整材料覆盖的必要量词，再核同一入出编码、固定原端口和窗口。参考资源增加后不得将旧材料输出接口悄悄改名；'
 '不以未纳入的更强闭合或普遍自治G扩大阶段门槛。精确公理和材料实现关系仍分列，整个目标尚未验收。',text,count=1)
set_text(p,text)

p=STAGE/'文件索引.md'
text=read(p)
anchor='1106—1111六份报告完成并验收，共39组数值检查。1106按原预留编号补齐，最新完成1111。下一项共同模型闭合只记准入，不将准备或审计冒充完成轮次。'
assert anchor in text
rows='''|1112|[完整逻辑装卸与原CO端口编译](research_note_1112.md)：完成，独审PASS|[证明](1112/proof.md) · [结果](1112/results.json) · [验收](1112/mainline_review.md) · [下一步](1112/NEXT.md)|
|1113|未立项候选，尚无正式报告|[固定窗口准入](1113/admission.md) · [原端口接口](1113/boundary_port_admission.md)|
|1114|[整网参考闭合与全部影响边界](research_note_1114.md)：完成，独审PASS|[证明](1114/proof.md) · [结果](1114/results.json) · [验收](1114/mainline_review.md) · [下一步](1114/NEXT.md)|'''
position=text.index(anchor)
text=(text[:position].rstrip('\r\n')+'\n'+rows+'\n\n'
      +'已完成1106—1112及1114共八份报告，53组构造检查；1113候选不计入已完成轮次。最新完成1114，完整共同模型验收继续。'
      +text[position+len(anchor):])
assert '|1114|' in text
set_text(p,text)

p=STAGE/'_shared/dependency_ledger.json'
ledger=json.loads(read(p))
ledger.update(completed_rounds=[1106,1107,1108,1109,1110,1111,1112,1114],
              new_completed_rounds=8,numerical_check_groups=53)
ledger['progress'].update(meaningful_new_result=1114,no_progress_streak=0,goal_complete=False,
 previous_goal_turn='progress: 1112 logical type/port compilation and 1114 whole-preparation reference/influence closure independently verified')
for c in ledger['open_connections']:
    if c['id']=='REFERENCE':
        c.update(status='conditional_whole_preparation_family_closed',completed_delta='1114在明确整网Lorentz预部署和有限接线权限下，全部事件/端口/钟/四动量/完整输入输出及逆/有限词保持，真实任务菜单AP=P')
    elif c['id']=='INFLUENCE':
        c.update(status='declared_material_family_future_support_proved',completed_delta='1114对同一声明材料族的全部开放局部干预、完整记录及反馈证明未来支持；1106单源完整任务给c*=1及跨全准备并集SR-F，不靠噪声精确放大',next='新纳入材料若扩大干预菜单，须同时核其未来局部支持')
    elif c['id']=='JOINT':
        c.update(status='complete_material_port_interface_audit',next='区分阶段必需覆盖与更强全母体闭合；同一完整入出字典是否在参考资源增加后保持原端口和固定窗口',completed_delta='1112原CO逻辑仪器及真实装卸已接入1114整网参考/影响共同族；一般非标量完整材料的实现范围仍待核验')

for ident, rnd, motivation, scope, permission, witness in [
 ('K1112',1112,'共同资料和记录必须从运行主体完整移交到原定外部端口',
  '同一HFK、简并逻辑工作槽、全部未知逻辑输入及任意参考、有限原CO仪器和同窗接续',
  'K1—K4真实vac载体、正时隙dressed门、完整装卸及事前有限线路权限',
  '全864维装卸、完整55296乘8会合链；8组构造检查和独审'),
 ('R1114',1114,'参考比较须保持完整操作的兑现能力和所有内部资源',
  '1106/1111/1112已认证实际部门的有限网络族，全部开放干预/记录/反馈，任意有限Lorentz词',
  '整网主动预部署权限及原端口/未来因果支持相容的有限闭合',
  '全事件/轨道/钟相位/四动量运输，单源Q轨道与跨全族SR-F；6组校准和独审')]:
    ledger['new_conditions'].append(dict(id=ident,source_round=rnd,cognitive_motivation=motivation,
       atomic_compositional_or_macro_layer='明确有限有效材料及共同实现关系',quantifiers=scope,
       effective_domain=scope,adoption_status='explicit_effective_material_permission_not_new_axiom',
       necessity_or_sufficiency='该连接的充分条件；未证明必要、最小或唯一',dependencies=permission,
       excluded_models='仅有内部逻辑矩阵、边缘态或被动坐标而冒称实际完整端口/任务权限的实现',
       joint_witness=witness,alternatives='其它保持同一完整接口的局部材料与参考准备',
       finite_error_consequences='完整过程界保同一输出字典；守恒和因果支持结构核验，非由范数小推出',
       all_physical_SR_proved=False))
ledger['next_admission']=dict(path='research_cognition_physics/archive_1106_/1113/boundary_port_admission.md',
 result='逻辑原端口与整网参考接口已条件完成；一般完整材料的必要覆盖和同一固定接口仍待审计',
 next='最小必要量词审计及两态非共移材料、两种bank规模的保原端口筛查',
 completed_round=False,additional_device_expansion_priority=False,
 complete_C_of_all_physical_mother_objects_proved=False)
set_json(p,ledger)

# Perform all optimistic checks before the first mutation.
for path,raw in originals.items():
    assert (path.read_bytes() if path.exists() else None)==raw, f'Concurrent change: {path}'
for path,raw in updates.items():
    path.write_bytes(raw)
compat=MotionLayout().verify_motion_layer()
receipt=dict(schema='rounds1112_1114_publication_v1',date='2026-10-11',rounds=checks,
            compatibility=compat,historical_navigation_tails_unchanged=True,
            existing_scientific_code_not_modified=True,completed_reports=8,check_groups=53,
            candidate1113_not_counted=True,goal_complete=False)
out=STAGE/'1114/integration_checks.json'
assert not out.exists()
out.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'PASS','completed':[1112,1114],'stage_checks':53,'compatibility':compat},ensure_ascii=False))
