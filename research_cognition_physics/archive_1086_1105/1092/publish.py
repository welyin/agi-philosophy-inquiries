"""One-time root acceptance and reversible navigation publication of round1092."""
from pathlib import Path
from difflib import SequenceMatcher
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
BASE = ROOT/'research_cognition_physics'
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import relocation_original_bytes


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    acceptance = HERE/'acceptance.json'
    if acceptance.exists():
        raise SystemExit('Already published; refusing to overwrite acceptance.')
    reviewed = {
      '1092/proof.md':'315eeb02b080bd2e45d91258494d8998e0767e401fa20f571af4839f78c7c9ba',
      '1092/check.py':'cbc8782a8b92b40a436507ac9c1565194c25cc8e931113ac46b3bee097fd7f74',
      '1092/results.json':'db283c5afc90d3e68f5a72e6b5781ec4e7c56052e29b19b077958ba4feb04237',
      '1092/NEXT.md':'231286149a1473804fe1830828bd5adb3a15fd49e007c82381646740ca3817aa',
      'research_note_1092.md':'2ff1251d459e4b21c1d66dc0b453b0bb1ec31308be4859c194d728973134022e'}
    for path, value in reviewed.items():
        assert sha((STAGE/path).read_bytes()) == value, path

    layerpath = BASE/'_migration/active_documents/layer.json'
    original_layer = layerpath.read_bytes()
    layer = json.loads(original_layer)
    changes = []
    common = '''## 当前目标：验证可扩展共识与共同洛伦兹结构（1092起）

新目标active；[猜想](../猜想/可扩展共识与共同洛伦兹结构猜想.md)的有限能力与真实组合条件继续检验。F＋U＋C＋P及CO1—CO5共同保留，G未采用，目标内容未改。

**最新正式1092，累计3860。** [本轮报告](archive_1086_/research_note_1092.md)证明：统一维数有界、同组成有向且保全部资料的自然等距类型族，可以表示到一个最大维有限载体。由此直接复用1088紧性，定位1090全域目录何时必须真实扩容。B1—B4是本轮形式化分支，未自动成为新公理；尚未认证1091满足新猜想，更没有推出Lorentz。[证明与复算](archive_1086_/1092/proof.md) · [验收](archive_1086_/1092/acceptance.json) · [下一步](archive_1086_/1092/NEXT.md)。

用户新增“有效距离不能退化、因果先后须相容”的思路已进入[准入分析](archive_1086_/_admission/nondegenerate_distance_and_order.md)。下一步优先厘清独立输入选择时刻、已备资源及最快传递时间的量词，并与真实增员接口共同检查。390、282、385／412、1041相关旧结论直接复用；未重复立1093。速度变化本身不意味着因果次序颠倒。

研究继续在[archive_1086_](archive_1086_/README.md)保存，报告在阶段根目录、资产在编号目录。[新目标启动回执](archive_1086_/_shared/consensus_goal_start_20261010.json)保留；以下1086—1091结项属于原目标的限定反证，不代表新猜想已判定。

'''
    stagebody = '''## 当前目标：可扩展共识与共同洛伦兹结构猜想

新目标active，[启动回执](_shared/consensus_goal_start_20261010.json)保留。最新正式[1092](research_note_1092.md)、累计3860，新公理0。固定能力的有界相容类型族稳定到单有限载体，再复用1088紧性；新增的是这条连接，未证明或反驳整个[新猜想](../../猜想/可扩展共识与共同洛伦兹结构猜想.md)。

[1092证明](1092/proof.md) · [结果](1092/results.json) · [验收](1092/acceptance.json) · [下一步](1092/NEXT.md)。两代理独立终审通过工具消息交付，主线验收如实记录；没有伪称存在代理撰写的审阅文件。1091对新条件的联合资格仍未签署，原四原则与五公理不变。

新增[非退化测距与因果次序准入](_admission/nondegenerate_distance_and_order.md)接用户最新推测。优先区分已备资源下的信号时间与含准备的总任务时间，复用390、282、385／412、1041；不重复弱共识实验，不预占1093。

'''
    for entry in layer['entries']:
        name = entry['destination']
        p = ROOT/name
        old = p.read_bytes()
        before = relocation_original_bytes(old, entry)
        enc = 'utf-8-sig' if old.startswith(b'\xef\xbb\xbf') else 'utf8'
        current = old.decode(enc)
        new = current
        nl = '\r\n' if '\r\n' in current else '\n'
        if name in [f'research_cognition_physics/{f}' for f in
                    ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]:
            start = current.index('## 当前目标：')
            end = current.index('## 狭义相对论连接', start)
            new = current[:start]+common.replace('\n',nl)+current[end:]
        elif name.endswith('archive_1086_/README.md'):
            start = current.index('## 当前目标：')
            end = current.index('## 阶段目标', start)
            new = current[:start]+stagebody.replace('\n',nl)+current[end:]
        elif name.endswith('archive_1086_/文件索引.md'):
            addition = '''
- [1092报告](research_note_1092.md)：统一有界的自然类型族稳定到单载体；条件性扩容资格连接。
- [1092证明](1092/proof.md)、[复算代码](1092/check.py)、[结果](1092/results.json)、[主线验收](1092/acceptance.json)、[只读核验](1092/verify.py)与[下一步](1092/NEXT.md)。独立审阅通过工具消息交付，未创建代理审阅文件。
- [非退化测距与因果次序](_admission/nondegenerate_distance_and_order.md)：用户最新猜想的准入与旧结果复用，未占用1093。
'''.replace('\n',nl)
            pos = current.index('\n')+1
            new = current[:pos]+addition+current[pos:]
            new = new.replace('暂未占用1092。','现接1092条件性接口结论；未签整个新猜想。')
        if new == current:
            continue
        baseline = before.decode(enc)
        a, b = baseline.splitlines(keepends=True), new.splitlines(keepends=True)
        offsets = [0]
        for line in a:
            offsets.append(offsets[-1]+len(line))
        edits = []
        for tag,i,j,k,l in SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
            if tag != 'equal':
                edits.append({'start':offsets[i], 'end':offsets[j],
                              'before':''.join(a[i:j]), 'after':''.join(b[k:l]),
                              'kind':'round1092_conditional_bridge_and_user_distance_steering'})
        raw = new.encode(enc)
        entry.update(current_sha256=sha(raw), current_bytes=len(raw), edits=edits)
        assert relocation_original_bytes(raw, entry) == before
        changes.append((p,old,raw))
    assert len(changes) == 6
    assets = dict(reviewed)
    for name in ['verify.py','publish.py']:
        assets['1092/'+name] = sha((HERE/name).read_bytes())
    prior = {}
    for relative in ['1088/finite_carrier_proof.md','1090/proof.md','1091/proof.md']:
        p = STAGE/relative
        prior[p.relative_to(ROOT).as_posix()] = sha(p.read_bytes())
    evidence = {'round':1092, 'date':'2026-10-10', 'status':'PASS_CONDITIONAL_BRIDGE',
       'scope':'bounded directed isometric type family stabilizes; prior compactness reused',
       'entire_conjecture_decided':False, 'Lorentz_derived':False,
       'joint_FUCP_CO_and_new_conjecture_model_certified':False,
       'new_adopted_axioms':0, 'scientific_count_increment':0, 'scientific_count_total':3860,
       'assets_sha256':assets, 'prior_sources_sha256':prior,
       'independent_reviews':[{'agent':'/root/'+agent, 'status':'PASS_CONDITIONAL_BRIDGE',
           'delivery':'read_only_collaboration_message', 'reviewer_files_created':False,
           'source_sha256':reviewed, 'actual_check_run':{'groups':10,'maximum_residual':2.370726887964196e-15}}
           for agent in ['motion_source_review','typed_model_adversarial']],
       'review_delivery_note':'Root records actual read-only review messages. A reviewer file-write was rejected by automatic review citing the earlier read-only assignment; no retry or proxy review file is represented here.',
       'known_remaining':['real accession with whole motion and subsequent-task transport',
          'capacity to propagation bound', 'cross-resource uniformity', 'actual observer invariance'],
       'goal_status':'active'}
    assert layerpath.read_bytes() == original_layer
    for p,old,raw in changes:
        assert p.read_bytes() == old, str(p)
    acceptance.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    for p,old,raw in changes:
        p.write_bytes(raw)
    layerpath.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS', 'navigation_documents_updated':6,
                      'layer_sha256':sha(layerpath.read_bytes()), 'goal_status':'active'},ensure_ascii=False))


if __name__ == '__main__':
    main()
