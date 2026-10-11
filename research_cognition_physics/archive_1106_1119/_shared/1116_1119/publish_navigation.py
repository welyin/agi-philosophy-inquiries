"""Refresh one live navigation block using the existing reversible layer."""
from pathlib import Path
import hashlib
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'research_cognition_physics'
STAGE=BASE/'archive_1116_'
sys.path.insert(0,str(ROOT/'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes
from run_research_motion_active import MotionLayout

def sha(raw): return hashlib.sha256(raw).hexdigest()

def main():
    lp=BASE/'_migration/motion_active_1106/layer.json'
    old_layer=lp.read_bytes(); layer=json.loads(old_layer)
    ledger=json.loads((STAGE/'_shared/dependency_ledger.json').read_text(encoding='utf-8'))
    done=ledger['completed_rounds']; active=ledger['in_progress_rounds']
    lines=['## 当前目标（1116起）：认知动力学与质能关系', '',
      '**目标active。** 按用户授权，从[1116起阶段](archive_1116_/README.md)审查PT来源，补上实际惯性及完整能动量连接，完成明确适用范围的狭义相对论基础推导。'
      '继承四原则、CO及SR；有理据的扩展可公开采用，未找到推导不等于已证明不可推导。'
      '不把质能公式直接写成认知公理，不把正质量部门变成全部物质的要求。', '',
      '上一阶段的[有限有效域共同实现](archive_1106_/_shared/stage_completion.md)保持结项。'
      '本阶段[目标与顺序](archive_1116_/_shared/goal.md)、[当前索引](archive_1116_/文件索引.md)及'
      '[依赖账](archive_1116_/_shared/dependency_ledger.json)给最新状态；以下旧“当前／完成／未启动”属于其原研究时点。', '']
    if ledger['new_adopted_axioms']:
        ready=ledger['goal_status']=='complete'
        lines=['## 当前目标（1116起）：认知动力学与质能关系', '',
          ('**1116—1119阶段已完成，目标complete。** ' if ready else '**四轮研究、公理采用和论文已完成，目标正在最终结项。** ')+
          '正式采用一组认知动力学公理D-P（四款）；F＋U＋C＋P、CO及SR继续共同使用，SR-S仅明示扩展，G未采用。', '',
          '[阶段论文](认知动力学公理与狭义相对论基础动力学_阶段论文.md)给出稳定正能有限静止材料部门内的质能关系、真实速度—动量、低速惯性和完整反作用。'
          '同一有限有效实验域的实际钟、原端点、物质、反冲及全部允许后续已有共同仿射证书。'
          '本次签收明确范围的狭义相对论基础动力学；绝对质量、材料谱、具体作用及更广物质覆盖仍有自由，不宣称经验真实性或全部物理完成。', '',
          '[正式采用与依赖](archive_1116_/_shared/adoption_and_dependencies.md) · '
          '[逐项目标验收](archive_1116_/_shared/stage_completion.md) · '
          '[阶段目录](archive_1116_/README.md)。'
          '本阶段4份报告、22组构造复算，科学校准增量0、累计3860；不重复旧器件，不启动下一目标。', '',
          '以下原阶段“当前／候选／未采用／未完成”保留其原研究时点，不覆盖本栏。', '']
    if done:
        items='、'.join(f'[{n}](archive_1116_/research_note_{n}.md)' for n in done)
        lines.append('**已完成并验收：'+items+'。** 标准代数及数值校准不计经验实验，历史科学计数保持3860。')
    if active:
        lines.append('正在推进：'+ '、'.join(map(str,active))+'；未完成轮次不计为已验收成果。')
    lines.extend(['', '---', ''])
    block='\n'.join(lines)
    updates={}; before={}
    for e in layer['entries']:
        path=ROOT/e['destination']
        if path.parent!=BASE or path.name not in {'README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md'}: continue
        raw=path.read_bytes(); before[path]=raw
        baseline=relocation_original_bytes(raw,e)
        txt=raw.decode('utf-8-sig'); nl='\r\n' if '\r\n' in txt else '\n'
        normalized=txt.replace('\r\n','\n')
        marker='## 当前目标（1116起）：认知动力学与质能关系'
        if marker in normalized:
            pos=normalized.index(marker); end=normalized.index('\n---\n',pos)+len('\n---\n')
            changed=normalized[:pos]+block+normalized[end:]
        else:
            split=normalized.index('\n')+1
            changed=normalized[:split]+'\n'+block+normalized[split:]
        old_marker='## 1106起阶段结项：有限有效域的实际狭义相对论共同连接'
        assert normalized.split(old_marker,1)[1]==changed.split(old_marker,1)[1]
        new=changed.replace('\n',nl).encode('utf-8')
        e.update(current_sha256=sha(new),current_bytes=len(new),
            edits=text_edits(baseline.decode('utf-8-sig'),new.decode('utf-8')))
        assert relocation_original_bytes(new,e)==baseline
        updates[path]=new
    layer.update(scope='One live navigation layer: completed 1106-1115 and dynamics 1116-1119',
                 active_stage='archive_1116_',goal_status=ledger['goal_status'],
                 new_formal_reports=9+len(done),latest_completed_round=max([1115]+done))
    for p,b in before.items(): assert p.read_bytes()==b,('Concurrent edit',str(p))
    assert lp.read_bytes()==old_layer,'Concurrent layer edit'
    for p,b in updates.items(): p.write_bytes(b)
    lp.write_text(json.dumps(layer,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    result=MotionLayout().verify_motion_layer()
    result.update(completed_rounds=done,in_progress_rounds=active,old_navigation_tails_unchanged=True)
    (STAGE/'_shared/navigation_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
