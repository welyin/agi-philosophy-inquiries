"""Maintain one reversible live-navigation layer for the active field stage."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

BASE = ROOT / 'research_cognition_physics'
LAYER = BASE / '_migration/field_active_1120/layer.json'
NAMES = ('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')
NOTICE = '''## 当前目标（1120起）：局域相对论有效场与共同能动量来源

**目标状态active，研究目录为[archive_1120_](archive_1120_/README.md)，最新正式轮次1152。** 继承F＋U＋C＋P、CO、SR及D-P，保原任务和材料范围；SR-S仅明示扩展，G未采用。前阶段1106—1119及五篇编号论文保持结项。

本阶段连接原位置与实际钟、区域操作和局域演化、因果传播、局部能动量与完整反作用，为旧规范／引力条件定理提供同一对象。允许明示有限有效材料与误差；不把全局守恒、表示或边界散射直接当成局域场与应力。

[阶段目标](archive_1120_/_shared/goal.md) · [依赖与复用账](archive_1120_/_shared/dependency_ledger.json) · [轮次索引](archive_1120_/文件索引.md)。1120—1125给必要接口和限定障碍；1126—1129连接正常连续交接与自治来源；1130接原钟校准，1131筛选两准备的梯度方案。[1132](archive_1120_/research_note_1132.md)使正常有限质量支撑承担真实反冲，并装入完整旧D；[1133](archive_1120_/research_note_1133.md)补同源自治接收、两项不对易记录及继续，限明示区域材料与固定g／P读口。

**[1134](archive_1120_/research_note_1134.md)** 以原正常连续过程的静止纤维定义正内部质量，构造真正相互作用Poincaré表示；同一空间平移和旋转精确保持。原有限矩准备的完整输出差小于3.000004×10^-8，未知交接小于0.001258186，因果有限菜单小于0.000120170。新H的分项资源有限性与相反总冲量独立证明；局部应力或1133定向接收器尚未并入。

**[1135](archive_1120_/research_note_1135.md)** 补同一实际Y位置与原裸旧钟的非零boost联合端点。原288维D、正常源、LLP纠缠、内部角动量及空间相位共同保留；联合概率和完整平滑仪器后态误差均小于3.126437×10^-8，错误实验室时间字典仍有大于0.07453608的读数分离。正常运动预部署和终末读取是明示材料权限；末读装置的自治做功与四P闭合未证明。

**[1136](archive_1120_/research_note_1136.md)已通过独审与主线验收：内部密度形变有同一正常控制器的有限可读响应。** 按[事前合同](archive_1120_/_shared/after1135_dilation_admission.md)保完整来源、场余部和Q/R；实际余弦读数差大于0.01215034，对应效果概率差大于0.00607517。真实内部力由完整H推得，资源域与空间反冲另证。未形变支路保旧交接，活跃支路为明示新控制任务；终末仪器许可和测量做功边界保留。新增公理0、科学校准0，历史累计3860。

**[1137](archive_1120_/research_note_1137.md)已验收：正常方向资料、两项实际接收仪器与继续，共同接入同一相互作用四P。** 全Q/R接收误差小于0.001446465，正确与错误方向的同裸P概率差大于0.49421414；非零boost下原裸P／实际Y完整终末仪器误差小于0.001447。转子资源、真实角反作用、空间反冲及后态域另证。限b=0、明示正常准备和终末读取，未领取T=4旧Mr或探测器自治做功。本阶段18轮共72组复算，新增公理0、科学校准0，累计3860。

**[1138](archive_1120_/research_note_1138.md)已验收：同一实际H上的有限分辨率区域能源三点证书。** 完整Fock区域图界、各比较过程独立资源界及有效J弱读通过。三点误差小于1.8×10^67，同一壳流均值大于2.89×10^68；固定读口完整扰动小于2.233×10^-7，原接收总误差仍小于0.001447。定量读口只J；端点E测读、实际H连续历史流及自治探测器未签收。本阶段19轮、77组复算，新增公理0、科学校准0，累计3860。

**[1139](archive_1120_/research_note_1139.md)已验收：原区域接收器的有限局域细化存在匹配下界。** 保原时限，再加入迟发远翼干预及小区域实际读口时，任何严格因果父过程的单项概率匹配误差必大于0.014605，超过1%目标。原正常海、Y尾和插入资源分别控制；新增读写的自治实现未认领，原粗区域合同没有被反证。本阶段20轮、80组复算，新增公理0、科学校准0。

**[1140](archive_1120_/research_note_1140.md)已验收：原g两种记录与继续可接到有限相干汇集边界桥。** 原模式误差小于1e−5，完整CP操作差小于6.286e−5，保正常CAR及任意Q/R。新增时限和全盒部署明列；逐胞元保持与旧保持的生成元差至少1，不能用操作误差替代完整背景匹配。本阶段21轮、83组复算，新增公理0、科学校准0。实际Y/钟、自治来源及原完整动力学局域化未签收。

**[1141](archive_1120_/research_note_1141.md)已验收：原H_o交互区的相对源能、独立交换功及相对波流分别受控。** 保原参数、完整Q/R、正常Y/O及全部接收量子尾；源能差和功率期望绝对值积分均小于10^23，相对流期望绝对值积分小于10^24。独立能源形式、退迟力与时间弱收支共同核验；净增量来自共同过去相消，未要求简并交接必须正净功。本阶段22轮、87组复算，新增公理0、科学校准0。实际BT迁移、未分辨源能与材料的完整局部账、能源仪器仍开放。

**[1142](archive_1120_/research_note_1142.md)已验收：原固定交互区的源能、独立波流和交互功共同接入实际BT动力学。** 原正常准备、量子尾及参数保持；弱缺陷绝对期望积分小于2×10^11，同一裸功率期望绝对值积分小于10^23。实际低涨落、二阶动量域、完整Fock图界及半群弱身份分别证明，两份资格审查与[正式独审](archive_1120_/1142/independent_review.md)通过。本阶段23轮、91组复算，新增公理0、科学校准0。实际源能/流的旧单项精度、完整材料局部账与能源仪器未自动获得。

**[1143](archive_1120_/research_note_1143.md)已验收：活跃区域响应由同一正常来源、正程序材料和实际四P共同实现。** 固定裸末读概率差大于0.0335，两个程序各自完整输出误差小于0.0015；控制器、源及Y/O反作用保留。两种部署是明确不同的新增物理操作，没有因此选出唯一局部应力。本阶段24轮、96组正式复算，新公理0、科学校准0。

**[1144](archive_1120_/research_note_1144.md)已验收：独立材料记录的量子中介及有限误差方向判据。** 正常自治单元全输出误差小于0.000709；连续有质量场固定剖面精确回收有障碍，有限期限记录误差仍小于0.000865。两实例尚未共同连接，实际位置/钟、完整局部来源及四P保持开放。冻结资格稿能源绘景错误已在[正式证明](archive_1120_/1144/proof.md)勘误；7组复算和两份交叉审阅通过。阶段25轮、103组正式复算，新增公理0、科学校准0。

**[1145](archive_1120_/research_note_1145.md)已验收：同一外给记录源的共同运动运输及强度字典筛查。** 保完整源相位、Q/R、正常场和KG四冲量；错误额外归一造成同裸指针概率差大于0.59533，即使末态场能差小于5.445×10^-20。4组复算与独审通过，实际自治材料、钟位置和全装置四P仍开放。阶段26轮、107组正式检查，新公理0、科学校准0。

**[1146](archive_1120_/research_note_1146.md)已验收：同一正常运动材料的实际钟、标量操作及量子来源共同连接。** 保全CAR、正常来源后态和任意Q/R；完整输出误差小于7.38×10^-6，实际钟字典差大于2.4×10^-5，源型概率差大于0.0934。5组复算及独审通过，能源域与真实源力另证。原CO/288维钟、独立两站记录、全装置四P和点态Ward仍开放。阶段27轮、112组正式检查，新公理0、科学校准0。

[1146后联合资格](archive_1120_/_shared/after1146_joint_screen_report.md)与实际材料筛选已复用到1147；普通域、衣着及二阶工具不重复计数。

**[1147](archive_1120_/research_note_1147.md)已验收：同一正常束缚Dirac材料的实际钟、独立记录与初始短窗整站因果。** 保全部海、正常中介及未知Q/R；完整记录误差小于0.0025（含明示末读容差），同一效果概率差大于0.995，实际能源及反作用另证。4组复算及[独审](archive_1120_/1147/independent_review.md)通过。固定壁、区域耦合、末读权限仍为输入，谱常数未数值认证为器件参数。阶段28轮、116组正式检查，新公理0、科学校准0。

**[1148](archive_1120_/research_note_1148.md)已验收：同一正常束缚材料的实际钟、记录、正常支撑与整体Y共同接入非零Lorentz运动。** 完整仪器误差小于0.00252，记录差大于0.9949，实际钟与错误计时字典分离大于0.035；保全海、未知自旋及Q/R。3组复算、[独审](archive_1120_/1148/independent_review.md)和[专项审阅](archive_1120_/_shared/after1147_spectral_boost_review.md)通过。本阶段29轮、119组正式检查，新公理0、科学校准0。

**[1149](archive_1120_/research_note_1149.md)已验收：原实际材料的局部密度共同匹配二阶记录与外区静态反作用。** 一个正耦合幅度使全部不重叠间距上的跨系数、外场与相反广义力一致；仅保积分荷的约4.25%差只作数值诊断。3组复算及[独审](archive_1120_/1149/independent_review.md)通过。阶段30轮、122组正式检查，新公理0、科学校准0。

**[1150](archive_1120_/research_note_1150.md)已验收：原长任务的实际中途干预与整站动态能源收支。** 保原耦合、时限、全正常材料/场及未知自旋/R，实际非对易中插后完整历史差小于2.21×10^-4；场不重置。同一H给关联态上的整站短窗因果及含实际功、交互存能和独立流的区域弱收支。3组复算及[内部交叉审阅](archive_1120_/1150/independent_review.md)通过，已披露审阅者参与前置资格。阶段31轮、125组正式检查，新公理0、科学校准0。

**[1151](archive_1120_/research_note_1151.md)已验收：实际材料操作的共同运动与有界单调用检验。** 同一正常材料及相互作用四动量下，新增有限相干平移测试使裸U与正确运输V的概率差大于0.499897，保未知A味、spin及R；原参数不变，资源域另证。3组复算及两份[内部审阅](archive_1120_/1151/independent_review.md)、[范围审阅](archive_1120_/1151/scope_review.md)通过，披露前置参与。阶段32轮、128组正式检查，新公理0、科学校准0。

**[1152](archive_1120_/research_note_1152.md)已验收：运动操作的空间控制来源限制。** 保1151同一实际正常对象与原参数，仅增平移平凡内部钟/程序的实现类有大于0.499897的完整误差下界；独立来源达到0.001误差需相对平移态距大于0.99799997。源可有反冲及后态变化，共享支撑/相关准备未被普遍排除。3组复算和两份[内部审阅](archive_1120_/1152/independent_review.md)、[范围审阅](archive_1120_/1152/scope_review.md)通过，均披露前置参与。阶段33轮、131组正式检查，新公理0、科学校准0。

[下一步](archive_1120_/1152/NEXT.md)先固定实际空间参与者、相对准备和旧系统嵌入，再构造共同局域过程。停止纯内部计时/电池参数优化；必要来源条件不等于自治实现、共同局部应力或旧规范/引力前提已完成。1153未分配，阶段目标active。

停止质量、波包或读口参数优化，不把修好一个候选当成原始目标。保全部历史，不预设微观连续或最小尺度，不以全部QFT、标准模型、GR或全UV完成为前置门槛。历史复算使用[当前只读入口](_migration/field_active_1120/README.md)。以下旧“当前／下一步／不启动”保留原研究时点，由本栏说明当前状态。

---

'''


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    frozen_path = BASE / '_migration/merge_1106_1119/manifest.json'
    frozen = json.loads(frozen_path.read_text(encoding='utf-8'))
    prior = {e['destination']: e for e in frozen['entries']}
    existing_bytes = LAYER.read_bytes() if LAYER.exists() else None
    existing = {e['original']: e for e in json.loads(existing_bytes)['entries']} if existing_bytes else {}
    entries, changes = [], []
    for name in NAMES:
        path = BASE / name
        relative = path.relative_to(ROOT).as_posix()
        current = path.read_bytes()
        original = relocation_original_bytes(current, existing[relative]) if existing else current
        assert sha(original) == prior[relative]['current_sha256']
        assert len(original) == prior[relative]['current_bytes']
        encoding = 'utf-8-sig' if original.startswith(b'\xef\xbb\xbf') else 'utf-8'
        before = original.decode(encoding)
        newline = '\r\n' if '\r\n' in before else '\n'
        split = before.index('\n') + 1
        after = before[:split] + newline + NOTICE.replace('\n', newline) + before[split:]
        new = after.encode(encoding)
        entry = dict(original=relative, destination=relative,
                     original_sha256=sha(original), current_sha256=sha(new),
                     original_bytes=len(original), current_bytes=len(new),
                     edits=text_edits(before, after))
        assert relocation_original_bytes(new, entry) == original
        entries.append(entry)
        changes.append((path, current, new))
    support_names = ['scripts/run_research_merged.py', frozen_path.relative_to(ROOT).as_posix()]
    layer = dict(schema='field_active_navigation_v1', date='2026-10-12', entries=entries,
                 unchanged_support_sha256={name: sha((ROOT / name).read_bytes()) for name in support_names})
    if args.apply:
        assert (LAYER.read_bytes() if LAYER.exists() else None) == existing_bytes
        for path, before, after in changes:
            assert path.read_bytes() == before
        LAYER.parent.mkdir(parents=True, exist_ok=True)
        LAYER.write_text(json.dumps(layer, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        for path, before, after in changes:
            if before != after:
                path.write_bytes(after)
    print(json.dumps(dict(status='APPLIED' if args.apply else 'CHECKED', files=4,
                          frozen_merge_unchanged=True, original_navigation_recoverable=True)))


if __name__ == '__main__':
    main()
