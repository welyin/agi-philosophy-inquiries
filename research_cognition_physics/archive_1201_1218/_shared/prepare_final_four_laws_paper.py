"""Build the final paper from the frozen, reviewed pre-adoption manuscript."""
from pathlib import Path
import hashlib
import os
import re
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SOURCE = HERE / 'four_laws_stage_paper_draft_after1218.md'
OUT = BASE / '08_认知热学角色与热力学四律_阶段论文.md'
SOURCE_SHA = '06cc5177ef656ae2899c4ea5edc5d220b1e14a7147680d1151c13af311e4e486'


def build():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    text = SOURCE.read_text(encoding='utf-8')
    replacements = [
        ('# 共同参考、完整资源与认知操作热力学四律', '# 认知热学角色与热力学四律'),
        ('## ——1201—1218阶段论文预审稿', '## ——共同参考、完整资源与有限操作的统一推导'),
        ('2026-10-15。本文综合已验收研究，不新增科学轮次。**热参考角色TR尚待采用；TW是可选充分工作实现，TC是资源范围登记。本文是待定稿论文，不是四律目标完成声明。** [现行目标](goal_four_laws.md) · [角色提案](cognitive_thermal_roles_proposal_after1218.md) · [原联合定理](four_laws_joint_theorem_after1218.md)',
         '2026-10-15。1201—1218阶段论文。**用户已正式采用分层热学定义：TR为共同热参考的部门要求；TW为可选充分工作实现；TC为冷却资源及定理适用域登记。** 本文在这一认知操作体系内给出四律的明确命题、证明、共同范围及有限过程证书。[正式采用](thermal_roles_adoption_after1218.md) · [逐项目标验收](four_laws_goal_completion_audit.md) · [原目标](goal_four_laws.md)。'),
        ('|新增角色候选|TR明确可扩展共同参考的资格；TW明确一种纯度中性工作支持方式。|认知动机本身不证明这些条件唯一，也不表示它们已正式采用。|',
         '|已采用的热学角色层|TR规定可扩展共同参考；TW定义可选充分工作子类；TC登记资源及适用域。|角色不普遍限制所有认知状态，不自动供应实际装置或操作权限；认知动机不证明其唯一性。|'),
        ('本文的采用状态与[角色提案](cognitive_thermal_roles_proposal_after1218.md)一致。提出“在热学部门采用TR”是一项理论定位决定；它不替代实际参考的准备与资格证明，也不自动完成四律结项。',
         '本阶段依据用户明确决定，正式采用[分层热学定义](thermal_roles_adoption_after1218.md)。它与四原则、CO、SR、D-P及CC共同构成本文的认知操作基础。理论中的角色、数学共同类的非空性、实际有限模型的资格和目标验收分别有证据；它们不能相互替代。原[候选提案](cognitive_thermal_roles_proposal_after1218.md)及预审稿保留历史状态。'),
        ('### 3.1 TR的操作化候选', '### 3.1 TR：可扩展共同热参考的操作定义'),
        ('本节未证明允许任意纯库存、任意反馈及无限活动资源的所有有限时间机器都不能冷却。目标终审必须明示采用哪一量词；新增角色被采纳也不会自动补齐更强量词。',
         '本阶段的动态第三律采用每份有限合格协议的精确不可达，并给固定资源预算族的统一余量；固定时间整族界在另有时间—资源控制时成立。允许任意纯库存、任意反馈及无限活动资源的普适有限时间版本不在已证范围中。所用操作类及q的独立来源进一步明确如下。'),
        ('|TR跨组合共同参考|扩大参与者后基准意义稳定，辅助代价在内部交代。|有选择力和删除检验的候选，尚待明确采用。|',
         '|TR跨组合共同参考|扩大参与者后基准意义稳定，辅助代价在内部交代。|已正式采用为理想热参考的部门要求；数学选择力与删除检验由1209—1210等承担。|'),
        ('阶段范围及逐律审计另由[第零／第一律审计](paper_audit_zeroth_first_after1218.md)、[第二／第三律及共同过程审计](paper_audit_second_third_after1218.md)承担。',
         '采用前逐律审计见[第零／第一律审计](paper_audit_zeroth_first_after1218.md)、[第二／第三律及共同过程审计](paper_audit_second_third_after1218.md)；采用后的目标逐项验收见[最终验收](four_laws_goal_completion_audit.md)，最终稿数学与范围复核见[论文终审](final_four_laws_paper_review.md)。'),
        ('## 10. 结论与尚未完成的验收', '## 10. 阶段结论、剩余自由与后继问题'),
        ('这份成果使剩余问题可以准确表述，而无需继续重复优化某个能梯或钟参数。当前不能宣布目标完成：TR的理论采用仍待明确；实际有限资格不代替全族权限；第三律的最终量词必须与目标对齐；论文与逐项终审还需独立核准。即使角色被采用，也仍须按这些证据完成验收，不能以名称替代证明。',
         '本阶段的结果是认知操作体系在所声明热学部门中的四律推导。第零律使用同对象、非标量能源及TR；第一律继承D-P的完整生成元身份；第二律使用完整熵账及所声明循环条件；第三律区分正温谱族的Nernst极限、有限合格协议的精确不可达与固定资源族的共同余量。采用后的逐项目标验收沿原目标核查这些范围及共同证据，不把一个有限实例当成全族证明。'),
        ('后续采用与定稿应保留全部历史版本。本文不改目标，不改变18轮科学计数，不自动采用新公理，不把普通热力学完成等同于广义相对论完成。',
         '具体Hamiltonian、材料参数、参考数值、正温分支、强耦合分割、纯库存和时间—资源增长仍按相应任务登记。由旧公理进一步生成或简化TR，扩大实际比较权限及冷却机器范围，属于后继基础问题。热力学—几何的同对象映射、Unruh、独立面积响应和有效Einstein动力学沿第9节接续。历史证明、提案与预审稿完整保留，科学计数仍为18轮、36组正式验收；普通热力学的阶段成果不等于广义相对论或经验真实性已获证明。'),
    ]
    for old, new in replacements:
        assert text.count(old) == 1, old[:100]
        text = text.replace(old, new)
    marker = '\n## 7. 工作角色、共同非空范围与条件的来源\n'
    corollary = r'''
### 6.4 推论：安全工作子类中的有限操作不可达

为使第三律的选择力不止于给定q的事件恒等式，考虑以下已证明存在的操作类。完整有限活动材料$X=S\otimes B$包括目标及本次参与的全部有限材料，联合输入$\Omega_X>0$；B可含热参考、库存或机器，其初始纯度和关联必须登记。与XR独立的工作支持按TW的完整振幅及全路径安全合同工作，保留其关联；末记录只使用1215已证的固定末QND copy、独立blank、无反馈且不后选合同。工作支持和这份记录辅助与材料X分开登记，其具体作用仍属于同一完整过程。如此得到的完整X通道$\Phi$保幺。设目标S的基态秩$g_S<d_S$。

此时同一实际末端非基态事件的效应为$E=\Phi^*(P_{\rm bad}^S\otimes I_B)$。由于$\Phi(I_X)=I_X$，事前已知

$$
q=\frac{\operatorname{Tr}E}{d_Sd_B}
=\frac{d_S-g_S}{d_S}>0,
\qquad
p_{\rm bad}\ge L_\Omega((d_S-g_S)d_B)>0.
$$

这是[1214的保幺冷却界](../1214/proof.md)、[1215的有限完整振幅实现](../1215/proof.md)、[1217的事件谱界](../1217/proof.md)及[1218的安全跨段闭包](../1218/proof.md)的直接合成推论。证明中的q来自可独立核验的振幅／路径规则，不是把“不能冷却”预写进TC。目标单体可以冷却且通道可以非保幺，受约束的是包含全部活动有限资源的X。任意被动R仍保留。

若实际整X通道相对该理想通道的half-diamond距离不超过$\varepsilon$（或有给出该界的全未知XR等距误差证书），令$q_0=(d_S-g_S)/d_S$并要求$\varepsilon<q_0$。对同一实际末事件、同一辅助准备及同一比较器，最大混合输入上有$q\ge q_0-\varepsilon>0$，可接定理七；固定$\Omega_X$上也可将理想谱尾下界减去相应事件误差。严格正值须由预算真正分辨。每份实际可在有限时间完成的合格协议都受此界；不要求所有协议有统一长度或统一时间成本。

本推论覆盖一个有非平凡能源交换的共同操作族，未把TW改成所有认知操作的普遍限制。开放可耗纯基态库存、输入相关辅助或新反馈时，应扩大活动域或重新认证通道；不得保留原q资格却删除新增资源。更一般的合格通道仍可直接应用定理七。
'''
    assert text.count(marker) == 1
    text = text.replace(marker, '\n' + corollary + marker)

    def relink(match):
        label, target = match.group(1), match.group(2)
        if re.match(r'^[a-zA-Z]+:', target) or target.startswith('#'):
            return match.group(0)
        name, sep, anchor = target.partition('#')
        resolved = (SOURCE.parent / unquote(name.strip('<>'))).resolve()
        path = Path(os.path.relpath(resolved, OUT.parent)).as_posix()
        return '[' + label + '](' + path + (sep + anchor if sep else '') + ')'

    text = re.sub(r'\[([^\]\n]*)\]\(([^)\n]+)\)', relink, text)
    return text


if __name__ == '__main__':
    assert not OUT.exists(), 'Final paper already exists; review changes explicitly'
    OUT.write_text(build(), encoding='utf-8')
    print(OUT)
