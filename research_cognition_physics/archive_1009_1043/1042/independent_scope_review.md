# 1042 独立操作范围、历史去重与准入审阅

审阅者：独立物种／操作范围代理 `/root/round1036_species_selection`。日期：2026-10-08。

**状态：冻结后的独立操作范围、历史去重与准入审阅通过。** 本文件是唯一新增审阅资产；未修改作者报告、代码、结果、导航或历史文件。审阅者未参与本轮代码与科学报告的撰写；此前只向作者提供准入及范围意见。本审阅不冒称另一份一般 Harlow 定理证明。

## 1. 阅读范围与结论

已读 [正式报告](../research_note_1042.md)、[科学代码](complementary_entropy.py)、[结果](complementary_entropy_results.json)、[准入](selection_audit.md)、[依赖增量](input_dependency_update.md) 与 [后继](NEXT.md)。成熟来源直接核对 [Harlow 原文定理 1.1、§5](https://arxiv.org/html/1607.03901)。

本轮可准入一项新的共同对象接口：给定全码态上的互补恢复合同，使两侧熵、相对熵与同一中心余项相容；随后在恢复代数、代表、解码与各编码内逻辑相对熵均保持的族中，展示中心余项数值仍可改变。一般成熟定理、随机态次数、既有中心熵分解不算本轮原创增量。

未见阻碍验收的科学或范围错误。冻结前，作者已补入两项历史明示：391 的零中心余项互补特例，以及 839 的旧精确相对熵运输；已回读确认。它们不推翻上述增量，但使“首次”边界更准确。

## 2. 操作合同与物理范围

逐项检查结论如下。

1. **全码恢复是额外合同。** 报告明确有限复量子结构、等距编码、A/B 划分、全码输入、恢复权限与熵迹约定。A4/A7 只提供动机，未被写成自动保证这份全称恢复合同。代表满足交织及伴随关系；不是只证明压缩矩阵元相等。
2. **中心不是基本超选。** 码空间容许跨扇区相干，等距编码保留全输入及被动参考。局部迹因另一侧的正交标签而去掉相干，不是整体禁止相干。解码只输出代数限制态，未声称单侧恢复完整 8 维未知态。
3. **同编码与跨编码分开。** 相对熵式比较固定同一个编码内的两逻辑态；变更辅助参数后的物理相对熵没有被宣称保持。报告还给固定辅助效果的概率变化，证明全部物理仪器通道并不相同。
4. **模算符变分固定对象。** 投影式及有限相对熵余项固定编码、参考及代数迹。参考的扇区限制须有正支持；改变编码另记编码、中心算符和资源的变化。未将有限模算符直接称作物理能量或几何 boost。
5. **资源确实入账。** 辅助纠缠是准备输入。另声明辅助 Hamiltonian 后，不同参数改变能量期望；该 Hamiltonian 也没有被说成由代码生成。反例范围是明确较弱合同，未认证完整 A1_D—A7 共同物理实现。
6. **固定编码的连续界有条件。** 有限维数、给定到固定编码的输出距离及其任务域是前提；近似恢复本身没有被替换成这份距离前提。未以有限精度认证精确可恢复代数，保留 391 的噪声边界。
7. **没有几何越界。** 无量纲中心熵算符尚不是物理面积。单一二分区不产生区域族、极值曲面、面积测量、Newton 系数、几何模流或 Jacobson 平衡。此条件工具不宣布整体生成目标完成。

## 3. 历史去重与实际增量

- [391](../../archive_370_428/research_note_391.md) 已从通道求最大可纠正代数与共同尖锐中心。其原 CZ 等距在 AB|CD 切分可读作互补特例：父侧为 M₂⊕M₂，读者联合为其交换子 C²；分支内辅助为乘积态，中心熵余项为零。故“首次有互补实例”不成立。1042 新在非零可变辅助谱下统一验收两侧熵／模／相对熵及剩余自由。
- [305](../../archive_301_341/research_note_305.md) 已有固定参考熵第一定律和精确相对熵余项；[306](../../archive_301_341/research_note_306.md) 与 [312](../../archive_301_341/research_note_312.md) 的几何、参考和模流条件仍独立。1042 不重复认领它们。
- [636 §4](../../archive_629_652/research_note_636.md) 已给实际 Gauss 切分的标签 Shannon 熵、表示指标熵和迹约定；也已有同电来源一二阶矩而不同熵。1042 不将一般中心熵分解或少数摘要不足定熵再算一次。新反例保留的是全部所列逻辑通道和各编码内的任意逻辑态对相对熵，范围严格强于仅保低阶来源矩。
- [837](../../archive_819_853/research_note_837.md) 已给区域可读见证及相对熵下界；[839 §7](../../archive_819_853/research_note_839.md) 对共同独立余部准备已给记录内容在实际自由区域的相对熵精确等式。故“一族相对熵可精确运输”也不是首次。1042 使用的是全码、双侧互补恢复与可变中心资源的一致接口，未将 839 的真实区域物理字典免费移植到此有限码。

以上为准入理由及限制，不把本轮有限模型等同现有完整物理父对象。

## 4. 已实际运行的独立核算

未导入作者模块，也未复用其大矩阵编码。下列脚本以有理数核两扇区辅助概率、半迹距离及能量；以形式独立振幅 c_(α,k) 核每个基矢的矩阵单位交织，因而验证与参数无关，而非仅对三个浮点样本检查。

执行环境为既有 Python，参数 `-B -X utf8`。代码只输出，不写任何文件。

~~~python
from fractions import Fraction as F
from math import log
import json

ps = ((F(1, 2), F(1, 4)), (F(1, 4), F(1, 2)))
rows = []
for alpha in range(2):
    p, q = ps[0][alpha], ps[1][alpha]
    distance = (abs(p-q) + abs((1-p)-(1-q))) / 2
    ep, eq = 2*(1-p), 2*(1-q)
    assert distance == F(1, 4)
    assert abs(eq-ep) == F(1, 2)
    rows.append([alpha, str(p), str(q), str(distance), str(ep), str(eq)])

def encoded_term(alpha, a, b, k):
    # Last pair labels the same formal amplitude c_(alpha,k).
    return (alpha, a, k, alpha, b, k, alpha, k)

def physical_unit(side, sector, i, j, term):
    t = list(term)
    label, logical = (0, 1) if side == 'A' else (3, 4)
    if t[label] != sector or t[logical] != j:
        return None
    t[logical] = i
    return tuple(t)

cases = 0
for side in ('A', 'B'):
    for sector in range(2):
        for i in range(2):
            for j in range(2):
                for alpha in range(2):
                    for a in range(2):
                        for b in range(2):
                            for k in range(2):
                                left = physical_unit(side, sector, i, j,
                                    encoded_term(alpha, a, b, k))
                                logical = a if side == 'A' else b
                                right = None if alpha != sector or logical != j else (
                                    encoded_term(alpha, i if side == 'A' else a,
                                                 i if side == 'B' else b, k))
                                assert left == right
                                cases += 1
assert cases == 256
# Both (i,j) and (j,i) are included, so adjoints are covered.
h = lambda p: -float(p)*log(float(p))-(1-float(p))*log(1-float(p))
gap = h(F(1, 2))-h(F(1, 4))
assert 3**3 > 2**4
assert abs(gap-(3*log(3)/4-log(2))) < 1e-15
print(json.dumps(dict(passed=True, rows=rows,
    matrix_unit_representatives=16, coefficient_cases=cases,
    entropy_gap=gap, cross_encoding_relative_entropy=log(4/3)/2)))
~~~

实际结果：

| 扇区 | 辅助零结果概率 p/q | 半迹距离 | 总辅助能量 p/q，单位 Δ |
|---|---|---|---|
| 0 | 1/2，1/4 | 1/4 | 1，3/2 |
| 1 | 1/4，1/2 | 1/4 | 3/2，1 |

两侧共 16 个矩阵单位、256 个形式振幅系数身份全部通过。扇区 0 的熵差为 0.130812035941137，跨编码辅助相对熵为 0.14384103622589042。熵差严格为正的解析检查等价于 27>16；浮点值只用于显示。能量的有符号差随扇区反向，绝对差均为 Δ/2，与报告指定扇区 0 的见证一致。

## 5. 冻结后只读签核

作者首份冻结通知后，已运行 [verify_round1042.py](verify_round1042.py) 默认只读入口；该入口包含科学代码的默认只读复算，返回 passed=true、8 份作者资产、9 份历史输入、9 个作者报告／导航本地链接。另用独立标准库脚本逐项重新计算收据内的 8+9 份 SHA-256，全部一致；本审阅的本地链接亦逐项存在。没有再次生成或改写作者结果。

签核所对应的 SHA-256：

- 作者收据 research_round_1042_checks.json：295b7dbeeeca06c76adcc4b0f16a730f38d842e789cedf8baa7ec92e541d7ef6。
- 正式报告 research_note_1042.md：5e538633726ed3dac6c756d70b3d2887c7b213cf9a56ca04c5dbbe40af89ff5d。
- 科学代码 complementary_entropy.py：ab4e1e720470480c0bcc6721b8bda8ebf1430d269f9073968f111a0566c1e70b。

冻结前 Kraus 秩记录从 tuple 改为 JSON 可原样比较的 list，未改变计算值或数学对象；默认复算已验证结果一致。本独立审阅仅在本文件签收，不修改作者冻结收据中的独立审阅标志，也不自行增加科学累计。

**最终结论：通过。** 本轮成立的是明确额外操作合同下的结构派生与剩余数值自由；其范围、历史继承、资源差异和未完成物理字典均已明示。无待解决的实质修订意见。
