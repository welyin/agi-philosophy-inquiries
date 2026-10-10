# 绝对弱归一成熟采用：独立审阅

2026-10-09。审阅者：`next_structural_selection`，独立于作者。**结论：通过；无待修项。** 只签下列固定版本和范围，成熟采用计0，不增加科学／经验组或认知公理，不立1063。

## 1. 已核内容

已完整阅读[正文](electroweak_absolute_normalization_adoption.md)、[来源](electroweak_absolute_normalization/sources.md)、[代码](electroweak_absolute_normalization/check.py)与[结果](electroweak_absolute_normalization/results.json)，并实读以下一手原文：

- [MuLan 1010.0991v2](https://arxiv.org/html/1010.0991v2)：寿命提取的是指定修正下的μ衰变归一；不是普适性自身的证明。打印寿命及耦合、误差转录正确。
- [Awramik等 hep-ph/0311148v3](https://arxiv.org/html/hep-ph/0311148v3)：展开的$1+\Delta r$、参数化系数、更新的$G_\mu$、质量定义及误差范围核对一致。$\Delta q$不再重复加进匹配；正根支有明确物理选择。式(9)的参考值没有与表2遗留绝对中心混用。
- [CMS 2412.13872v2](https://arxiv.org/html/2412.13872v2)：直接质量读数及running-width约定一致。稿件保留了$M_Z$闭合／分辨率、产生模型、宽度关系和两电荷质量相等的采用条件；未将其与预测端强称完全独立。

0.25 MeV参数化精度与约4 MeV遗漏高阶估计性质不同，均未被升级成全阶硬界或自动独立Gaussian标准差。top的MC质量与匹配质量方案不能无误差等同，正文已保此义务。

## 2. 独立计算与只读复算

实际运行默认`check.py`，退出成功；运行前后四份作者文件SHA256一致。15个本地引用全部存在。

另用70位Decimal直接展开原式(9)，**没有导入作者函数**。保持其他原参考输入，仅取$M_H=125$ GeV，得到

$$
M_W=80.3656181142090731303414916289132\ldots\ {\rm GeV}.
$$

独立代数转换原$M_Z=91.1875$、$\Gamma_Z=2.4952$ GeV得到
$\bar M_Z=91.1533805818080055515979620378729\ldots$ GeV，与保存结果一致。手算式(5)五个导数及固定$\Delta r$时的$\partial M_W/\partial\log G_\mu$符号、单位均一致。

这些只是出版参数化／代数诊断。没有重新计算完整两圈作用、现实输入协方差或CMS似然，也没有新增现实$M_W$预测、残差或显著性。浮点／高精度校准没有替代原文适用域。

## 3. 净覆盖与范围

[1052](../research_note_1052.md)的轻子比值消去了共同绝对归一，[1061](../research_note_1061.md)以其作CKM抽取输入；本稿明确补出$\Delta q\to G_\mu\to\Delta r\to M_W$的同参数责任，具有成熟共同接口价值。

正文只采用指定可重整化SM分支。其余旧参考保持不动的125 GeV示范不是现实输入集合；原作者的CMS相容判断也没有被包装为本项目独立检验。未知高维系数、共享误差、完整仪器与认知生成均未免费取得。不要求为此继续搭建新的拟合或器件。

## 4. 签收快照

|作者资产（相对本目录）|SHA256|
|---|---|
|`electroweak_absolute_normalization_adoption.md`|`a6be16f19d031c44320cc17968647a73b1a1bdcd4d624ee45b74f3fcad3389c9`|
|`electroweak_absolute_normalization/sources.md`|`fd430ca31aee7a0d0a20f0f58e819e47879b28fdff5fa40ae3b114ea67f30a26`|
|`electroweak_absolute_normalization/check.py`|`1cd7ac49b1eec958605930d740de2650ae24dcff1187714b8d9bdb62131add1d`|
|`electroweak_absolute_normalization/results.json`|`0d67d1e9a8290fdd52ed7cdeffb07c2e20b8fb76f8c0dec99b452ade0e59124c`|

本审阅只新增当前文件；作者资产、历史和导航均未修改。整体ROADMAP未由此次采用完成。
