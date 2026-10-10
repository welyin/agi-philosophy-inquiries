# 氢谱共同匹配：来源与版本

2026-10-09。仅支撑相邻[采用正文](../hydrogen_spectral_matching_adoption.md)，成熟采用计0。没有下载大型实验档案或读取图像数字。

## 一手来源及实际读取范围

1. **Maisenbacher等，Nature 650, 845–851 (2026)，DOI 10.1038/s41586-026-10124-3。** [开放原刊](https://www.nature.com/articles/s41586-026-10124-3)。实读主文式(1)、(7)—(11)，Methods重心／HFS段及式(17)—(18)、相关误差段、Data/Code availability。式(17)给共同能差比；式(18)仅是指定领先S态修正的灵敏度。原刊采用CODATA2022能级，并用2024缪氢理论半径作输入；新2S–6P派生的半径、Rydberg常数和Lamb位移不另算独立资料。出版的0.53 kHz为标准不确定度口径，本件不赋予新分布／覆盖／显著性。
2. **Pachucki、Lensky、Hagelstein、Li Muli、Bacca、Pohl，Rev. Mod. Phys. 96, 015001 (2024)。** [作者原稿2212.13782v3](https://arxiv.org/html/2212.13782v3)，[期刊DOI](https://doi.org/10.1103/RevModPhys.96.015001)。实读Table I、§IV、§V及有限尺寸／两光子交换定义。表中μH中心量为$E_{\rm QED}=206.0344$ meV、$E_{\rm NS}=0.0289$ meV、半径项$-5.2259r_p^2$ meV、实验Lamb位移202.3706 meV，由此提取约0.84060 fm。代码仅核出版舍入中心式；未重算核结构积分或重新推出完整半径误差。
3. **Mohr等，CODATA recommended values 2022，Rev. Mod. Phys. 97, 025002 (2025)。** [2022调整作者原稿2409.03787v1](https://arxiv.org/html/2409.03787v1)，[NIST原刊PDF入口](https://physics.nist.gov/cuu/pdf/RevModPhys.97.025002.pdf)。实读HTML §III.1—2及§IV有关核结构和相关不确定度。共用的$n^{-3}$误差与按态区分的误差不应混同；核尺寸与极化等也不是仅一份无理论半径。该稿是固定历史理论／调整来源，本件未采用其全局调整输出作另一个独立谱实验。

## 本次可复算范围

- `check.py`以有理数核共同能级比的灵敏度、有限比值恒等式、严格有界分母余项及秩1相关误差；这些是代数查错，不是物理余项证明。
- 以Decimal核原刊公开误差分组归账；原始分组及不相关近似作为采用输入。0.00 kHz出版残差没有被转为p值。
- 以μH表内舍入中心核“半径为理论提取”的代数；不输出新半径估计或新不确定度。
- 未取得完整实验分析代码。公开数据登记[Zenodo 15874385](https://zenodo.org/records/15874385)在原准入已核，本次不下载62 MB主包，不据摘要声称重建原始拟合。

## 项目复用与量词

[原准入](../../_admission/hydrogen_overlap_after1062/selection.md)和[完整输入综合](../overall_evidence_after1062.md)提供任务与范围；[1059](../../research_note_1059.md)保共享标定的独立性边界；[1049](../../research_note_1049.md)保同一理论变化须同步运输标定与来源的原则。此处没有将旧材料模型或原子反冲结果改成氢谱预测，也没有替它们修正参数。
