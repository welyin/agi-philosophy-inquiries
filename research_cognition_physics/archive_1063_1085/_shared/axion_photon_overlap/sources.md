# 来源、历史排重与采用边界

2026-10-09。对应[作者完成稿](../axion_photon_overlap_adoption.md)，已实际复算，作者提交时待独审及主线签收；成熟采用计0，不计新轮次。

## 直接读取的原始论文

|来源与版本|实际核对位置|本件复用|
|---|---|---|
|Matthew Reece，*Axion-Gauge Coupling Quantization with a Twist*，[2309.03939v1](https://arxiv.org/abs/2309.03939v1)，2023-09-07；[所读HTML](https://arxiv.org/html/2309.03939)页首为v1|式(2)—(6)，§4.2式(39)—(40)，§4.6及附录A|相关系数格、异常比归一、最小色系数时最近值8/3、重定义与适用边界。|
|Yichul Choi、Matthew Forslund、Ho Tat Lam、Shu-Heng Shao，*Quantization of Axion-Gauge Couplings and Non-Invertible Higher Symmetries*，[2309.03937v2](https://arxiv.org/abs/2309.03937v2)，2024-03-22；[所读HTML](https://arxiv.org/html/2309.03937)页首为v2|§2的G₆格及§3式(3.1)—(3.4)、脚注15|单QCD轴子、QCD主导质量假设，质量—光子比及0.15(1)历史结果。原文域壁简称在本件拆成明确K₃输入，不无条件推用。|
|Giovanni Grilli di Cortona、Edward Hardy、Javier Pardo Vega、Giovanni Villadoro，*The QCD axion, precisely*，[1511.02867v2](https://arxiv.org/abs/1511.02867v2)；[所读HTML](https://arxiv.org/html/1511.02867)|式(1)、(3)—(9)、(21)、(26)、§2.3式(37)、(43)—(46)|局部手征换基、质量归一、介子及NLO匹配、导数修正的低能抑制。固定采用1.92(4)与5.70(7)历史输入，不声称最新精度。|

直接阅读正文及公式，未看图、未提取曲线、未重做原文强子拟合或实验显著性。没有将论文的宇宙选择论证当作已证明的最小K₃来源。

## 直接复用与排重

1. [629](../../../archive_629_652/research_note_629.md)提供G₆全束的整数基、共同质量相位与测度字典。
2. [已签收全局周期采用§4](../../../archive_1046_1062/_shared/axion_global_period_adoption.md)已给$2K_3/3+k_\gamma\in\mathbb Z$、$E/N=2k_\gamma/K_3$，明确未据裸系数认证完整低能振幅。其[独审](../../../archive_1046_1062/_shared/axion_global_period_adoption_review.md)复算合法8/3例，但未签软振幅下界。
3. [局部QCD采用](../../../archive_1046_1062/_shared/qcd_axion_common_response_adoption.md)已给$m_a^2f_a^2=\chi$及EDM／局部态关系，没有光子匹配。原有效常数并不选择基本周期。
4. [π异常采用](../../../archive_1046_1062/_shared/pion_anomaly_adoption.md)只给π异常顶点与具体宽度匹配；不能单凭该采用给轴子全部NLO低能常数。
5. [上一整合](../../../archive_1046_1062/_shared/integration_branch_period_after1062.md)明确“完整光子幅度下界”未取得；[当前准入审计](../../_admission/generative_after1062.md)区分成熟接入与未完的生成义务。本件不会填平后者。

对全库正式报告及近期共享材料检索“QCD轴子／QCD axion”“axion-photon”“photophobic”“轴子—光子”“E/N”“1.92”；排除无关数字和普通分式命中后，未找到已签收式(2)—(5)这一共同接入。不是新物理定理：三个原论文已经承担关键物理结论，本文只把旧两份合同的交集和拒收对象写清楚。

## 复算口径

- **已生成固定`results.json`并实际运行复算通过。** 使用指定Python加`-B -X utf8`；默认入口输出PASS，运行前后四份作者资产SHA256不变。
- `check.py`只读，不提供保存或覆盖开关；默认入口重新计算并对比`results.json`。只用Python标准库，不联网、不安装依赖、不生成缓存。
- 有理数核G₆格、$K_3=1$最近点、误差界代数和$K_3=3$边界例。全整数格最小性的证明在正文，有限代入不替代该证明。
- 一份整数全局费米重定相核系数格保持，裸比值可改变；两种局部手征基的LO拆分核同一总括号的运输。两者均为已有匹配的记账复算，不重新推导强子有效作用，不声称导出1.92。
- 数字0.152仅用声明的$\alpha=1/137$、$m_af_a=0.00570\ \mathrm{GeV}^2$显示单位及中心数量级。0.04区间是条件式的代入，不是严格QCD误差证书；没有合并不独立的误差或把它当新置信区间。
- 不改任何已归档文件或根导航；正式1062／累计3838与各原标准保持。
