# 致密物质—潮汐采用：一手来源与核验位置

2026-10-09。为[正文](../neutron_star_tidal_adoption.md)提供定位。本表中的式号属于各自原文，不混用$\lambda$、$\Lambda$及两种加权记法。未读取／重算原始GW时间序列，未运行数值EOS。

|一手来源与读取版本|本次核过的位置|采用边界|
|---|---|---|
|[Hinderer，0711.2420v4，2009-03-07](https://arxiv.org/html/0711.2420v4)|页首Erratum；§III式(4)—(5)；§IV式(15)及更正式(23)|原式(20)、(23)和部分表值有更正；使用v4。采用静态偶宇称四极与正确k2，不引用旧表作新校准。|
|[Hinderer、Lackey、Lang、Read，0911.3535，PRD81,123016](https://arxiv.org/html/0911.3535)|§II式(5)—(15)：TOV、响应、界面；§IV式(20)—(22)；附录A|$f=d\epsilon/dp$是该慢变闭合。表面非零能密度须匹配跳跃。$\widetilde\lambda$有量纲，正文的$\widetilde\Lambda=32\widetilde\lambda/M^5$。本次未复制EOS数值表或其旧预报。|
|[Flanagan—Hinderer，0709.1915，PRD77,021502](https://arxiv.org/html/0709.1915)|式(5)—(10)，式(12)的Fourier相位约定及其后加权定义|轨道与总四极辐射共同进入领先相位；慢变、微小形变范围。原特定模型的频率误差不推广为统一界。|
|[Alford—Han—Prakash，1302.4732，PRD88,083013](https://arxiv.org/html/1302.4732)|§II—III的混合星稳定支分类|只用来定位分离稳定支的边界；不采用其参数模型作本项目EOS或拟合。|
|[LIGO/Virgo，1805.11581v2，2018-10-15](https://arxiv.org/html/1805.11581v2)|§II.1低自旋／资料／频段；II.2 PhenomPNRT；II.3关系误差；II.4谱EOS；III公布区间；IV半径下界的解释|同EOS、EOS族和脉冲星质量要求属于条件。90%后验非硬界；同事件数据和派生量不可独立累乘。固定论文版本，不是最新全球复合分析。|

## 本地检查做了什么

[check.py](check.py)用精确有理数检查两质量权重、两种$\widetilde\lambda/\widetilde\Lambda$约定、领先相位系数、质量整体缩放、交换及等质量极限，并用量纲指数检查SI换算。测试质量只是代数代入，不是EOS解或拟合参数。默认运行只读比较[results.json](results.json)，首存仅允许不存在的结果文件。

上述检查不能确认TOV数值积分、星体稳定性、近似余项或公开后验；那些部分分别是成熟理论采用与原论文条件性推断。正文和[准入审计](../../_admission/after1060_dense_matter/selection.md)的范围独立保留。
