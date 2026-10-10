# 来源、访问状态与复用边界

2026-10-09。对应[采用正文](../page_geilker_record_source_adoption.md)。固定1981原实验及本次实际可核的争论范围，不称最新实验汇总。

## 一、实验一手来源：全文已读

Don N. Page and C. D. Geilker, *Indirect Evidence for Quantum Gravity*, **Physical Review Letters 47, 979–982 (1981)**，收到1981-06-09，发表1981-10-05。

- [原刊DOI与摘要](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.47.979)。
- [公开四页原刊PDF](https://bpfoundations.org/plans/page1981.pdf)：本次通过文本提取完整阅读，无图像取数。
- APS harvest入口返回401，未尝试登录或绕过访问限制。

|位置|本次采用内容|
|---|---|
|p.979，式(1)—(2)|无条件期望应力方程及测量解释／守恒讨论；不将任意坍缩直接代入的困难扩大为所有随机经典几何不可能。|
|p.980，实验段|Co60双计数、决策分组、实验者接续、来源与扭秤几何。原计数均值／标准差为1509.1±31.0、887.6±23.0；这些不独立认证每轮Born权重恰好1/2。|
|p.981，程序与数据段|AB／BA各含两个30分钟段，预定开始时刻；710秒振荡、极值拟合、十个已发表平衡差。程序只录该十数，不重建原时序。|
|p.982，统计与结论|原相关密度及粗G核查：$G=(6.1\pm0.4)\times10^{-8}\,\mathrm{cm^3\,g^{-1}\,s^{-2}}$，不确定度为均值的一个标准差。换成SI是$(6.1\pm0.4)\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}$。没有独立重算这个G估计；它与同一扭秤资料共源。|

原文PDF文本的极小概率指数OCR损坏。本件使用直接可读的十数和密度式独立积分，得到约$4.29\times10^{-7}$，不把OCR缺字本身当新资料。所有高精度位数属于代数复算，不代表1981测量精度。

## 二、1982评论／回复：访问范围明确留账

|原始条目|本次真正读到的范围|
|---|---|
|Bruce Hawkins, *Indirect Evidence for Quantum Gravity?*, PRL **48, 520**；[DOI](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.520)|书目信息及原刊镜像的搜索索引文字可读；涉及测量解释及选择子系综。下列镜像全文抓取超时，作者本次不称完整读过该页；不将其意见当已证定理。|
|Leslie E. Ballentine, *Comment on “Indirect Evidence for Quantum Gravity”*, PRL **48, 522**；[DOI](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.522)|元数据已核，原正文需要授权或未能从公开镜像取得；不引用其具体论证或二手转述的语句。|
|Don N. Page, *Page Responds:*, PRL **48, 521**；[DOI](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.521)|元数据已核，收到1981-11-17、发表1982-02-15；正文未取得。|
|Don N. Page, *Page Responds:*, PRL **48, 523**；[DOI](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.523)|元数据已核，收到1981-12-01、发表1982-02-15；正文未取得。与p.521是两篇单署名回复，不写成Page–Geilker合著的一篇。|

[Hawkins原刊镜像入口](https://electronicsandbooks.com/edt/manual/Magazine/P/Physical%20Review%20Letters/Physical%20Review%20Letters%201982-1986/root/data/PhysRevLett%201982-1986/pdf/PRL/v48/i7/PRL_v48_i7_p520_1.pdf)仅记录访问；不以存在URL代替全文实读。并行独立代理亦核了原刊条目与访问限制，随后从索引补读Hawkins；本作者不借此宣称自行读完另外三页。

因此本件**没有完成1982争论的全文裁决**。正文所用的窄结论与排除边界均可直接定位于1981已读全文；若要新增“某条回复已经解决某个反对”的历史主张，必须另取得相应原文，不能由本件授权。

## 三、项目内复用与自己的分析

- [973](../../../archive_956_989/research_note_973.md)是原实际热准备的有限分支诊断，并未完成Page–Geilker经验采用。其“均值相同不保证联合任务相同”结论直接复用，不重计科学。
- [1053](../../research_note_1053.md)、[1056](../../research_note_1056.md)分别处理电子Newton源与弱衰变总角源，不当作1981宏观来源／读数的完整量子仪器证明。
- [共同范围校正](../joint_scope_after1062.md)及[公共分支审计](../common_branch_scope_after1062.md)禁止把局部成功拼成全局存在，也禁止把全装置闭合作为每份有限证书的前置要求。
- 正文式(1)—(3)是本作者显式写出的弱场、线性读口概率重述。它在共同核的条件下区分联合律，未从原十数识别完整量子态或未知噪声核。
- Pearson复算仅从出版摘要出发；连续零相关分布的噪声假设与条件置换的交换性假设各自列明。后者是本件诊断，不归给1981作者，也不新增独立经验资料。

## 四、交付与执行

[check.py](check.py)使用标准库Fraction及65位Decimal：精确求相关平方、检查零相关密度归一、积分单／双尾；枚举全部252个五五标签分配。首次保存用排他创建，默认只读与[results.json](results.json)逐项比较。不做图像处理、不拟合原始事件、不重建扭秤或量子控制器。
