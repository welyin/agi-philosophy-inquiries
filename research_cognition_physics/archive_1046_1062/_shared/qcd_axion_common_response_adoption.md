# QCD轴子：同一θ来源的质量、核子响应与暗态振幅

2026-10-09。**成熟采用计0；不新增1063或经验试验组。** 这是一个可替换的物理扩展分支，不把轴子指定为唯一暗物质或强CP解释，也不修改P981原字段菜单。

## 1. 与历史的准确连接

[629](../../archive_629_652/research_note_629.md)已证明指定全局规范群／质量相位的共同商，并明确未决定强CP参数。本次局部QCD有效字典不覆盖该全束周期／边界结论。[930](../../archive_923_934/research_note_930.md)和[1040](../../archive_1009_1043/research_note_1040.md)的电磁“axion型”响应系数本身不证明存在动力学粒子；[1008的P14](../../archive_990_1008/1008/overall_operation_hypothesis_v2_2.md)及[1043](../../archive_1009_1043/research_note_1043.md)采用的是另一可替换暗标量分支。

本次新增共同预测关系：若选择标准QCD轴子，则势曲率、振荡频率、同θ引起的核子CP奇响应，以及声明局部暗态的振幅，必须使用同一规范化。它们不能被当作四个独立拟合量。

## 2. 一个θ字典，两个不同的QCD匹配量

取正则归一的实字段$a$，用有效常数$f_a>0$吸收色反常系数，定义

$$
 \theta_{\rm eff}=\bar\theta+\frac{a}{f_a},\qquad
 \bar\theta=\theta_0+\arg\det M_q,\qquad
 \mathcal L\supset\theta_{\rm eff}\frac{\alpha_s}{8\pi}G\widetilde G . \tag{1}
$$

这是所选低能CP／符号约定。它不选择PQ全局群、域壁数或629的全局商；有文献另保$C_G/f_a$时，本文$1/f_a$对应那个比值。

采用QCD主导势、无另加显著显式PQ破缺，并忽略本合同外的弱CP残余位移。围绕所选CP守恒极小值令$\vartheta=\theta_{\rm eff}$，小振幅时

$$
 V_{\rm QCD}(\vartheta,0)-V_{\rm QCD}(0,0)
 =\tfrac12\chi(0)\vartheta^2+O(\vartheta^4),\qquad
 m_a^2 f_a^2=\chi(0).                                  \tag{2}
$$

$\chi$是同一QCD真空的拓扑易感率；式(2)使用正则动能及领先低能匹配。一般零温势并非精确单余弦。成熟推导见[1511.02867v2，式(1)、(14)、(22)](https://arxiv.org/html/1511.02867v2)。

同一小、慢变$\vartheta$对自由中子产生的θ型EDM写为

$$
 d_n^{(\theta)}(t)=\kappa_n\vartheta(t)+O(\vartheta^3)
 +\text{非绝热／高阶匹配修正}.                          \tag{3}
$$

$\kappa_n$是需独立QCD／强子计算的响应系数，**不能从$\chi$一个数推出**。这里也不把轴流导数耦合、其他CP来源或整个核／原子的EDM等同自由中子EDM。原[2017实验论文式(2)](https://arxiv.org/pdf/1708.06367v1)采用约$2.4\times10^{-16}\ e\,\mathrm{cm}$的θ响应估计；本文保留$\kappa_n$符号，不将该估计当成无误差常数。[CASPEr原始提案§I、V](https://arxiv.org/html/1306.6089v2)已明确同QCD来源的振荡EDM关系。

## 3. 固定局部暗态时，振幅与频率如何联动

以下用$\hbar=c=1$。限定低温、非相对论、近谐振且在一个有限相干时空域内可用单模近似的局部态：

$$
 a(t)-a_{\min}=A\cos(\omega t+\varphi),\quad
 \omega\simeq m_a,\quad
 \rho_a^{\rm loc}\simeq\tfrac12m_a^2A^2 .             \tag{4}
$$

零动量谐振近似中这些等式成立；速度、梯度、非谐项及相干丢失分别限制其精度。将(2)—(4)相接，θ型EDM的峰值振幅为

$$
 \vartheta_A=\sqrt{\frac{2\rho_a^{\rm loc}}{\chi(0)}},\qquad
 D_n=|\kappa_n|\sqrt{\frac{2\rho_a^{\rm loc}}{\chi(0)}},\qquad
 \omega\simeq\frac{\sqrt{\chi(0)}}{f_a}.               \tag{5}
$$

因此固定$\rho_a^{\rm loc},\chi,\kappa_n$后，改变$f_a$不独立压低这个裸响应振幅，却改变频率、相干时间和仪器可见度。这是成熟关系，不是本项目新定理。若取$\rho_a^{\rm loc}=\xi\rho_{\rm DM}^{\rm loc}$，振幅按$\sqrt\xi$变；不能把$\xi=1$免费输入后宣称已识别暗物质。

普通ALP若有独立质量来源、不同CP算符或多场混合，未必满足(2)及(5)。这提供一个限定分支内的相容性判据，不能由它排除所有轴子／暗物质候选。

## 4. 裸振幅不等于有限时间读数

在声明的理想自由自旋Ramsey合同中，恒定共线电场的θ型作用为$H_{\rm EDM}=-d_n(t)E\sigma_z$。它给出的两自旋分量相对相位为$2E\int d_n(t)dt/\hbar$，相位正负随分量排序而变。对矩形自由演化窗$[t_0,t_0+T]$，

$$
 \bar d_T=\frac1T\int_{t_0}^{t_0+T}d_n(t)dt
 =D_n\,\operatorname{sinc}\!\left(\frac{\omega T}{2}\right)
 \cos\!\left[\omega(t_0+T/2)+\varphi\right],\quad
 \operatorname{sinc}x=\frac{\sin x}{x}.                \tag{6}
$$

这里吸收了$\kappa_n$符号到相位中。单次相位未知或窗长恰在零点，裸振幅非零仍可给零读数。实际采样时刻、脉冲宽度、电场换向、去趋势、噪声和相干谱给不同滤波，不能用(6)认证整个实验；超过相干域应改用多模／随机谱任务。核与原子材料还需Schiff屏蔽等响应，不能把$d_n E$直接搬进任意中性材料。

作为已有实际读口的固定历史采用，[Abel等1708.06367v1，§II—III](https://arxiv.org/pdf/1708.06367v1)用ILL／PSI中子—${}^{199}Hg$频率比寻找振荡，读到的是$d_n-(\mu_n/\mu_{Hg})d_{Hg}$，并保磁场系统项。作者在其超低质量搜索范围未见满足信号判据的信号，给频率依赖限制；耦合解释采用局部密度$0.4\ \mathrm{GeV/cm^3}$。这不是所有频率等敏感，也不是本项目重新验出的QCD轴子排除。我们不读取曲线重造限制、不合并旧静态EDM数据作独立证据，也不把2014提案的灵敏度当成已完成实验。

## 5. 宇宙丰度与局部准备仍须另给

在均匀、忽略额外耗散的宇宙近似中，同一规范化给

$$
 \ddot\vartheta+3H\dot\vartheta
 +f_a^{-2}\partial_\vartheta V_{\rm QCD}(\vartheta,T)=0. \tag{7}
$$

初始位移、$\chi(T)$及势形、膨胀／熵历史、PQ破缺时序和可能的缺陷决定产生；星系形成与聚集还影响局部密度／速度。故(7)不从质量关系单独给出$\Omega_a$，更不能自动把宇宙平均丰度换成实验室$\rho_a^{\rm loc}$。采用[1511.02867v2 §3.3](https://arxiv.org/html/1511.02867v2)的这项依赖结构，不重现其年代相关参数窗口或声称其热QCD误差代表当前精度。有限温势与热背景之间的能源交换也不能凭(7)当成已闭合的整体来源解。

## 6. 验收范围与停止线

[代码](qcd_axion_common_response/check.py)／[结果](qcd_axion_common_response/results.json)只用有理缩放核质量—振幅恒等，并将矩形窗解析式与有限求积比较。数值用无量纲代入，不是实际轴子参数、原始资料或严格QCD误差证书。[来源表](qcd_axion_common_response/sources.md)保留版本及复用位置。

本次增加M3B／M5中同一θ参数的物理相容关系，以及M6中源振幅—真实读口滤波的区分。没有生成PQ机制、$f_a$、核匹配常数或暗态；未交付全量子未知输入仪器、完整引力反作用或宇宙产生模拟。不将本分支并入D991或声称统一目标完成。**到此停止，不扫描参数、不建探测器、不以候选UV完成作为原目标门槛。**
