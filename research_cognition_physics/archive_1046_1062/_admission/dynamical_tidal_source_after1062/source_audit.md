# 动态潮汐原文的符号与量纲审计

2026-10-09。独立来源审计，计0，不立新轮，不修改既有静态潮汐采用或原文。对象为 Apostolidis 等 *Dynamical tidal response of neutron stars: From effective field theory to gravitational waveforms*，[arXiv:2606.19446v2 PDF](https://arxiv.org/pdf/2606.19446v2)与[同版HTML](https://arxiv.org/html/2606.19446v2)，版本日期2026-09-22。实读相关定义及公式；以下是本项目的代数核对，**不是作者发布的勘误，也不是对该理论的否定**。

## 1. 已确认的两处不一致

PDF印刷页9、12—14的式(56)—(71)与HTML一致，不是网页转换差异。记紧致度为C，另用不同名字保留两种符号约定：

$$
\beta_{\log}=2C^{-8}c_{\dot E}^{(\log)}
=-\frac{64}{45}\left(1+\frac{321}{112}\Lambda_0\right).
$$

式(56)、(59)给$\Lambda_0=16b_0/(45a_0)$，故上述式(66)系数与式(57)相符。但是直接代入$\mu=\gamma/r$必须得到

$$
\Lambda_2(r)=\bar\Lambda_2+
\beta_{\log}\log(\gamma r_s/r).
\tag{A}
$$

式(65)却使用相同β的减号；原文明确把式(66)的β定义为对数系数，因而不能靠默认相反的β定义解消。有限代数见证取$b_0/a_0=1$、$\bar\Lambda_2=10$、$\log(\mu r_s)=-1$，两式差为$27136/4725$。这是公式身份的见证，**未声称该组系数来自已积分的恒星解**。

另外，解式(63)得到

$$
c_{\dot E}\simeq\frac{GM}{R_\star^3}
\frac{c_E}{\omega_f^2}.
\tag{B}
$$

式(64)印作$GM/R_\star$，缺少分母$R_\star^2$。式(6)声明两个c无量纲，式(62)使用有频率量纲的$\omega_f$；未见此处另行定义$R_\star\omega_f$。正确右端长度次数为0，印刷右端为2。

## 2. 波形段自身可保持一致，但不能无修改地跨接

把式(65)本身当作定义，另命名其参数为$\beta_{\rm wave}$：

$$
L(r)=\bar\Lambda_2+\beta_{\rm wave}\log[r/(\gamma r_s)],
\qquad rL'(r)=\beta_{\rm wave}.
\tag{C}
$$

独立的保留阶圆轨道代数核对式(67)→(68)，以下取$G=M_{\rm tot}=1$：一体动态项写作$A L(r)\omega^2/r^6$。其直接能量与圆轨道半径改变量共同给

$$
\frac A{r^9}\left[L+\frac23(6L-\beta_{\rm wave})\right]
=\frac A{r^9}\left(5L-\frac23\beta_{\rm wave}\right).
$$

再以式(68)、(69)为输入，在线性潮汐保留阶用$\psi''=-2\mathcal E'/\mathcal F$核对式(70)、(71)。设$G=M_{\rm tot}=1$、$a=m_1/M_{\rm tot}$、$b=m_2/M_{\rm tot}$，一体修正为

$$
\delta\psi''=\frac{5a^6\omega^{5/3}}{48\eta}
\left[-3(8+147b)L+105b\beta_{\rm wave}\right].
$$

式(70)—(71)二阶求导给同一式；特别是$38/11$、$1253/44$的常数项及对数项通过精确分数核验。本检查没有重新独立推导全部辐射通量式(69)，也没有验收全套波形或Fisher计算。

因此，式(65)—(71)的这一内部代数与式(57)、(66)跨接时，需$\beta_{\rm wave}=-\beta_{\log}$。这是本次审计提出的一致记号；**不能声称原文已经如此定义**。另一等价处理是改式(65)的符号并同步翻转后续所有β项；仅改单个式子不足。

## 3. 本项目采用边界

- 可保留指定方案下零频与$\omega^2$匹配所需的独立系数、式(56)—(59)的上述代数及式(B)的尺度关系。原文的非旋转完美流体保守分支不提供任意内部弛豫或耗散匹配。
- 本批**不同时采用原样的式(57)、(65)、(66)**，不领取完整8PN／Fisher数值、严格有限频带物理余项或实际探测证书。系数符号审计本身不产生新的认知选择结论。
- [既有静态TOV／Love采用](../../_shared/neutron_star_tidal_adoption.md)不受本次局部审计撤销。是否新增动态共同接口须另作准入，不能仅因发现排印不一致就替代该工作。

## 4. 可复算资产

[check.py](check.py)仅用Python标准库精确分数，输出见[results.json](results.json)。默认只读重算并逐项比较；首次仅以`--write-exclusive`创建结果，文件已存在即拒绝覆盖。未重算恒星、波形、实验或径向积分，也未修改历史资产。
