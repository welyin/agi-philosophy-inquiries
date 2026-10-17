# 静态 de Sitter 参考：正常相干态、固有时与完整共形来源荷

日期：2026-10-17。非编号成熟定理接入；不改目标、认知公理、正式计数或旧证据。单位取 $c=\hbar=k_B=1$，度规号差为 $(-+++)$。本稿证明固定背景上的同对象身份，不认证实际准备器、有限温度计、完整反作用或认知原则对 de Sitter 的选择。

## 1. 复用、净增和命题范围

- [318，§2—5](../../../archive_301_341/research_note_318.md)：改进能流、Wald 接触和两端项必须共同运输；其例子为经典平直类光波形。
- [1004，§4](../../../archive_990_1008/research_note_1004.md)：相干位移的两点函数与应力差、准备成本必须保留；该段明确未另认证任意全局 Fock 表示中的酉实现。
- [既有 de Sitter 有限字典，§3、5—6](../thermal_mechanism_admission_after1270/de_sitter_finite_dictionary_after1270.md)：给物理参数、完整改进源及有限面积端点；明确尚未从其 $C^2$ 类光数据领取正常相干态/完整模 Hamiltonian。
- [既有候选，§4.2—5](../thermal_mechanism_admission_after1270/actual_thermal_reference_candidate_after1270.md)：BD 参考、中心固有时和非零激发已经被列为候选。

本次新增的项目内接入是：给出一份**紧支光滑 Cauchy 初值 → patch 内紧支 Weyl 酉 → BD 表示中的正常相干态 → 同一静态时间的模荷 → 物理共形 Hilbert 来源荷**的明确链。其核心模流定理来自成熟文献，不能重报为新物理定理；边界运输及下面有限例子是本稿的直接计算。没有重新证明一般相干两点函数或平均来源。

采用 [Fröb–Much–Papadopoulos，arXiv:2310.12185v2](https://arxiv.org/pdf/2310.12185) 的自由场楔区结果：式(66)、(70)—(78)把正常相干态相对熵写成同一 Killing 流的经典来源荷；§II 的楔定义与 §III 的模流固定其区域和归一。仅调用这一自由场结果，不调用其对相互作用场的展望。

## 2. 完整楔、静态 patch 与中心固有时

四维完整 de Sitter 的一个右楔与静态 patch 等距。取 $N(r)=\sqrt{1-H^2r^2}$，

$$
ds^2=-N^2d\tau^2+N^{-2}dr^2+r^2d\Omega_2^2,\qquad 0\le r<H^{-1}.
\tag{1}
$$

可用嵌入坐标直接核验：

$$
X^0=\frac{N}{H}\sinh(H\tau),\quad
X^j=\frac{N}{H}\cosh(H\tau),\quad
\boldsymbol X_\perp=r\boldsymbol n .
\tag{2}
$$

它满足 $X^j>|X^0|$，即原文式(13)的右楔。原文 boost 参数 $s=H\tau$，未来向量 $\xi_{\rm FMP}=-M_{0j}=H^{-1}\partial_\tau$。本稿采用 $\zeta=\partial_\tau$，它在 $r=0$ 的范数为 $-1$，所以 $\tau$正是中心测地线固有时。完整楔不是原文某一 Poincaré 坐标片中截下的有限部分。

定义正号物理荷

$$
Q_T[\Sigma;\zeta]=\int_\Sigma n^\mu\zeta^\nu T_{\mu\nu}\,d\Sigma ,
\tag{3}
$$

其中 $n^\mu$未来向。令 $\omega_{\rm BD}$为同一自由场的 BD 参考；对 patch 内的局部 Weyl 酉 $W(f)$，记 $\omega_f$为其相干态、$\varphi=Ef$为均值。原文式(66)、(70)—(78)的归一给

$$
S_{\mathcal A(\mathrm{patch})}(\omega_f\Vert\omega_{\rm BD})
 =2\pi Q_{T^{\rm eff}}[\Sigma;\xi_{\rm FMP}]
 =\frac{2\pi}{H}Q_{T^{\rm eff}}[\Sigma;\zeta].
\tag{4}
$$

因此 $\beta_0=2\pi/H$。这里 $Q$为相对同一参考的荷，式(4)不是全空间真空总能量。局部静态观察者的逆温度为 $\beta_{\rm loc}(r)=N(r)\beta_0$；不把它与中心归一混用。恢复 $\hbar$后，式(4)右侧为 $2\pi Q/(\hbar H)$。

## 3. 同一场方程不等于同一度规来源：完整运输

实际选择的物质作用是无质量共形标量

$$
I_{\rm phys}
=-\frac12\int\sqrt{-g}\left[(\nabla\Phi)^2+\xi_c R\Phi^2\right]d^4x,
\qquad \xi_c=\frac16 .
\tag{5}
$$

在固定 de Sitter 上 $R=12H^2$，其场方程与有效最小耦合质量 $m_{\rm eff}^2=2H^2$的 KG 方程相同。其辛结构和 BD 自由场表示也相同，可以调用式(4)。但度规变分时，不能把 $\xi_cR$当作一个固定质量再声称应力完全相同。

对 $\varphi=Ef$，令 $F=\varphi^2$。零均值 BD 参考与相干态按同一固定背景、同一重整化处方相减，几何独立于状态的项抵消，物质来源差为相应经典 Hilbert 应力：

$$
T^{\rm phys}_{\mu\nu}
=\nabla_\mu\varphi\nabla_\nu\varphi-\frac12g_{\mu\nu}(\nabla\varphi)^2
+\xi_c\left(G_{\mu\nu}F+g_{\mu\nu}\Box F-\nabla_\mu\nabla_\nu F\right),
\tag{6}
$$

$$
T^{\rm eff}_{\mu\nu}
=\nabla_\mu\varphi\nabla_\nu\varphi
-\frac12g_{\mu\nu}\left[(\nabla\varphi)^2+2H^2F\right].
\tag{7}
$$

故差是

$$
D_{\mu\nu}:=T^{\rm phys}_{\mu\nu}-T^{\rm eff}_{\mu\nu}
=\xi_c\left(R_{\mu\nu}F+g_{\mu\nu}\Box F-\nabla_\mu\nabla_\nu F\right).
\tag{8}
$$

对任一 Killing 向量 $\zeta$，定义反对称量

$$
\boxed{\mathcal B^{\mu\nu}
=\xi_c\left(\zeta^\mu\nabla^\nu F-\zeta^\nu\nabla^\mu F
+F\nabla^\mu\zeta^\nu\right).}
\tag{9}
$$

直接微分：第一、第三项含一阶 $\nabla F\,\nabla\zeta$的部分由 Killing 方程相消；$\nabla_\nu\zeta^\nu=0$，且
$\nabla_\nu\nabla^\mu\zeta^\nu=R^\mu{}_\rho\zeta^\rho$。因此

$$
\boxed{D^{\mu\nu}\zeta_\nu=\nabla_\nu\mathcal B^{\mu\nu}.}
\tag{10}
$$

这个式子的 Ricci 项符号与原文采用的 $R_{\mu\nu}=3H^2g_{\mu\nu}$一致。它是两个独立来源定义之间的表面运输，不是将 $T^{\rm eff}$重新命名为物理来源。

在静态切片的球 $0\le r\le b$上，若数据径向，式(9)给

$$
\mathcal B^{\tau r}
=\xi_c\left(N^2F'+H^2rF\right).
\tag{11}
$$

这里 $\nabla^\tau\zeta^r=H^2r$，$n_\tau=-N$，$d\Sigma=r^2N^{-1}dr\,d\Omega$。于是正荷约定(3)给**负的外边界项**

$$
\boxed{Q_{\rm phys}(b)-Q_{\rm eff}(b)
=-4\pi\xi_c b^2\left[(1-H^2b^2)F'(b)+H^2bF(b)\right].}
\tag{12}
$$

一般非径向数据有相应表面积分。若一份 Cauchy 切片上的初值及导数紧支于 patch 内部，则完整切片外边界项为零，故式(4)真正接成

$$
\boxed{S_{\mathcal A(\mathrm{patch})}(\omega_f\Vert\omega_{\rm BD})
=\beta_0 Q_{\rm phys}[\Sigma;\zeta].}
\tag{13}
$$

只有完整区域及这一支撑条件才得到式(13)。有限子球的荷不是该子球的完整模 Hamiltonian；不能把式(13)原样限制到任意有限球。有限类光截面也必须保留同作用的 Wald 接触项 $-2\pi\xi_c\int F\,dA$及端点运输，具体字典直接复用既有 de Sitter 有限字典§5，不能只留面积项。

## 4. 明确的正常紧支准备：不只指定类光波形

现在取 $H=1$，在 $\tau=0$静态 Cauchy 切片给定

$$
\varphi(0,r)=h(r)=
\begin{cases}
\exp\!\left(1-\dfrac1{1-4r^2}\right),&r<\dfrac12,\\
0,&r\ge\dfrac12,
\end{cases}
\qquad \partial_\tau\varphi(0,r)=0 .
\tag{14}
$$

它在中心是 $r^2$的光滑函数，在 $r=1/2$各阶导数为零，因而是紧支 $C^\infty$初值；在全局 Cauchy 面的其余部分取零也光滑。线性 KG Cauchy 演化给唯一光滑解，其有限时间支撑由有限传播速度控制。静态径向光速 $|dr/d\tau|=1-r^2\le1$；在 $|\tau|\le1/16$，该初值的因果支撑包含在 $r\le9/16<1$内。

令 $\chi(\tau)$光滑地由过去的0变为未来的1，导数支撑于 $(-1/16,1/16)$。令 $P=\Box-2$，在该管中置

$$
f=P(\chi\varphi)=[P,\chi]\varphi .
\tag{15}
$$

其支撑紧含于 patch；延拓为零得到 $f\in C_c^\infty(\mathrm{patch})$。采用 $E=E_{\rm ret}-E_{\rm adv}$时，$E_{\rm ret}f=\chi\varphi$、$E_{\rm adv}f=(\chi-1)\varphi$，故 $Ef=\varphi$。这些式子按 patch 的传播子理解；因果兼容的全局延拓给同一局部初值。这是标准 time-slice 构造，不是新造时空依赖控制器。

因此 $W(f)$是 BD 表示中的有界局部 Weyl 酉，$W(f)\Omega_{\rm BD}$给正常向量态，也保 Hadamard 奇性。相对参考的静态荷有限。这里的“正常准备”只说明这份态在指定表示中实际存在；没有证明由许可控制器在有限时间无代价制备。把式(15)直接当作外加 $j\Phi$的脉冲时，脉冲的应力、交换功与后态必须另外保留。

## 5. 两个可复算的有限结果

在切片(14)，完整有效静态荷为

$$
Q_{\rm eff}
=4\pi\int_0^{1/2}r^2
\left[\frac12(1-r^2)(h')^2+h^2\right]dr .
\tag{16}
$$

在 $0\le r\le1/4$上，$h^2\ge e^{-2/3}$，故

$$
Q_{\rm phys}=Q_{\rm eff}
\ge4\pi e^{-2/3}\int_0^{1/4}r^2dr
=\frac{\pi e^{-2/3}}{48}>0 .
\tag{17}
$$

这是一份正常态的非零完整来源荷；式(13)因此给严格正相对熵。它不是参考存在性的重复声明。

再取有限子球 $b=1/4$。由(14)直接得

$$
F(b)=e^{-2/3},\qquad F'(b)=-\frac{64}{9}e^{-2/3},
\tag{18}
$$

代入(12)有

$$
\boxed{Q_{\rm phys}(1/4)-Q_{\rm eff}(1/4)
=\frac{77\pi}{288}e^{-2/3}>0.}
\tag{19}
$$

因此删除该有限边界项会产生这个明确非零误差；乘中心 $\beta_0=2\pi$可比较带相同时间归一的荷，但不把它称为子球相对熵。完整球中 $h$及导数在外边界为零，所以式(19)与完整荷相等没有矛盾。

该初值还确有非零完整类光来源。中心的局部惯性坐标中，
$h(0)=1$、$\partial_i h=0$、$\partial_i\partial_jh=-8\delta_{ij}$；KG 方程给 $\partial_\tau^2\varphi=-26$。对中心 $k=(1,1,0,0)$，

$$
T^{\rm phys}_{kk}(0)
=-\xi_c\,\partial_k^2(\varphi^2)
=\frac{34}{3}>0 .
\tag{20}
$$

这是固定背景的相减应力，不是已经求出的自洽量子曲率。将初值乘 $\eta$时，式(17)—(20)的荷、来源及相对熵均乘 $\eta^2$。

## 6. 准入判断与停止线

本稿给一项项目内此前未完整写出的共同身份：同一正常 Weyl 相干态、同一静态模流/中心固有时、与同一共形作用的完整 Hilbert 来源荷，在明确支撑条件下相容。它是成熟结果加表面运输的**非编号接入证明**，尚不足以单独构成下一正式研究轮。

以下范围保持：

1. BD、de Sitter、$\xi_c=1/6$及自由作用是本候选的输入；没有从治理机制或认知原则选出它们。原1270的 FLRW 材料、实际钟、端口和资源证书不随更名迁入。
2. 固定参考相减不解决背景变化后的完整真空极化，也未给线性 Einstein 或 CT 条件选择后场方程的有限反作用误差。CT 没有被本稿采用。
3. $W(f)$的正常性与有界性不等于实际准备、观测或交互源已认证；完整模荷也不是一份任意有限探针读数。
4. 有限切片/类光段应另保区域剩余荷、侧通量及式(9)的边界运输；不得将有限子任务直接等同整个 patch 的模任务。
5. 相干族的源和相对熵从 $\eta^2$开始，正负 $\eta$给同一二阶值；不能以此宣称满足317所需的双向非零一阶饱和，或已经选定面积系数。有限或二阶选择问题也没有因此被否定。
6. 对社会机制的用途仅是：同一对象两侧独立定义、交接时边界责任不能删除。它没有检验“完整治理机制迫使特定几何响应”的强猜想，更不能用这份物理例子代替该选择力论证。

本稿不写设备方案，不改旧材料，不追加完整 QFT 或全场 PDE 工作。后继若继续，必须说明实际准备/完整来源和有限几何比较如何在**这同一对象**上成立，或明确停止；不能把(13)重复当作这些身份的替代。
