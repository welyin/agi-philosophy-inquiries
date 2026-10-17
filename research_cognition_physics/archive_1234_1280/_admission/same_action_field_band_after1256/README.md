# 同一作用的场涨落有限频带：1256后有界准入筛查

日期：2026-10-16。**非编号方向筛查，不是新轮次准入、完整量子模型或热参考验收。** 本目录只增加这一项有限结论，不修改冻结材料、导航和阶段目标。

## 1. 对象、去重及所固定的东西

使用[1247原作用与两支实际历史](../../1247/proof.md)，在其已证明的保持窗 $I=[2,9/4]$ 上计算 $p,X,\phi^A$ 的二次变分。原[1250](../../1250/proof.md)沿同一历史核算完整来源及有限 null 余额；这里沿用其保持期方程，不重新积分背景。

本筛查**固定原度规和两流体背景，只取场—场 Hessian 子块**。二次参考近似使用原场及原耦合，未增加物质；但这不代表把完整原模型的其他变量删去后仍得到一份满足 Einstein 约束的物理解。一般场扰动改变来源，并激发度规、第二流密度/作用量及机械自由度。它们没有在本筛查被量子化或消元；不能由该子块的正性推论完整约束扰动正性。

已查重：

- [1159](../../../archive_1159_1160/1159/proof.md)和[1163](../../../archive_1161_1200/1163/proof.md)给有限实际读口与共同主部判据；未给本材料的 $F_p$—绕数杆混合矩阵。
- [312](../../../archive_301_341/research_note_312.md)已有 boost/球形流及热参考条件的区分；本筛查不重证 KMS/Unruh。
- [522](../../../archive_467_530/research_note_522.md)是另选的 $R/M$ 谐振材料和热参考。其 Gaussian 方法可复用，不能把其非紧 $R$ 零模结论搬到原 compact 杆。
- [1249后原 $X$ 参考筛查](../quantum_reference_bridge_screen_after1249.md)已给同背景 SLE 方向，并指出实际准备、有限探测响应及重整化来源尚缺。这里的净增仅是原 $p$—杆混合在实际时变背景上的明确正频带；没有再次发现“可以选某个 Gaussian 参考”。
- [1251后同标签指针](../local_collective_record_connection_after1251.md)已完成经典指针、钟补偿及完整源的桥，不把它重新算作未完成连接。

## 2. 同源二次作用及共同主锥

原场材料为

$$
\mathcal L_f=-\frac12(\nabla p)^2-\frac12p^2
-\frac12(\nabla X)^2-\frac12X^2
-\frac12F(p)\sum_{A=1}^3(\nabla\phi^A)^2,
\qquad F=1+\frac p5.
\tag{1}
$$

背景为 $ds^2=-dt^2+a^2d\mathbf x^2$、$u=a^{-1}$、$\bar\phi^A=kx^A$、$k^2=1/10$；$\bar p(t),\bar X(t)$ 取原两支实际解。三杆是周期 $kL$、$L=10$ 的 compact shift 场，原绕数没有改变。令 $s=\delta p$、$x=\delta X$、$\pi^A=\delta\phi^A$。固定原流体背景时，保持期第一流不耦合 $p,X$，第二流质量对 $p$ 线性，所以没有额外的纯场二次势 Hessian。原线性来源项仍由背景方程抵消。

二次作用精确为

$$
\begin{aligned}
S_f^{(2)}=\frac12\int dt\,d^3x\,a^3\bigg[
&\dot s^2-a^{-2}|\nabla s|^2-s^2
+\dot x^2-a^{-2}|\nabla x|^2-x^2\\
&+\bar F\sum_A\big(\dot\pi_A^2-a^{-2}|\nabla\pi_A|^2\big)
-\frac{2k}{5a^2}s\sum_A\partial_A\pi_A
\bigg].
\end{aligned}
\tag{2}
$$

$F_{pp}=0$，故杆背景不贡献 $s^2$ 项；$\partial_t\bar\phi^A=0$，故没有 $s\dot\pi$ 的背景混合。混合项只有一阶空间导数，主二阶符号为

$$
\operatorname{diag}(1,1,\bar F,\bar F,\bar F)
\,g^{\mu\nu}\zeta_\mu\zeta_\nu .
\tag{3}
$$

在 $F>0$ 内，该场子块共享原 $g$ 的 null 主锥。这不是从认知原则新推出共同度规，也不是证明所有尘流/仪器均有同一波动频率。

## 3. 冻结系数检验及其限制

取非零共动波矢 $\boldsymbol\xi=2\pi\mathbf n/L$，记 $\mathsf q=|\boldsymbol\xi|/a$ 为物理波数；它不是原钟内部的振子变量 $q$。将杆扰动按 $\boldsymbol\xi$ 分成一个纵向和两个横向方向，冻结 $a,F$ 后，规范化纵向杆与 $s$ 的频率矩阵为

$$
\overline\Omega_{pL}^2=
\begin{pmatrix}
\mathsf q^2+1&i\gamma\mathsf q\\
-i\gamma\mathsf q&\mathsf q^2
\end{pmatrix},
\qquad \gamma=\frac{ku}{5\sqrt F}.
\tag{4}
$$

Fourier 符号约定相反只同时翻转两个非对角号。其本征值为

$$
\overline\omega_\pm^2=\mathsf q^2+\frac12
\pm\sqrt{\frac14+\gamma^2\mathsf q^2},
\qquad
\det\overline\Omega_{pL}^2
=\mathsf q^2(\mathsf q^2+1-\gamma^2).
\tag{5}
$$

横杆为 $\mathsf q^2$，$X$ 为 $\mathsf q^2+1$。由后述背景界 $\gamma^2<1/175<1$，所有非零波数的冻结子块正定。但冻结本身不产生实际时变背景的 KMS 态，故还必须核下一节的时间项。

## 4. 实际时变背景的频率矩阵

令 $b=a^{3/2}$，取

$$
Q_p=bs,\qquad Q_X=bx,\qquad Q_A=b\sqrt F\,\pi_A.
\tag{6}
$$

删除明确的时间全导数后，动能为单位矩阵。定义

$$
U_p=\frac{\ddot b}{b}=\frac32\dot H+\frac94H^2,
\qquad
U_r=\frac{(b\sqrt F)^{\cdot\cdot}}{b\sqrt F}.
\tag{7}
$$

原保持期方程给 $H^2=\rho$、$\dot H=-3(\rho+P)/2$、
$\ddot p=-3Hv-p-3u^2/100-u^3/100$。因此

$$
\boxed{U_p=-\frac94P,\qquad
U_r=U_p-\frac{p+3u^2/100+u^3/100}{10F}
-\frac{v^2}{100F^2}.}
\tag{8}
$$

特别地，$3H\dot F/(2F)$ 与 $\ddot F/(2F)$ 中的 $Hv$ 项抵消，不能只保留一个。实际时间相关矩阵为

$$
\Omega_{pL}^2(t)=
\begin{pmatrix}
\mathsf q^2+1-U_p&i\gamma\mathsf q\\
-i\gamma\mathsf q&\mathsf q^2-U_r
\end{pmatrix};\quad
\Omega_X^2=\mathsf q^2+1-U_p,\quad
\Omega_{T_1}^2=\Omega_{T_2}^2=\mathsf q^2-U_r.
\tag{9}
$$

这是式(2)的实际时变二次系统，不是令 $H=\dot F=0$ 的近似。使用时变正则变量会改变瞬时 Hamiltonian 的解释；这里的“频率正定”只指已明确选择的式(6)—(9)。

## 5. 不重跑ODE的有理数正定界

直接复用1247保持窗证明：$\rho\le111/100$、$|p|<3/2$、$F>7/10$、$0<u\le1$，原完整尘质量正。正能标量、杆和尘的压力满足 $|P|\le\rho$。所以在两支整个 $I$ 上，

$$
U_p\le\frac{999}{400},\qquad
U_r\le\frac{999}{400}+\frac3{14}<\frac{339}{125},\qquad
\gamma^2<\frac1{175}<\left(\frac2{25}\right)^2.
\tag{10}
$$

对全部 $\mathsf q\ge2$，Hermitian 行和下界给五场矩阵共同的严格余量

$$
\lambda_{\min}\big(\Omega^2(t)\big)
>\mathsf q^2-\frac{339}{125}-\frac2{25}\mathsf q
\ge\frac{141}{125}>0.
\tag{11}
$$

最后一步可写为 $\mathsf q=2+y$：右边减去 $141/125$ 等于 $y^2+98y/25\ge0$。固定 $\mathsf q=2$ 的纵向块也以式(10)分别核主对角及行列式正性；代码保存精确分数。

这不是依赖无限波数才存在的域。原结果的两支 $u(2)<3/8$，保持窗又有

$$
\frac{29}{40}\frac{364929449}{10^9}<u(t)<\frac38.
$$

固定一份非零模式对 $\mathbf n=\pm(13,0,0)$，用 $3.14<\pi<22/7$ 即得整个窗口

$$
2<\mathsf q_{13}(t)<\frac{429}{140}<\frac{13}{4}.
\tag{12}
$$

所以一份固定有限模式对已经有严格正频带，不需要改空间维数、积分新背景或送到 UV 极限。模式对指实场的正弦/余弦自由度，不能把 $\pm\mathbf n$ 当两个独立复实场。

## 6. 它允许什么参考，尚未得到什么

对任意事先固定的有限模式集，式(9)是连续系数的有限线性 Hamilton 系统。在 $t=2$ 可数学地选其正定瞬时频率矩阵对应的 Gaussian 协方差，再由同一式(9)的辛基本解运输；有限窗系数有界，协方差有限。这只是明确的有限参考候选，实际初准备尚未给出。正定不排除时变参数产生放大，不证明绝热不变、平衡或温度。

以下均未由本筛查获得：

1. 完整 Einstein 约束的量子化、度规与两流相对模的共同稳定性；固定背景的场 Hessian 不是封闭的全物理涨落扇区。
2. 原实际钟、局部 $K$ 末读、局部仪器的准备与操作误差。有限 Fourier 带本身也不是紧支局域器件。
3. 完整参考态的 Hadamard/重整化来源证书、反作用回写或同一参考下的有限探测响应。二次场贡献不能未经核算直接加到旧经典源后保持旧几何不变。
4. KMS、Unruh、模 Hamiltonian 或 Jacobson 热—几何身份。冻结系数和正频带均不授予这些权限。
5. compact 杆的整个零模扇区。局部切空间 Gaussian 不能替代原转子拓扑；522非紧零模的无正规 Gibbs 结论也不能直接用于这里。

**方向判断：** 原 $p$—杆低阶混合在实际保持窗没有迫使这条有限场参考方向立即失败。该候选仍需先跨过已有1249后的实际参考/来源连接；仅把更多二次模式计算出来，不足以宣布新的完整研究轮或共同模型准入。

## 7. 复算、绑定及冻结

`check.py` 只用 Python 标准库 `Fraction`，核式(8)的系数抵消、式(10)—(12)的有理数界、固定 $\mathsf q=2$ 主子式及 $\mathsf q\ge2$ 的多项式余量；从冻结1247结果读两支 $u(2)$，不调用旧积分器。默认只读重算并与 `results.json` 完全比对；`--save` 仅可首次排他创建结果，已存在时拒绝覆盖。

绑定的原文件 SHA256 为：

| 文件 | SHA256 |
|---|---|
| `1247/proof.md` | `f3a79f8f43aa235d2fa13cc87a232e2b64643fc1bade375bfee96c204ab4bf3d` |
| `1247/results.json` | `dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06` |
| `1250/proof.md` | `5e9f5e1495156f6bcd81c48a5046ac538ea72b2823a6cb7a9e3b7a59ff499d9f` |
| `1250/results.json` | `84fa683b1e343f94ba2f5a4f2e26653dfa49c7d6c04d3833a436daa8cb880284` |

来源证明承担作用量变分、背景界和范围量词；有限分数检查只审计其具体代数应用，不冒充完整QFT证明。
