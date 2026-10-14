# 有限局部多副本背景到完整探针过程的均场误差界

2026-10-13。计算前准入证明，不编号，不采用新公理。本文只处理有限复制模型的严格误差，未签收实际空间运输、共同度规、仪器实现或Einstein动力学。

**结论：在以下有界、交换对称及乘积背景准备条件下，可以得到 $K^{-1/2}$ 的完整准备等距误差；不必退到 $K^{-1/4}$。比较过程保留精确背景联合演化，实际过程保留探针反作用。**

## 1. 对象、输入和有效域

固定有限空间管内的一份完整cluster Hilbert空间 $\mathcal H_c$。它可包含每站的关系自旋、完整有限CAR物质、控制及记录。cluster内的准备允许这些因子彼此纠缠；只要求不同复制之间初始独立：

$$
\Phi_K(0)=\phi_0^{\otimes K},\qquad \|\phi_0\|=1.
\tag{1}
$$

复制标号不是空间坐标。每份cluster包含同一有限空间管；不同复制的相互作用可以只连接同站或既定相邻站。所有空间支持必须在具体模型中事先给出，本证明不会由交换对称性自动推出空间局域性。

在普通复制张量积上，取

$$
H_K=\sum_{r=1}^K h_r+\frac1K\sum_{r<s}W_{rs},
\qquad
h=h^\dagger,\quad W=W^\dagger,\quad
\mathsf S W\mathsf S=W,\quad w=\|W\|.
\tag{2}
$$

$\mathsf S$ 交换两个完整cluster。$h,W$ 和有限空间管不随 $K$ 增长。CAR内容必须用完整有限Fock空间处理；为直接使用这个普通复制张量积证明，可取所有跨复制耦合在每一cluster上均为偶算符、保其奇偶扇区。若另允许复制间费米子交换，则须另核分次交换表示和准备，不能直接套用本页。这里的复制对称性不把cluster内部的费米子改成玻色子。

令

$$
\Phi_K(t)=e^{-itH_K}\phi_0^{\otimes K},\qquad
w^\phi=(I\otimes\langle\phi|)W(I\otimes|\phi\rangle),
\qquad
i\dot\phi=(h+w^\phi)\phi.
\tag{3}
$$

以下取 $\hbar=1$。有限维保证所有演化、导数和算符均定义良好；$\|w^\phi\|\le w$。若 $h$ 随时间连续变化但每份相同，证明亦成立；当前无需外部时变控制。

另取完整有限探针 $P$ 及任意参考 $R$。$P$ 可包含内部钟、位置端口、信号、控制、记录，只要具体模型确实给出这些身份。对任意未知且可与 $R$ 纠缠的 $\chi_{PR}$，背景准备固定为式(1)，并与 $PR$ 独立。实际哈密顿量为

$$
H_{\rm act}=H_K+H_P+\sum_a A_a\otimes\overline B_a,\qquad
\overline B_a=\frac1K\sum_r B_{a,r},
\quad A_a=A_a^\dagger,\quad B_a=B_a^\dagger.
\tag{4}
$$

实际探针—背景耦合和反作用均在同一 $H_{\rm act}$ 中。以下比较器是证明工具：

$$
H_{\rm cmp}(t)=H_K+H_P+\sum_a A_a b_a(t),\qquad
b_a(t)=\langle\phi(t),B_a\phi(t)\rangle .
\tag{5}
$$

它的背景部分是**精确** $\Phi_K(t)$，不把其替成 $\phi(t)^{\otimes K}$。比较器中的 $b_a$ 不被宣称为额外免费物理驱动。

## 2. 激发数估计：一个自含有限维Pickl证明

本节复用Pickl的计数思路；原文第三节、Lemma 3.2给出连续Hartree模型的对应估计，也明确整体向量接近乘积态并非其结论。这里直接证明任意有界有限维 $W$ 的版本，不把原文对乘法势的假设直接照搬。成熟方法来源：[Peter Pickl, A simple derivation of mean field limits for quantum systems](https://arxiv.org/pdf/0907.4464)。

设

$$
p_r(t)=|\phi(t)\rangle\langle\phi(t)|_r,\quad
q_r=I-p_r,\quad
\widehat\alpha=\frac1K\sum_rq_r,\quad
\alpha(t)=\langle\Phi_K,\widehat\alpha\Phi_K\rangle
=\|q_1\Phi_K\|^2.
\tag{6}
$$

最后一个等式用背景交换对称性。取 $K\ge2$，并简写

$$
V_{12}=\frac{K-1}{K}W_{12}-w^\phi_1,\qquad
\|V_{12}\|\le2w.
\tag{7}
$$

对式(6)求导，单体 $h$ 精确消去。其绝对值不超过

$$
2\left|\operatorname{Im}
\langle\Phi_K,p_1V_{12}q_1\Phi_K\rangle\right|.
\tag{8}
$$

在第2因子的两侧插入 $p_2+q_2$。四项如下：

1. $p_1p_2V_{12}q_1p_2$：由
   $p_2V_{12}p_2=-K^{-1}w^\phi_1p_2$，包含式(8)因子2后的界为 $2w/K$。
2. $p_1p_2V_{12}q_1q_2$：用下面的加权Cauchy界。
3. $p_1q_2V_{12}q_1p_2$：$w^\phi_1$ 部分为零；$W$ 部分经交换1、2及取伴随不变，其期望为实数，不贡献虚部。
4. $p_1q_2V_{12}q_1q_2$：两侧分别含 $q_2$，故包含因子2后的界为 $4w\alpha$。

第2项在正激发子空间上插入
$\widehat\alpha^{1/2}\widehat\alpha^{-1/2}$，逆在零激发子空间定义为零。于是

$$
\begin{aligned}
2|\langle\Phi_K,p_1p_2V_{12}q_1q_2\Phi_K\rangle|
&\le
2\|\widehat\alpha^{1/2}V_{12}p_1p_2\Phi_K\|
\|\widehat\alpha^{-1/2}q_1q_2\Phi_K\|,\\
\|\widehat\alpha^{1/2}V_{12}p_1p_2\Phi_K\|^2
&\le4w^2(\alpha+2/K),\\
\|\widehat\alpha^{-1/2}q_1q_2\Phi_K\|^2
&\le\frac{K}{K-1}\alpha\le2\alpha.
\end{aligned}
\tag{9}
$$

第二行将激发数求和分成 $r=1,2$ 和其余项；其余 $q_r$ 与 $V_{12}p_1p_2$ 对易。第三行在恰有 $j$ 个激发的对称扇区上使用
$\langle q_1q_2\rangle=j(j-1)/[K(K-1)]$，再乘 $K/j$。因此第2项不超过
$4\sqrt2\,w\sqrt{\alpha(\alpha+2/K)}
\le4\sqrt2\,w(\alpha+1/K)$。

合并四项，取稍粗而方便的常数：

$$
|\dot\alpha|\le12w(\alpha+K^{-1}),\qquad
\alpha(0)=0,\qquad
\alpha(t)+K^{-1}\le\frac{e^{12wt}}K .
\tag{10}
$$

$K=1$ 时直接用 $|\dot\alpha|\le2w$ 亦满足同一粗界。$w=0$ 则 $\alpha=0$。

式(10)控制的是激发**比例**。允许背景内有有限数量或更多但占比受控的关联激发；没有断言 $\Phi_K$ 在整体向量范数中趋于Hartree乘积。

## 3. 直接控制平均算符涨落，避免损失四次根

固定某个Hermitian $B$，记 $\beta=\|B\|$、
$b=\langle\phi,B\phi\rangle$、$C=B-bI$。因为 $pCp=0$，

$$
\overline B-bI=S+S^\dagger+D,\qquad
S=\frac1K\sum_r p_rB_rq_r,\qquad
D=\frac1K\sum_rq_rC_rq_r.
\tag{11}
$$

三角不等式和式(6)给出
$\|S\Phi_K\|\le\beta\sqrt\alpha$、
$\|D\Phi_K\|\le2\beta\sqrt\alpha$。
不同因子的项对易，而同因子的交换子满足

$$
[S,S^\dagger]
=K^{-2}\sum_r(p_rB_rq_rB_rp_r-q_rB_rp_rB_rq_r)
\le\frac{\beta^2}{K}I.
\tag{12}
$$

所以
$\|S^\dagger\Phi_K\|^2
\le\beta^2(\alpha+K^{-1})$，最终得到

$$
\begin{aligned}
\|(\overline B-b)\Phi_K(t)\|
&\le\beta\bigl(3\sqrt{\alpha(t)}
+\sqrt{\alpha(t)+K^{-1}}\bigr)\\
&\le4\beta\sqrt{\alpha(t)+K^{-1}}
\le\frac{4\beta e^{6wt}}{\sqrt K}.
\end{aligned}
\tag{13}
$$

这里无需先把一体／二体约化态误差变成trace-norm界再开平方。因此没有人为降为 $K^{-1/4}$。

## 4. 保留实际背景纠缠和探针反作用的完整误差

定义固定背景准备等距
$E_K\chi=\chi\otimes\phi_0^{\otimes K}$。
以实际传播子作Duhamel左侧传播，比较态在每个中间时刻仍为
$U_P^{\rm mf}(t)\chi\otimes\Phi_K(t)$。因此

$$
\begin{aligned}
\|(U_{\rm act}(T)-U_{\rm cmp}(T))E_K\|_{\rm op}
&\le \sum_a\|A_a\|\int_0^T
\|(\overline B_a-b_a(t))\Phi_K(t)\|\,dt\\
&\le\frac{4L}{\sqrt K}\int_0^T e^{6wt}\,dt\\
&\le\frac{4LT e^{6wT}}{\sqrt K}
=:\epsilon_K(T),\\
L&=\sum_a\|A_a\|\,\|B_a\|.
\end{aligned}
\tag{14}
$$

$w>0$ 时中间积分可写成 $(e^{6wT}-1)/(6w)$；$w=0$ 时为 $T$。当然可将最终向量界截为2。

式(14)是完整输出的准备等距范数界。它对任意未知探针及任意维参考 $R$ 同时成立，因为残差估计只用 $\|A_a\|$，没有对探针态采样。背景、源、全部控制和记录留在输出中；实际的探针反作用没有删除，只是在这个复制尺度、指定准备和有限窗内，与比较过程的差被界住。

任意相同后续完整仪器下，联合结果总变差距离不超过 $\epsilon_K(T)$，通道diamond距离不超过 $2\epsilon_K(T)$。这些是非选择、包含失败记录的界；稀有后选择条件概率没有同样的无分母保证。若实际源／端点／钟仪器还有独立实现误差，应另相加，不能由均场估计替其取得。

## 5. 自作用计数修正不是被遗漏的整体相位

式(2)没有 $r=s$ 项，而式(3)使用没有 $(1-1/K)$ 系数的Hartree方程。这一差别已由第2节的 $2w/K$ 项明确处理。

若某具体完整模型另有
$K^{-1}\sum_r D_r$，其中 $D=D^\dagger$、$\|D\|\le d$，最直接的办法是把它完整并入 $h_K=h+D/K$，背景比较仍保留新模型的精确演化。相应Hartree解记为 $\phi_K$。与式(3)的解比较，由
$\|w^\phi-w^\psi\|\le2w\|\phi-\psi\|$ 可得

$$
\|\phi_K(t)-\phi(t)\|
\le\frac{dt\,e^{2wt}}K,\qquad
|b_{a,K}(t)-b_a(t)|
\le\frac{2d\|B_a\|t\,e^{2wt}}K.
\tag{15}
$$

若仍用未修正的探针比较器，式(14)另加
$LdT^2e^{2wT}/K$ 即可。这里的 $O(K^{-1})$ 是探针有效系数的修正；不能据此宣称整个背景向量只变了 $O(K^{-1})$。整个背景哈密顿量修正范数可为 $O(1)$，精确背景因子必须继续保留。

## 6. 对下一模型的真正作用与停止条件

1. 这是有限复制模型的条件误差工具。$K$ 份物质、关系体和初始准备具有实际资源代价；总能源通常为 $O(K)$，没有被删除。有限自伴 $H_{\rm act}$ 可用常数统一移到正谱，但这一选择不等于取得实际制备和读出的免费权限。
2. 空间管、cluster维数、$\|W\|$ 及有限时间 $T$ 必须先固定。本文不给无限空间体积统一界，也不给 $T$ 随 $K$ 任意增长的保证。若 $W$ 的范数因空间体积增长而增大，须公开相应常数。
3. 同站复制间有 $O(K)$ 个弱耦合邻项；每个复制的耦合强度总和仍为 $O(1)$。这可以保持原空间邻接支持，却不是“每节点固定度数”的证明。任何有限传播估计须由具体空间Hamiltonian另行核验。
4. 背景只对式(1)及明确可推广的相近准备取得结论；任意未知背景、跨复制猫态或任意与参考纠缠的背景不受此定理保证。一个全同方向猫态可能有 $O(1)$ 的平均算符涨落，不能被大 $K$ 自动经典化。
5. $\phi(t)$ 可保留完整单cluster内部物质—关系关联，$b_a(t)$ 可由同一来源和关系体的自主Hartree方程产生。证明不要求手工外给 $b_a(t)$，但也不证明其等于某种Lorentz几何系数。后者需要原钟、端点、信号的共同运输和未拟合留出。
6. 探针数量和 $L$ 固定时，反作用对背景平均量可随 $K$ 受抑；这不是有限源反作用被忽略的许可。若下一任务要求 $O(1)$ 的写回或宏观改变，必须把相应源纳入cluster／改变缩放，再重新证明，不能套本界宣称既可忽略又有宏观响应。
7. 即使式(14)很保守而导致所需 $K$ 很大，它仍给有限、可审核的条件预算。下一轮的净增应是同一具体空间过程的原资料、钟、传播和留出；不能仅把本页成熟均场工具再编号为物理生成成果。

本页没有进行数值实验，也没有修改任何冻结轮次。其交付是可逐行检查的有限维证明及适用范围。
