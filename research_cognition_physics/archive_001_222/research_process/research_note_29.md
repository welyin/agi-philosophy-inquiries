# 第二十九轮：把完整更新核的存在性化为凸序

日期：2026-09-15。接续第二十八轮提出的联合分布问题。固定原准备半宽 h=1/4、α=sinc(h)≈0.9896158370，研究第七轮 λ=1 的 fresh 族，0≤η<1。

## 1 结果与范围

对固定的仿射隐藏角编码

\[
F_z(x)=\frac{w+(u/\alpha)\cos x+(v/\alpha)\sin x}{2\pi},
\qquad \sqrt{u^2+v^2}\le\alpha w,
\]

完整正更新核存在，当且仅当一个明确的内圆分布小于均匀单位圆分布的**凸序**。凸序指所有凸函数在前者上的平均不超过后者；它比仅比较协方差更强。

本轮先证明这一等价，再用一个绝对值函数取得必要条件

\[
\boxed{\alpha\big[\sqrt{1-\eta^2}+\eta\arcsin\eta\big]\le1.}
\]

其等号正根约为 0.144739232095。此处排除仅针对固定 F 编码；推广到其他准备非情境编码需要独立证明，见接续的第三十轮。

## 2 从核到条件均值

记 n_x=(cos x,sin x)，先令读取方向 θ=0。对结果 s=±1，定义

\[
q=\sqrt{1-\eta^2},\quad
f_s(y)=\frac{1+s\eta\cos y}{2},\quad
n_{\psi_s(y)}=\frac{(\cos y+s\eta,\ q\sin y)}{1+s\eta\cos y}.
\]

原 fresh 矩阵为

\[
A_sz=\frac12(w+s\eta u,\ \alpha u+s\alpha\eta w,\ \alpha qv).
\]

直接展开整个输出密度，得到

\[
F_{A_sz}(y)=\frac{f_s(y)}{2\pi}
[w+u\cos\psi_s(y)+v\sin\psi_s(y)].
\]

设输入为均匀隐藏角，构造联合测度

\[
J_s(dx,dy)=\frac{dx}{2\pi}K_s(dy\mid x).
\]

比较 w、u、v 的系数，核归一化及 K_sF_z=F_{A_sz} 等价于：

1. 输入 x 的边缘为 dx/(2π)。
2. 输出 (s,y) 的边缘为 f_s(y)dy/(2π)。
3. 给定 (s,y)，有 E[n_x|s,y]=αn_{ψ_s(y)}。

这是对完整概率测度的条件，没有把高阶谐波忽略掉。

## 3 输出标签变为内圆分布

第二十六轮的换元计算给出：上述输出边缘经 (s,y)↦ψ_s(y) 推送后，角密度为

\[
R_\eta(\psi)=\frac{q^3}{2}
\left[(1-\eta\cos\psi)^{-2}+(1+\eta\cos\psi)^{-2}\right].
\]

R 平均为 1，且 π 周期，所以相应向量均值为零。令

\[
X=n_x\sim\sigma,\qquad
Y=\alpha n_\psi\sim\nu_\eta,\qquad
\nu_\eta(d\psi)=R_\eta(\psi)\frac{d\psi}{2\pi},
\]

其中 σ 是均匀单位圆分布。上一节条件要求存在 E[X|Y]=Y 的耦合。

反方向也成立。若有该耦合，先给定 Y 抽取 X，再按既定输出边缘的条件分布抽取 (s,y)，两者在给定 Y 后独立。由于 α>0，Y 唯一确定圆上的 ψ。更明确地，令

\[
b_s(\psi)=\frac{q^3}{2(1-s\eta\cos\psi)^2},\qquad
T_s(\psi)=\operatorname{atan2}(q\sin\psi,\cos\psi-s\eta),
\]

则以 b_s(ψ)/R(ψ) 选择 s，再令 y=T_s(ψ)。这恢复指定的 (s,y) 边缘及条件均值。最后给定 x 分解联合测度，得到正核 K_s。圆和有限结果集合都是标准 Borel 空间，因此这里的条件概率分解适用。

由 Strassen 鞅耦合定理，

\[
\boxed{\exists\,K_s\text{ 完整保持 }F_z
\iff\exists\,\operatorname{Law}(Y,X):E[X|Y]=Y
\iff\nu_\eta\le_{\rm cx}\sigma.}
\]

使用的是 Leskelä–Vihola 对该定理的陈述，定理 1.1、1.2；两边支撑有界，有限一阶矩条件满足。若一个强度区间内逐点满足凸序，该文定理 1.3 还给出随参数可测的耦合选择。其他读取方向可直接旋转构造。[原文：Conditional convex orders and measurable martingale couplings](https://arxiv.org/pdf/1404.0999)。

这个标准定理提供存在性接口，没有自动判断本模型何时满足全部凸函数不等式。

## 4 一个比二次矩更强的必要条件

取 φ(v)=|v_x|，利用输出换元：

\[
\begin{aligned}
E_{\nu_\eta}|Y_x|
&=\alpha\int\sum_s f_s(y)|\cos\psi_s(y)|\frac{dy}{2\pi}\\
&=\alpha\int\max\{|\cos y|,\eta\}\frac{dy}{2\pi}\\
&=\frac{2\alpha}{\pi}G(\eta),\\
G(\eta)&=\sqrt{1-\eta^2}+\eta\arcsin\eta.
\end{aligned}
\]

单位圆的对应平均为 2/π，于是 αG(η)≤1。G(0)=1，G′(η)=arcsin η>0，因此等号在 (0,1) 中有唯一正根 η_B。

数值定位 η_B≈0.1447392320951458。结果文件中的浮点夹值仅用于复算数值，不能当作最后一位的严格区间证书；第三十轮另给 η=0.145 的严格有理数排除。

作为对照，二次函数 φ(v)=v_x² 只要求

\[
\frac{\alpha^2}{2}(1+c_1)\le\frac12,\quad
c_1=t^2(1+2q),\quad t=\frac{\eta}{1+q}.
\]

在 η=0.15，它给出约 0.49796418≤0.5，仍然通过；绝对值测试却给出约 0.63710999>2/π，已经排除固定编码的核。因此仅解低阶矩的正性问题仍会漏掉障碍。

## 5 对接后续研究

现在可以从两侧推进：构造保持条件均值的正耦合，或寻找违反凸序的凸函数。有限扫描未找到违反，不能证明全部凸序成立；构造中的某个权重变负，也不能证明所有耦合不存在。

这轮没有把一般准备非情境模型强行替换为单位圆编码。一个一般编码可有模长小于 1 的仿射系数，且带有额外隐藏变量。下一轮先利用精确旋转处理这一差别，再改进正耦合。

## 6 复算与文献边界

代码：[convex_order_interface.py](convex_order_interface.py)。结果：[convex_order_interface_results.json](convex_order_interface_results.json)。7 项检查通过，包含边缘归一化、输出换元、独立角积分、绝对值与二次矩对照、已知可行点及根的数值定位。

~~~
python -X utf8 research_cognition_physics/convex_order_interface.py --write-results
~~~

本轮实际核对上述论文前言中的凸序定义、定理 1.1—1.3 及适用条件；没有宣称复现全文证明或提出新的 Strassen 定理。本模型到耦合的双向翻译及具体必要函数计算写在本笔记中。准备非情境性仍是新增假设，尚未从认知原则导出。
