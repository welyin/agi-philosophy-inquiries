# 第十七轮：完整单轴记录等价于一个随机翻转的隐藏比特

日期：2026-09-15。接续第十六轮，将“仍有经典解释”加强为明确、可复算的二状态等价。

## 1 任务与边界

固定读取方向 0。比较的是任意有限二元记录串的概率，以及由这些记录更新的 u；不是只匹配某两个均值，也不是拟合一条样本路径。

原模型含连续潜变量和三摘要。替代实现只保留一个隐藏比特 X=±1。这个比特是另一种生成记录的模型，不被解释为已经识别了原来的潜变量 x。

## 2 离散步的精确二状态实现

第 n 个实验取隐藏比特的初始概率

\[
p_0=\tfrac12(1+u_0,1-u_0)^T.
\]

每步先根据当前比特产生记录，再独立翻转比特：

\[
\Pr(s|X)=\frac{1+s\eta_nX}{2},\qquad
r_n=\Pr(X\text{ 翻转})=\frac{1-\alpha_n}{2}.
\]

这里 r_n 是每步概率，后面的 r 才是连续率。写

\[
T_n=\begin{pmatrix}1-r_n&r_n\\r_n&1-r_n\end{pmatrix},\quad
D_{s,n}=\tfrac12\operatorname{diag}(1+s\eta_n,1-s\eta_n),
\quad H_{s,n}=T_nD_{s,n}.
\]

用 (w,u)=(p_++p_-,p_+−p_-) 换基，H_{s,n} 恰好变成

\[
\frac12\begin{pmatrix}1&s\eta_n\\
\alpha_n s\eta_n&\alpha_n\end{pmatrix},
\]

即原 A_{s,n} 的整个相关子块。因此对**每个**记录串都有

\[
\mathbf1^T H_{s_k,n}\cdots H_{s_1,n}p_0
=e_0^T A_{s_k,n}\cdots A_{s_1,n}z_0.
\]

这是一项逐矩阵相同的解析结论，不受枚举长度限制。程序还分别枚举记录串与隐藏比特路径，独立复核这个等价。

根据一个结果先做 Bayes 更新再做隐藏翻转，条件均值为

\[
u'=\alpha_n\frac{u+s\eta_n}{1+s\eta_nu},
\]

也与原条件预测精确相同。因而它是整个固定方向记录模型的替代实现。

## 3 连续极限：电报过程的带噪观察

因为 nr_n→ν/12，定义连续隐藏比特以对称率

\[
r=\nu/12=a/2
\]

翻转，其生成矩阵为 Q=[[-r,r],[r,−r]]。注意比特均值衰减率是 2r=a，不能把翻转率与均值衰减率混为一谈。

观察定义为

\[
dY_t=\sigma X_tdt+dB_t,
\]

其中 B 与隐藏比特过程独立。令 U_t=E[X_t|Y_{[0,t]}]，则二状态 Wonham 滤波器给出

\[
dU_t=-2rU_tdt+\sigma(1-U_t^2)
(dY_t-\sigma U_tdt).
\]

这正是第十五轮 U、Y 的系统。核对来源为 [van Handel 博士论文，推论 1.2.1](https://minty2.stanford.edu/wp/wp-content/thesis/RVHandel_thesis.pdf)（印刷页 27，PDF 第 41 页），以及 [Chigansky、Liptser、van Handel，§2.1、式 (2.13)–(2.14)](https://web.math.princeton.edu/~rvan/handbook-new.pdf)。本轮读取了对应有限状态信号、独立布朗观察和滤波方程的正文，未声称复现其一般稳定性理论。

取相同初始均值，由第十五轮的唯一性，两个连续模型的 (U,Y) 路径分布相同。潜在过程本身可以不同；观察等价不意味着潜变量历史相同。

## 4 不只比较矩：得到终端记录分布的特征函数

特征函数 χ_T(ω)=E[e^{iωY_T}] 编码整个终端记录分布。对隐藏模型，带观察相位的 2×2 生成矩阵为

\[
K(\omega)=Q+i\omega\sigma\operatorname{diag}(1,-1)
-\frac{\omega^2}{2}I.
\]

因此 χ_T(ω)=1ᵀe^{TK(ω)}p₀。令 d=√(r²−σ²ω²)，直接计算 2×2 指数得到

\[
\boxed{\chi_T(\omega)=e^{-(r+\omega^2/2)T}
\left[\cosh(dT)+(r+i\omega\sigma u_0)\frac{\sinh(dT)}d\right].}
\]

d 可以为纯虚数；d=0 时 sinh(dT)/d 取 T。这里的复数仅用于经典特征函数，不是引入量子概率幅。

独立的离散复算使用

\[
M_n(\omega)=\sum_s e^{i\omega s/\sqrt n}H_{s,n},\qquad
\chi_{n,T}(\omega)=\mathbf1^TM_n(\omega)^{nT}p_0,
\]

其中 nT 为整数。展开 M_n=I+K(ω)/n+O(n^{-2})，便得上述连续表达式。对多个不相交区段使用不同 ω 并连乘，就得到联合增量特征函数的同类极限；终端公式本身不应被误称为全部路径统计。

ν=0 时比特不翻转，公式退化为

\[
\chi_T(\omega)=e^{-\omega^2T/2}
[\cos(\omega\sigma T)+iu_0\sin(\omega\sigma T)],
\]

即均值 ±σT、方差 T 的两个高斯分布的经典混合。对 ν>0，则是积分电报过程加独立高斯噪声。脚本用离散记录矩阵验证多个频率的收敛，也覆盖 d=0 的退化点。

## 5 二状态解释何处停止

这个比特只承载 u。不同 v 的原准备在固定轴上本来就不可分辨，因而舍去 v 不影响该协议类别。

但旋转探针会读出 u cosθ+v sinθ。对 u 相同、v 不同的输入，本轮隐藏比特初态完全相同，因此它无法同时复现这些旋转后的结果。不能把本轮结论扩大为“整个三摘要理论只有两个经典状态”，更不能由二状态实现失效反推多方向统计就是非经典的；原连续潜变量核仍在。

## 6 研究判断与下一轮入口

至此，固定方向的弱记录路线已被一个明确的经典隐藏比特模型完整解释。非线性条件更新、同一创新噪声、边界保持和可识别的噪声率，都不足以独立产生量子内容。

下一轮应扩大实际可检验的协议：在同一个已有经典核内加入预先声明的不同方向读取，检验二状态压缩在哪些记录上失败，并求这些协议所需的最小预测表示。进一步检验若保留所有连续可逆平移，有限经典隐藏状态是否还能精确实现。重点区分“二状态不够”“有限状态不够”与“任何经典解释都不可能”；只有最后一种才涉及更强的非经典结论。

代码：[two_state_filter.py](two_state_filter.py)。结果：[two_state_filter_results.json](two_state_filter_results.json)。新增 10 项检查通过，无新增依赖。

```text
python -X utf8 research_cognition_physics/two_state_filter.py --write-results
```
