# 正能初始约束面上的局部接续必要条件

## 1. 要减少的是哪个输入

351在整个未约束相空间要求强HDA身份，得到Hamiltonian迹系数λ=1/2。但那里具体见证不满足Hamiltonian约束，不能直接说原实际系统能准备它。本轮给一份同时满足两类初始约束的正源数据，将必要性缩到这份数据的约束保持。一般括号、λ-R与CMC机制均为旧工具，不是本轮新发现。

认知动机是：许可的共同初始关系，在同一接续中不应立即失去其资格。**哪些初始关系实际许可，仍须新增准备合同证明；不是由“认知”一词推出。**

## 2. 完整候选及约定

在周期2π的三维环面上，积分和正则配对均取同一归一化体积。沿用351泊松方向，时间导数为 $\dot F=\{F,H[N]\}$。取Λ=0，完整密度为

$$
\mathcal H_\lambda=
\frac{\pi^{ij}\pi_{ij}-\lambda\pi^2}{\sqrt q}-\sqrt q R
+\frac{p_\phi^2}{2\sqrt q}
+\frac{\sqrt q}{2}q^{ij}\partial_i\phi\partial_j\phi,
\qquad \pi=q_{ij}\pi^{ij},
$$

$$
\mathcal D_i=-2q_{ij}\nabla_k\pi^{jk}+p_\phi\partial_i\phi.
$$

标量有正常正动能；其度规变分必须进入几何动量方程。这是经典约束生成元，不是正量子Hamiltonian，不能靠常数谱移将两者等同。

直接复用351、352，有完整非消去异常

$$
\{H[N],H[M]\}=D[v]+2(1-2\lambda)
\int v^i\nabla_i\pi,
\qquad v^i=q^{ij}(N\partial_jM-M\partial_jN).
\tag{1}
$$

物质自身给 $D_m[v]$。两个引力—物质交叉项的局部lapse因子为NM，反对称相消；不是删除物质反作用。无导数势不改变这一机制，本轮为最短见证固定为零。

## 3. 同时满足两种初始约束的正源族

令 $z=\sin x^1$，$\lambda\in[2/5,3/5]$，规定

$$
q_{ij}=\delta_{ij},\qquad
\pi^{ij}=2\delta^{ij}+z\,\delta^{i2}\delta^{j2},\qquad
\phi=0,
\qquad p_\phi=\sqrt{-2h_\lambda(z)},
\tag{2}
$$

其中

$$
h_\lambda(z)=12-36\lambda+(4-12\lambda)z+(1-\lambda)z^2.
\tag{3}
$$

因为 $\partial_\lambda h=-(6+z)^2<0$，上界取λ=2/5的凸二次式端点，下界取λ=3/5的单调二次式端点，精确得到

$$
-\frac{62}{5}\le h_\lambda\le-1,
\qquad
1\le\rho_\phi=\frac{p_\phi^2}{2}\le\frac{62}{5},
\qquad
\sqrt2\le p_\phi\le\sqrt{\frac{124}{5}}.
\tag{4}
$$

所以源平滑、严格正、环面总量有限。$R=0$，故逐点 $\mathcal H_\lambda=0$。$\pi^{22}$ 沿x¹变化而非x²变化，故 $\nabla_j\pi^{ij}=0$；$\phi$ 常数，故 $\mathcal D_i=0$。但迹 $\pi=6+\sin x^1$ 不恒定。

不同λ使用不同的 $p_\phi$，不能把它们叫同一数值初态。该族证实每个候选最初两条约束都允许这类正源准备；没有声称加入全部次级／三级约束后它还许可。

## 4. 一份正lapse即可给必要性，所有正lapse也救不了该初值

取 $N=1,M=2+\sin x^1$，两者都严格正。$v^1=\cos x^1$，$D[v]=0$，归一化 $\int\cos^2x^1=1/2$。因此

$$
\{H[1],H[2+\sin x^1]\}=1-2\lambda,
\qquad
\left.\frac{dH[2+\sin x^1]}{dt}\right|_0=2\lambda-1.
\tag{5}
$$

若这份初值必须沿同一正unit-lapse生成元保持初始约束，必要条件就是λ=1/2。这里不需要错误λ的长期解存在；任何声称可微且保持约束的曲线在初刻已不相容。

更强地，对任意光滑严格正N，在该初始约束面式(1)给

$$
\dot{\mathcal H}
=2(1-2\lambda)(N\Delta\pi+2\nabla N\cdot\nabla\pi).
\tag{6}
$$

当λ≠1/2，保持全部局部H要求右侧为零。乘N后是

$$
\nabla\cdot(N^2\nabla\pi)=0.
$$

环面无边界，乘π积分得

$$
0=-\int N^2|\nabla\pi|^2.
\tag{7}
$$

而本数据 $\nabla\pi=\cos x^1e_1$，正连续N有正下界n₀，积分至少 $n_0^2/2>0$，矛盾。平滑shift只对已经为零的约束作Lie运输，不能抵消该初刻项。**改变正lapse不能救这份初值；可以救候选理论的方式是拒绝这份准备或改变约束／模型族。**

只含整体Hamiltonian约束的projectable模型不是这里要求全部局部H的合同；不能把同一常N误当成它已经满足本地约束传播。

## 5. 不靠异常公式的完整切向核对

初片令矩阵 $\Pi=\operatorname{diag}(2,2+z,2)$。unit lapse和零shift下，完整场方程给

$$
\dot q=2(\Pi-\lambda\operatorname{tr}\Pi I),\qquad
\dot\Pi=-2\Pi^2+2\lambda(\operatorname{tr}\Pi)\Pi,
\qquad \dot\phi=p_\phi,\quad\dot p_\phi=0.
\tag{8}
$$

第二式中，几何行列式变分的 $h_\lambda I/2$ 与源反馈 $p_\phi^2 I/4$ 相消。若省略源的q变分，此式错误。

动能加源的局部代数项沿式(8)正好抵消，剩下曲率的切向变化。在平直初片

$$
\dot R=\partial_i\partial_j\dot q_{ij}-\Delta\operatorname{tr}\dot q
=2(1-2\lambda)\sin x^1,
\qquad \dot{\mathcal H}=-\dot R,
\qquad \dot{\mathcal D}_i=0.
\tag{9}
$$

代码把 $q(t)=I+t\dot q,\Pi(t)=\Pi+t\dot\Pi,\phi(t)=tp_\phi$ 代入完整非线性密度，包括精确对角空间曲率、源梯度与行列式，再中心差分。**该直线是正确切向的诊断曲线，不是解，更不是有限时间物理预测。** 两种步长核对式(9)与式(5)，不以差分替代解析证明。

初片 $\dot\phi=Np_\phi>0$ 给一个标量计时角色的局部候选；存在读数单调性不等于实际量子钟、仪器或独立 $T_{21}$ 已建立。

## 6. 删弱与真实限制

1. 令非均匀振幅为零，即 $\Pi=2I$，仍有正源且两初始约束满足，但迹梯度与本次异常消失。这只显示检测的盲区，不证明任意λ完整理论已健康。
2. λ≠1/2的理论可增加CMC等条件，排除式(2)。经典λ-R文献早已说明额外约束和lapse地位；这不违反本轮定理。本轮证明的命题是“**额外许可式(2)的准备，并要求原完整生成元保持两类局部约束**”在本类中的必要系数，不是所有认知原则已选择该系数。
3. 每个非正确λ都有非零精确切向残差，不等于在固定实际误差下可区分任意近邻λ。有限回路须另证高阶余项和读出误差；本轮没有这份证书。
4. 尚无全量子物理码、未知输入/R保持、正能源父过程、实际钟尺或信号的共同实现。源只是选定经典标量，不是已经接回1178完整材料。G、Λ、其他模式、高阶项和暗部门均未被选定或删除。

原始文献范围核对：[Loll–Pires §2](https://arxiv.org/html/1407.1259)、[Bellorín–Restuccia](https://arxiv.org/html/1004.0055)。其λ通常是Lagrangian系数，应经351式(12)转换；本轮不声称发现已知CMC机制。

## 7. 下一项

不再扫描λ或扩大网格。优先判式(2)及其有限邻域能否来自同一量子物质—几何父过程的可许可准备，且保原未知资料/R。边界运输须由该过程自己的钟、端点和信号独立规定；仅在普通模拟器中重现式(5)不会补上这项桥。若实际准备只能是CMC，则诚实保留本次条件选择不适用；不能从认知名义禁止它。
