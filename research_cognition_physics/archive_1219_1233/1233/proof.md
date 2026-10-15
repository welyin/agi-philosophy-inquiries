# 1233：同Dirac自治记录过程的区域能源、边界通量与控制源

## 0. 净增与命题范围

本轮复用[1230](../1230/proof.md)的完整两味质量1 Dirac真空、同两模式、两只正能Gaussian程序钟和全部P/E寄存器。没有更换参考、删除Dirac海、改变门许可或把程序钟当加速固有时钟。结论是：**在原实际自治过程的每个时刻，正常序区域能源期望、独立边界通量和真实控制源存在，并满足同一弱收支；收入相应钟与交互存能后，源功抵消。** 全空间极限回到原完整能账。

一般守恒身份不是新发现。[1128 §7](../../archive_1120_1158/1128/proof.md)已有格点自治局部能源，[1138](../../archive_1120_1158/1138/proof.md)已有时间涂抹能源，[1150 §§3—6](../../archive_1120_1158/1150/energy_balance.md)已有另一KG有效过程的实际中插弱收支。当前增量是把该接口接到**1228—1230这一完整连续Dirac海、正常初态、实际双钟和实体记录的同一全过程**，解决尖时局部能源未必作用于真空的域障碍。[1228 §4](../1228/proof.md)的终态后自由源延拓和[1229 §6](../1229/proof.md)的终态boost能源不能直接替代这里的实际驱动历史。

本轮不授予区域能源测量仪器，不证明逐点门实现、装置的空间四动量、面积动力学、Einstein方程或阶段完成。遵守[阶段范围审计](../_admission/stage_scope_audit_after1230.md)。

## 1. 原材料、过程和边界先于能源

用同一单位与惯性截面，单粒子空间为两味 $L^2(\mathbb R^3,\mathbb C^4)$，

$$
h=\boldsymbol\alpha\cdot p+\beta,\qquad p=-i\nabla,
\qquad h^2=-\Delta+1,\qquad S=P_-(h).
\tag{1}
$$

真空为 $\Omega$，正常序自由能源 $H_F\ge0$ 且 $H_F\Omega=0$。1230的 $f_1=f,f_2=g$ 是同味、互相正交、归一的 $C^7\cap H^8$ 紧支包；两个菜单均允许横半宽、纵厚与中心距的 $10^{-8}$ 原误差。它们满足 $\|hf_j\|<400$，没有截去负谱。

实际Hamiltonian完全保留为

$$
H=H_F+H_{C1}+H_{C2}+4+V_1+V_2,
\quad H_{Cj}=P_j^2/(2M)+1,
\quad V_j=g_\delta(X_j)n_{f_j}K_j,
\quad n_u=c^\dagger(u)c(u),
\tag{2}
$$

$$
K_j=\frac\pi2(I_\ell-X_\ell)_{P_jE_j},\qquad
g_\delta(x)=\frac{15}{16\delta}
\left(1-\frac{(x-1)^2}{\delta^2}\right)_+^2,
\quad 0\le K_j\le\pi,
\quad\|g_\delta\|=\frac{15}{16\delta}.
\tag{3}
$$

原参数为 $M=10^{260},\sigma=10^{-120},\delta=10^{-24},T=2$；两钟初态均为均值 $(X,P)=(0,M)$ 的正常Gaussian，全部四个P/E寄存器初态保留。两门中心同时为1，场模式始终是该菜单的固定Schrödinger模式。完整自由初态与准备成本沿1230；有界正交互给自伴、半有界的同一 $H$，所有以下期望均取真实 $e^{-itH}\Psi_0$，$0\le t\le T$。

先选择实 $\chi\in C_c^2(\mathbb R^3)$，再定义区域；不从熵反推边界。取

$$
A_\chi=\tfrac12(\chi h+h\chi)
=\chi h-\tfrac i2\boldsymbol\alpha\cdot\nabla\chi,
\qquad
J_\chi=i[h,A_\chi]
=\tfrac i2[h^2,\chi]
=\tfrac12\{\nabla\chi,p\}
=\nabla\chi\cdot p-\tfrac i2\Delta\chi.
\tag{4}
$$

$A_\chi$ 是能源权重，$J_\chi$ 是指向区域内的边界流权重；两者是对称的 $H^1\to L^2$ 映射。若 $c_0=\|\chi\|_\infty,c_1=\|\nabla\chi\|_\infty,c_2=\|\Delta\chi\|_\infty$，则

$$
\|A_\chi u\|\le c_0\|hu\|+\tfrac12c_1\|u\|,
\qquad
\|J_\chi u\|\le c_1\|hu\|+\tfrac12c_2\|u\|.
\tag{5}
$$

式(4)来自实际Dirac方程，而不是先命名一个“流”再强加连续性。它在光滑场上等价于正常序能源密度的局部连续方程；允许相差无散度改进项的点态流表示，本轮固定其空间测试为(4)。只给期望差定义，不假设尖时 $Q(A_\chi)\Omega$ 存在。

## 2. 有界CAR截断与真实历史的无截断极限

令 $\Pi_R=1_{|p|\le R}$，只用它证明定义，**不是材料的物理高能截断**。对 $B=A_\chi$ 或 $J_\chi$，令 $B_R=\Pi_RB\Pi_R$，并定义

$$
Q_R(B)=d\Gamma_{\rm CAR}(B_R)-\operatorname{Tr}(SB_R)I.
\tag{6}
$$

这是海表示中的正常序CAR二次元。它在每个有限R实际有界：对任意紧支可积乘法 $b$，$|b|^{1/2}\Pi_R$ 为Hilbert—Schmidt，故 $\Pi_R b\Pi_R$ 为trace class；再用 $h\Pi_R,p_k\Pi_R$ 有界，得到 $B_R$ 为trace class。二次元可用绝对收敛的核分解定义，其全Fock范数由trace norm控制。这里没有使用错误的“Hilbert—Schmidt核就给有界全Fock二次元”；其界可随R发散。真空均值严格为0。

对归一 $u$，CAR给精确有界交换子

$$
[Q_R(B),n_u]
=c^\dagger(B_Ru)c(u)-c^\dagger(u)c(B_Ru),
\qquad \|[Q_R(B),n_u]\|\le2\|B_Ru\|.
\tag{7}
$$

记 $H_0=H_F+H_{C1}+H_{C2}+4$，$U_I(t)=e^{itH_0}e^{-itH}$。其有界强连续驱动为

$$
V_I(s)=\sum_j g_\delta(X_j(s))n_{f_j(s)}K_j,
\quad f_j(s)=e^{ish}f_j,
\quad X_j(s)=X_j+sP_j/M.
\tag{8}
$$

不同时间的钟因子不需要对易；下述式子使用完整时间序列。固定观测时刻t，设 $B_R(t)=e^{ith}B_Re^{-ith}$。Duhamel给

$$
q_{B,R}(t)
=\langle U_I(t)^*Q_R(B;t)U_I(t)\rangle_0
=\int_0^t\!ds\,
\langle U_I(s)^* i[V_I(s),Q_R(B;t)]U_I(s)\rangle_0.
\tag{9}
$$

右侧只有有界CAR交换子；正常序真空项为0。由于 $B:H^1\to L^2$ 连续、$f_j\in H^8$，

$$
B_R(t)f_j(s)\longrightarrow B(t)f_j(s)
\quad\text{在 }(s,t)\in[0,T]^2\text{ 上一致于 }L^2.
\tag{10}
$$

证明是 $\Pi_R\to1$ 在 $H^1$ 强收敛，对紧时间参数的 $H^1$ 紧集一致，再用B的连续性和第二次 $L^2$ 截断。由(7)—(10)，$q_{B,R}$ 一致收敛到有限、连续的 $q_B(t)$，与具体合格截断序列无关。等价定义是直接把(9)中的交换子用(7)的 $B(t)f_j(s)$ 替代。**这是实际态上的正常序期望差，不是新授予所有正常态都可用的局部能源算符。**

令 $b_A=c_0\,400+c_1/2,b_J=c_1\,400+c_2/2$。相同界适用于自由输送模式，故

$$
|q_A(t)|\le4\pi t\delta^{-1}b_A,
\qquad |q_J(t)|\le4\pi t\delta^{-1}b_J.
\tag{11}
$$

没有从C²借取 $O(R^{-1})$ 速率；若要求特定截断速率需另给更高导数。本轮结论不使用该额外条件。准备有H⁸足够提供这里所需的H²强导数，不需全部高阶矩。

## 3. 同一实际源、区域收支与时间测试

因 $\Pi_R$ 与h和S对易，$\partial_t Q_R(A;t)=Q_R(J;t)$，且正常序导数常数为0。式(10)也适用于J，故可在(9)内取强导数并通过极限：

$$
\frac{d q_\chi}{dt}(t)=j_\chi(t)+s_\chi(t),
\qquad q_\chi=q_A,\quad j_\chi=q_J,
\tag{12}
$$

其中真实源是同一完整态上的

$$
s_\chi(t)=\left\langle\sum_j
i\,g_\delta(X_j)K_j
\big[c^\dagger(f_j)c(A_\chi f_j)
-c^\dagger(A_\chi f_j)c(f_j)\big]\right\rangle_t,
\quad |s_\chi(t)|\le4\pi\delta^{-1}b_A.
\tag{13}
$$

没有把 $\langle gK\,c^\dagger c\rangle_t$ 换成边缘期望的乘积；钟—场—指针关联、完整控制存能和反作用都保留。源算符来自(2)，不是外给classical函数。真实未来继续使用(2)，并非终态后另换自由Hamiltonian。

外向通量定义 $\Phi_\chi=-j_\chi$；它由独立边界测试 $J_\chi$ 而非能源差定义。对任意 $\zeta\in C_c^1((0,T))$，

$$
-\int_0^T\zeta'(t)q_\chi(t)dt
=-\int_0^T\zeta(t)\Phi_\chi(t)dt
+\int_0^T\zeta(t)s_\chi(t)dt.
\tag{14}
$$

这既给实际全历史，也允许有限时间分辨率的测试。删掉源一般破坏(14)，即使末指针概率仍被别的近似拟合。源非零的有限CAR删除检验只说明必要性，不替代无限场的上述证明。

## 4. 把真实控制收入各区域：源功抵消而边界流保留

一个明确、事前确定的边界采用C²平台函数

$$
\theta(u)=\begin{cases}1&u\le0,\\
1-10u^3+15u^4-6u^5&0<u<1,\\0&u\ge1.
\end{cases}
\quad
\chi_1=\theta\!\left(\frac{|x|-0.3}{0.2}\right)
\theta\!\left(\frac{|y|-0.3}{0.2}\right)
\theta\!\left(\frac{|z-2|-0.1}{0.1}\right).
\tag{15}
$$

它对1230所有允许几何误差都恒1于 $f_1$ 的完整支撑邻域，恒0于 $f_2$ 支撑邻域；$\chi_2$ 取名义z中心3的同样固定函数，平台宽度已经覆盖D的不确定度，无需额外精确标定。二者支撑不交。边界平滑厚度是有效区域的明确选择，不是微观最小长度。可取 $c_0=1,c_1=24,c_2=900$，由 $|\theta'|\le15/8,|\theta''|\le6$ 得到；(11)—(13)给明确有限常数。

由于h局域且 $f_j\in H^1$，有

$$
A_{\chi_1}f_1=hf_1,\qquad A_{\chi_1}f_2=0,
\qquad s_{\chi_1}=s_{H,1}:=\langle i[V_1,H_F]\rangle_t.
\tag{16}
$$

同理交换1、2。同期正交模式密度对易，因而 $[V_1,V_2]=0$；注意这不要求不同时间相互作用表象门对易。完整实际演化中

$$
\frac{d}{dt}\langle H_{Cj}+V_j\rangle_t
=\langle i[H_F,V_j]\rangle_t=-s_{H,j}(t).
\tag{17}
$$

初态属于 $D(H_0)$，有界V保 $D(H)=D(H_0)$，$[H_F,V_j]$ 为(7)给出的有界算符；$H_{Cj}$ 的能源存在。在场的平滑有限激发与钟的Schwartz向量组成的共同核上，$H$ 与 $H_{Cj}+V_j$ 的形式交换子延拓为该有界源；按二次型Ehrenfest身份传到 $D(H_0)$ 的真实解即得(17)。不分别假定 $H_{Cj}$ 和 $V_j$ 的强时间导数各自有界。定义区域能源变化包含相应控制

$$
\Delta\mathcal E_j(t)
=q_{\chi_j}(t)+\Delta\langle H_{Cj}+V_j\rangle_t.
\tag{18}
$$

P/E自由常数可分配至对应区域，虽无变化仍保留原初能账。式(12)、(16)、(17)立即给

$$
\frac{d\Delta\mathcal E_j}{dt}=-\Phi_{\chi_j}(t).
\tag{19}
$$

这不是把所有控制能源定位成真实时空应力：钟在此按其被授权执行的有限区域门归属，未构造钟的物理世界管或空间动量。它准确闭合当前有效模型内的能源角色；更强四P归属另列。部分切过模式的边界仍用(12)—(14)，不把区域双局域门误拆成点源。

取包含两门和全部有限传播halo的更大紧支平台 $\chi_L$，并令 $\chi_{\rm sh}=\chi_L-\chi_1-\chi_2\ge0$。例如中心 $(0,0,2.5)$、三个方向半宽3的平台、向外宽1的过渡足够覆盖 $T=2$ 及全部几何误差。线性给

$$
q_L=q_1+q_2+q_{\rm sh},\qquad
j_L=j_1+j_2+j_{\rm sh},\qquad
s_L=s_1+s_2,
\tag{20}
$$

壳的真实控制源为0，壳通量保留。对更一般嵌套区域同样用差函数，不能把两条不同边界上的通量各自删掉。有限传播与(9)说明此大平台已给 $q_L(t)=\Delta\langle H_F\rangle_t$，$j_L(t)=0$；不需要先取无限体积热极限。

## 5. 全能账、同场来源与剩余接口

直接取已定义的全局 $H_F$，(7)给

$$
\dot E_F(t)=\sum_j s_{H,j}(t),
\qquad
\Delta E_F+\Delta E_{C1}+\Delta E_{C2}
+\Delta\langle V_1+V_2\rangle=0.
\tag{21}
$$

式(21)保真实初交互存能，不把端点Gaussian尾设为0。1230已经给 $|\Delta\langle V\rangle|<10^{-60}$ 和完整准备、实体记录、钟熵、能源及涂抹源误差；本轮原封复用，未重做仪器优化。

现在同一个连续Dirac材料和实际双钟过程同时具备：1230的独立有限面积与实际可读互信息响应；当前的区域/壳能源、真实源和完整控制反作用；1228—1229在各自同材料菜单的参考流及模能源接口。不同菜单允许不同准备和读口，不把1230相隔模式当1228的boost模式，也不把Gauss分支面积公式移植过来。原有限区域门许可与共同参考/单位/截面合同仍为输入。

现有一般引力重建接口还需具体参考区域的平衡/非平衡前提、几何响应、面积项系数和引力耦合选择等。**(19)—(21)不推出 $\delta A\propto\delta E$、$S=sA$ 或Einstein动力学。** 这些属于须明确接入的剩余条件，不能从目前的有效区域能源闭合中免签；同样不要求本轮先完成所有未来引力问题。

## 6. 可复算检查的资格

[check.py](check.py)有四组：有限CAR的正常序/源/通量恒等式；保全部有限控制与指针的自治矩阵全历史及区域抵消；遗漏源或遗漏交互存能的删除反例；原连续模式、事前区域支持和有限常数/旧证据绑定。有限矩阵只是代数诊断，**不是Dirac海或Gaussian钟的截断模拟**。无限场的资格来自第2—4节证明；数值不验证其UV极限，亦不授予尖时能源读仪。

默认与 `--verify` 只读，`--write-results` 只允许首次生成。[results.json](results.json)保源文件哈希与范围。所有历史稿、正式目标、导航和已验收结果不改。
