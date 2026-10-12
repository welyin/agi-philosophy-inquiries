# 1150：实际中插后的整站能源弱收支与因果接续

2026-10-12。独立证明附件；采用[1150准入](admission.md)的同一有限有效过程。本文给解析证明，不另增数值组数；最终验收由主线和另一代理共同核对。

范围是原裸全16维味/自旋码、任意正常被动参考R、原T0、τ=T0/2处的实际整站CAR味旋转，以及正常KG完整后态。未访问的材料模式保原自由Hamiltonian和正常真空，全部仍在Hilbert空间中，不迹掉。本文不以旧/新输出误差控制能源，也不领取等时局部应力算符或能源仪器。

## 1. 同一过程、匹配项与源的正则性

以μ=1的同一正常KG场为中介，采用单位制ℏ=c=1。记

$$
N_i=\sum_s(b_{1s,i}^\dagger b_{1s,i}
+d_{1s,i}^\dagger d_{1s,i}),\qquad
h_i(x)=a_*\rho_{s,i}(x),
$$

$$
K_i=H_{D,i}+\lambda^2cI+\lambda^2\delta m\,s_1N_i,
\qquad
V_i=\lambda N_i\chi(h_i),
$$

$$
H=K_A+K_B+H_\chi+V_A+V_B .
\tag{1}
$$

$H_{D,i}$是完整正常CAR自由材料Hamiltonian；$c,\delta m,a_*$采用主证明中的同一匹配值。质量匹配在此是已声明的有限带算符，不冒称完整连续Dirac质量变分。未访问高模式上的自由能源仍保留，其真空是本任务的明确不变准备。

各$N_i$有界，彼此对易，并与$K_A,K_B$对易。源只访问保留模式。场算符由

$$
\chi(h)=\frac{a^\dagger(\omega^{-1/2}\widehat h)
+a(\omega^{-1/2}\widehat h)}{\sqrt2},
\qquad\omega(k)=\sqrt{|k|^2+\mu^2}
\tag{2}
$$

定义；傅里叶变换取酉归一。有限个$N_i$本征部门中，(1)是可平方配方的线性正常场源加完整自由材料。因此得到同一自伴、半有界过程；本文不以形式三场乘积替代该定义。

球内$\rho_s$是光滑径向函数，球面为零，零延拓属于紧支撑$H^1(\mathbb R^3)$，并非$C^\infty$。只使用它所具备的有限正则性，不继承旧g的全部高阶能源矩。

给定实际旋转

$$
U_A=\exp\!\left[-i\theta\sum_s
(b_{0s,A}^\dagger b_{1s,A}+b_{1s,A}^\dagger b_{0s,A})\right],
\qquad\theta=\pi/4,
\tag{3}
$$

它作用于整站A，保spin，不与$N_A$对易。演化为

$$
e^{-i(T_0-\tau)H}\,U_A\,e^{-i\tau H},
\qquad\tau=T_0/2 .
\tag{4}
$$

实际场在τ不重置。无中插任务与(4)是两个菜单，不要求(4)仍执行无中插CZ。

## 2. 实际完整态的有限相干历史表示

无脉冲的一段中，对固定材料数值$n_A,n_B$，场态保持正常相干态，位移满足

$$
\alpha(t)=e^{-i\omega t}\alpha(0)
-\frac{\lambda}{\sqrt2}\omega^{-3/2}
(1-e^{-i\omega t})\widehat{(n_Ah_A+n_Bh_B)} .
\tag{5}
$$

完整相位由同一线性源演化保留，不删去自相位或历史干涉项。旋转(3)把有限材料分支作有限线性组合。故(4)全过程由有限个材料向量、正常相干态及R向量构成相干和；同一末材料分支可以包含多条不同场历史，不能用经典混合替换。

$\mu>0$及$h_i\in L^2$给出每条历史位移属于$D(\omega^{1/2})$。例如新产生的一段有

$$
\|\omega^{1/2}\Delta\alpha\|
\le\sqrt2|\lambda|\,
\|\omega^{-1}\widehat{(n_Ah_A+n_Bh_B)}\|
\le\frac{\sqrt2|\lambda|}{\mu}
\|n_Ah_A+n_Bh_B\|_2 .
\tag{6}
$$

有限段叠加仍有有限界。它保证相应经典场位形属于H¹、经典动量属于L²，因而场能和能流的空间积分有限。该界不随R维数改变。

更明确地，规范化相干向量$|\alpha\rangle,|\beta\rangle$之间的正常序场多项式矩阵元，等于重叠$\langle\alpha|\beta\rangle$乘“产生部分取$\bar\alpha$、湮灭部分取$\beta$”的复经典场积。经典动量与空间梯度均在L²，所有能流积由Cauchy–Schwarz控制。

任意纯码—R态最多以16个码基向量展开；一次(3)仅有限扩展历史数。其R系数的内积有统一有限矩阵界。混态再由纯化得到同一结论。因此以下等式覆盖任意正常R及未知spin纠缠，不要求已知的材料输入。

## 3. 区域量仅作可控制矩阵元定义

取$f\in C_c^\infty(\mathbb R^3)$，$0\le f\le1$，满足

$$
f=1\quad\text{在A完整源支撑的邻域},\qquad
f=0\quad\text{在B完整源支撑的邻域}.
\tag{7}
$$

球半径1、中心距3，故可选择A外、B外的无源过渡壳。f把整个A材料收入区域，不切割A内部自由模式。

对第2节的实际历史定义

$$
\mathcal E_\chi(f;t)
=\frac12\int f(x)
\langle:\pi^2+|\nabla\chi|^2+\mu^2\chi^2:\rangle_t\,d^3x ,
$$

$$
\langle S(x)\rangle_t
=-\frac12\langle:\{\pi(x),\nabla\chi(x)\}:\rangle_t,\qquad
\Phi_f(t)=-\int\nabla f(x)\cdot\langle S(x)\rangle_t\,d^3x .
\tag{8}
$$

正常序始终针对同一自由正常真空，保相干历史的交叉矩阵元。$\Phi_f$符号对应外向通量。第2节证明上述期望有限，且在有限时间区间可积。

**这里没有断言等时空间涂抹的局部能源是作用于真空的算符。** 3+1维中的该算符域问题不能由有限期望绕过。本文只需要所声明相干历史的正常序双线性矩阵元与紧时间测试，它们有独立定义；没有新增能源仪器权限。

需要严格弱化时，可先用球内光滑紧支源逼近$h_i$于H¹。线性源解及其相干位移在任意固定有限时间区间按能源范数收敛，源相位也收敛。由L²积界，(8)及其时间测试收敛；下列光滑弱身份遂传到真实$h_i$。该步骤不需要源是$C^\infty$，也不改变原固定T0。

## 4. 两段自由演化中的精确弱收支

脉冲时刻以外，由实际(1)得到

$$
\dot\chi=\pi,\qquad
\dot\pi=(\Delta-\mu^2)\chi
-\lambda(N_Ah_A+N_Bh_B).
\tag{9}
$$

材料与场的等时对易关系由共同酉演化保持。局部场能的弱导数为

$$
\partial_t\langle e_\chi(x)\rangle_t
+\nabla\cdot\langle S(x)\rangle_t
=-\lambda\sum_i h_i(x)
\langle N_i\pi(x)\rangle_t .
\tag{10}
$$

等价地可用反对易乘积书写；等时$N_i$与$\pi$对易，二者一致。对f积分，使用(7)：

$$
\frac{d}{dt}\mathcal E_\chi(f;t)
=-\Phi_f(t)-\lambda\langle N_A\pi(h_A)\rangle_t .
\tag{11}
$$

另一方面，$[H,N_A]=[H,K_A]=0$，而$[V_B,V_A]=0$，所以

$$
\frac{d}{dt}\langle K_A+V_A\rangle_t
=+\lambda\langle N_A\pi(h_A)\rangle_t .
\tag{12}
$$

从同一H来的真实源功精确抵消。令

$$
\mathcal E_A(t)
=\langle K_A+V_A\rangle_t+\mathcal E_\chi(f;t),
$$

即有每段中的弱恒等式

$$
\frac{d}{dt}\mathcal E_A(t)=-\Phi_f(t).
\tag{13}
$$

这是材料、有限能源匹配、交互存能和传播流的同一账，不是单独自由壳流，也不是由终末任务接近推得无界量接近。共同常数c在固定支持下无时间导数，但保其区域归属；本文未从常数不可见性推出它对引力来源无贡献。

## 5. 实际旋转的功及全Hamiltonian身份

取τ前的真实完整态$\rho_-$，$\rho_+=U_A\rho_-U_A^\dagger$。KG整个边缘、B及未访问高模式在该瞬间不改变。定义

$$
W_A=\operatorname{Tr}\rho_-
\left[U_A^\dagger(K_A+V_A)U_A-(K_A+V_A)\right].
\tag{14}
$$

因其余完整H项与$U_A$对易，

$$
W_A=\operatorname{Tr}(H\rho_+)-\operatorname{Tr}(H\rho_-)
=\mathcal E_A(\tau+)-\mathcal E_A(\tau-).
\tag{15}
$$

本任务恰有每站一个粒子、无反粒子及高模式真空，这个部门被全部分段演化及(3)保持。其上

$$
K_A=k_0I+\nu_{\rm eff}P_{1,A},\qquad
N_A=P_{1,A},
\qquad
\nu_{\rm eff}=\nu+\lambda^2\delta m\,s_1 ,
$$

故实际工作量可写为

$$
W_A=\left\langle
(\nu_{\rm eff}+\lambda\chi(h_A))
\bigl(U_A^\dagger P_{1,A}U_A-P_{1,A}\bigr)
\right\rangle_- .
\tag{16}
$$

材料—场关联的交互功不可省略。令$\Delta P=U_A^\dagger P_1U_A-P_1$，对该实际态的Cauchy–Schwarz给

$$
|W_A|\le
\|\Delta P\|\left(
|\nu_{\rm eff}|+
|\lambda|\sqrt{\langle\chi(h_A)^2\rangle_-}
\right)<\infty .
\tag{17}
$$

此二阶矩由第2节及真空涂抹方差$\|\omega^{-1/2}\widehat h_A\|^2/2$有限而得，不靠迹距。对于输入00，τ前源为零、KG为真空，$\theta=\pi/4$给$W_A=\nu_{\rm eff}/2$；其正号随主证明的频差界核验。一般未知输入可以吸收或释放能量。

有效脉冲权限需要执行系统承担(14)的功，本文把责任算清，但没有自治制造脉冲装置，也没有闭合其空间反冲或全装置四P。

## 6. 时间测试包含真实工作跳变

对$\zeta\in C_c^\infty((0,T_0))$，在τ两侧分别用(13)并保一侧极限：

$$
-\int_0^{T_0}\zeta'(t)\mathcal E_A(t)\,dt
=-\int_0^{T_0}\zeta(t)\Phi_f(t)\,dt+\zeta(\tau)W_A .
\tag{18}
$$

对于$a<\tau<b$，等价积分身份是

$$
\mathcal E_A(b)-\mathcal E_A(a)
=-\int_a^b\Phi_f(t)\,dt+W_A .
\tag{19}
$$

不跨τ时删去工作跳项。有限持续控制可另由其实际控制功率积分替代；不能由瞬时权限免费推得它的执行过程。若换成选择性仪器，需连同结果权重及完整后态一起处理，复用624，不用单条条件态取代联合身份。

## 7. 已相关状态上的中插短窗因果

设脉冲后经过$s<d-2=1$。把这一时段相互作用表象的起点取为τ，只是计算方式，不物理重置ρ_-或场。交互项为

$$
V_{I,A}(t)=\lambda N_A\chi_t(h_A),\qquad
V_{I,B}(u)=\lambda N_B\chi_u(h_B),
\qquad 0\le t,u\le s.
$$

两源支撑的最短空间距离为1，而$|t-u|<1$，KG的Pauli–Jordan支持给

$$
[V_{I,A}(t),V_{I,B}(u)]=0 .
\tag{20}
$$

跨站材料也对易，故该短窗的相互作用传播子因果分解为A和B两个可交换因子。这一模型中每因子可由Weyl线性源定义，不以无界Dyson形式级数代替证明。

对任意B材料有界效果E_B，A因子与其自由演化效果对易，故完整Heisenberg效果可只含B因子及场，进而

$$
\left[
U_A,\ e^{isH}(E_B\otimes I)e^{-isH}
\right]=0,\qquad 0\le s<1 .
\tag{21}
$$

因此在**任意已经演化且可相关的正常联合态**上，τ施加或不施加(3)，B在此短窗的无条件统计完全相同。被动R可一并保留；A的迹保持材料仪器同理。局部记录若留A，不视为已广播给B。

(21)不要求场回到真空，不排除预先存在的相关性，也不把条件筛选的相关性误作信号。它只认证整站菜单；不认证A内部子球的微因果，不把延伸模式操作缩成点，也不继承1148的共同boost。

## 8. 交付范围与去重

本附件不重计线性源、一般Noether流或1147初始无信号工具。增量是同一原长任务候选中，实际不对易插入、全历史、包含交互存能的弱能源收支及实际脉冲功同时成立。主证明仍须给原λ/T0的全输出误差；本文没有独立签收其数值。

该有限分辨率结果不是连续局部Dirac父过程、完整海虚对、任意子区应力算符、能源仪器、自治壁或全装置四P。未完成这些较强接口不抹去(18)—(21)对当前实际菜单的结论。
