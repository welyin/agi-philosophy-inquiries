# 动态四极的正常场父过程与有限耦合记录余项

状态：未编号候选；不增加正式轮次、公理或科学组。本文把 [来源证明](source_proof.md) 的二阶系数接到一份明确的正常 Hilbert 空间演化，并给有限耦合余项。**固定世界管的数学父过程已定义；真实 COM、支撑、有限触发和最终仪器的匹配仍未证明。** 后者不是以下概率差的输入替代品，也不能因本数学构造成立就登记为已满足。

## 1. 共同对象与明确的父 Hamiltonian

采用两份七维头、设置／读取 qubit，以及连续、未作动量或占据截断的两个 TT 偏振的对称 Fock 空间。两头均有 $H_c=3I-J_x$，$P=|m_z=0\rangle\langle m_z=0|$；每个完整静止比较物的质量算符取

$$
M_a=I+H_{c,a},\qquad 1\le M_a\le7,
\qquad H_0=M_A+M_B+d\Gamma(\omega),\quad\omega(k)=|k|.
\tag{1}
$$

此处的正质量基线、固定身体标架及涂抹中心是有效材料输入，不是已经含有机械定位和反冲的微观构造。设置 qubit 自由能简并；其内部驱动能源没有被从 $M_a$ 中扣除。

固定 $q_*=1$、$e=\operatorname{diag}(1,-1,0)$、$r=1/4$、$R=4$。形状为归一的 $f(x)\propto(1-|x|^2/r^2)^4_+$。Fourier 变换采用 $\widehat f(k)=\int e^{-ikx}f(x)d^3x$、积分测度 $d\mu(k)=d^3k/(2\pi)^3$。在 Schrödinger 表示中令

$$
\begin{aligned}
\rho_a(k)&=e^{-ikX_a}\widehat f(k)
 [M_a-(e:kk)Z_aP_a],\\
v_{a,i}(k)&=-i e^{-ikX_a}\widehat f(k)e_{ij}k_jZ_a\dot P_a,\\
T_{a,ij}(k)&=e^{-ikX_a}\widehat f(k)e_{ij}Z_a\ddot P_a,\\
\dot P&=i[H_c,P],\qquad\ddot P=-[H_c,[H_c,P]].
\end{aligned}
\tag{2}
$$

这是自由表示下完整单极加四极的矩阵值来源；不把 $P$ 替成其期望。它的空间迹为零，且自由 Ward 身份成立。相互作用后的精确总应力与自由插入不是同一概念；本页只匹配弱耦合领先交换，不声称已完成所有阶 Ward 或非线性引力。

设 $\epsilon_s^{ij}(k)$ 是标准正 TT 偏振，$\sum_s\epsilon_s\epsilon_s^*$ 给 TT 投影。令 $f_{a,s}^{\mathrm{TT}}(k)=e^{-ikX_a}\widehat f(k)(e:\epsilon_s)/\sqrt{2\omega}$（单粒子空间测度为 $d\mu$），$\phi(f)=a(f)+a^\dagger(f)$。定义

$$
\begin{aligned}
V_1&=\sum_{a=A,B}Z_a\ddot P_a\otimes\phi(f_a^{\mathrm{TT}}),\\
W_\xi&=-\frac\xi2\int\frac{d\mu(k)}{k^2}
\operatorname{Sym}\left[\frac12\rho_{\mathrm{tot}}^\dagger\rho_{\mathrm{tot}}
-2v_{\mathrm{tot}}^\dagger\!\cdot v_{\mathrm{tot}}
+\frac12v_{\mathrm{tot},L}^\dagger v_{\mathrm{tot},L}\right],\\
H_\xi(g)&=H_0+gV_1+g^2W_\xi,\qquad g=\kappa/2,\quad 0\le\xi\le2.
\end{aligned}
\tag{3}
$$

所有同站自项、共同质量、质量与四极的交叉项都保留。$\operatorname{Sym}$ 固定 Hermitian 排序；分离站点的因子相互对易。式(3)采用 [940 的弱场传播／约束归一](../../../archive_935_955/research_note_940.md)，但不预置非线性 Einstein 作用。它是一个明确的比较族：$\xi\ne1$ 时允许非局部瞬时相互作用，正是待检验的候选性质，不先假设它符合因果组合公理。

若加入另外的局部同阶 seagull 或器件修正，应把它们放进同一 $W$ 并计其范数。只有已证明 A 设置差／B 交换子为零的局部项才不改变下面的二阶见证；本页数值常数只覆盖式(3)，不把未匹配修正的误差设成零。

## 2. 无场截断的自伴性、正性和显式常数

用 $H_c-3I=-J_x$ 的范数为3，得到 $\|\dot P\|\le6$、$\|\ddot P\|\le36$。设

$$
I_0=\|f\|_2^2=\frac{776160}{4199\pi},\qquad
I_2=\|\nabla f\|_2^2=\frac{10644480}{221\pi},
\qquad I_{-2}=\int\frac{|\widehat f(k)|^2}{k^2}d\mu(k).
\tag{4}
$$

由于 $|\widehat f|\le1$，在 $k_0=1/r=4$ 处分开积分，

$$
I_{-2}\le\frac{k_0}{2\pi^2}+\frac{I_0}{k_0^2},\qquad
\|f_a^{\mathrm{TT}}\|^2\le
2\left[\frac{k_0^2}{8\pi^2}+\frac{I_0}{2k_0}\right],\qquad
\|f_a^{\mathrm{TT}}/\sqrt\omega\|^2\le I_{-2}^{\rm upper}.
\tag{5}
$$

这些包括 $k\to0$ 和全部高动量尾，不是忽略 UV 或 IR 模。由 Parseval 及有限矩阵范数得到统一界

$$
L:=\sum_a\|Z_a\ddot P_a\|\|f_a^{\mathrm{TT}}\|\le300,
\qquad
S:=\sum_a\|Z_a\ddot P_a\|\|f_a^{\mathrm{TT}}/\sqrt\omega\|\le150,
$$

$$
\|W_\xi\|\le w:=2\,[49I_{-2}^{\rm upper}+194I_0+I_2]<60000.
\tag{6}
$$

例如 $\|\rho_{\rm tot}\|\le2(7+k^2)|\widehat f|$、$\|v_{\rm tot}\|\le12k|\widehat f|$ 给出式(6)；因此质量和全部自项确实进入常数。代码用 $\pi>3$ 的有理上界验证300、150、60000，无浮点常数承担证明。

线性场在自由场定义域上相对有界且相对界任意小，$W$ 有界，因此式(3)自伴。形式估计为

$$
|g\langle V_1\rangle|\le2|g|S\|d\Gamma(\omega)^{1/2}\psi\|\|\psi\|,
\qquad
H_\xi(g)\ge2+\tfrac12d\Gamma(\omega)-g^2(2S^2+w).
\tag{7}
$$

固定 $g=10^{-15}$ 时右侧严格大于 $I+\tfrac12d\Gamma(\omega)$。总生成元与其能源谱分布守恒；此正性是明确父模型的性质，不是完整实际物体已经匹配的证明。

## 3. 同一过程中的非选择 B 记录

准备映射 $E$ 将两头放在 $|-3\rangle$、场放在正常 Fock 真空，保留全部设置／读取 qubit 和任意参考 $R$。演化为同一 $e^{-iTH_\xi(g)}E$，$T=\pi$；输出仍含全部场、头和 qubit。B 效果为 $F_B=I_{\rm others}\otimes(I+Y_B)/2$，没有头末态、A 设置或场的后选。

见证取 A 的两已知 $Z=\pm1$、B 的 $+X$。两者具有相同头、质量、场和其他初始资料。自由单插入满足

$$
\langle P(t)\rangle=a(t)=\frac5{16}\sin^6t,\qquad
i\langle[Z_BP_B(t),F_B]\rangle=a(t).
\tag{8}
$$

此期望只在求二阶系数时使用。$P(t)$ 在整个七维头空间并不是紧支算符。真正的端点身份是 $\operatorname{ad}_{H_c}^{j}(P)|\pm3\rangle=0$（$j\le2$）；它保证当前准备下所用跨站单插入边界消失。没有将全算符强行截窗或删掉开关来源。对实际反馈后的末态也不声称头恰回到 $|+3\rangle$。

在 interaction picture 中直接展开实际概率，TT 交叉项为

$$
-g^2\int_{t>s}dt\,ds\,
\langle[V_A(s),[V_B(t),F_B]]\rangle.
\tag{9}
$$

用 $[Z_B,F_B]=-iX_B$ 及 $[\phi_A(s),\phi_B(t)]=iD^{\rm TT}_{BA,\rm ret}(t-s)$，该项的 A 设置差为 $-2g^2 I_{\rm TT,ret}[a,a]$。$W$ 的跨站四极项是 $-\xi Z_AZ_B(c_0P_AP_B+c_1\dot P_A\dot P_B)$，故 $ig^2\int[W(t),F_B]dt$ 的设置差为 $-2g^2\xi I_C[a,a]$。

其他项的领先相消有明确量词：本模型中 A 的自由质量来源不含 $Z_A$，在两个设置下相同；B 的质量来源与 $F_B$ 对易；同站 A 项与 $F_B$ 对易；同站 B 项对 A 设置无依赖。因此不能产生这一设置差，但仍在实际 H 和余项中保留。这个结论不自动覆盖尚未定义的真实触发器或记录器来源。

整个支撑窗类空分离：$T=\pi<7/2=R-2r$。连续 TT 的完整退迟核与所保约束核按940的身份共同使用，才有 $I_{\rm TT,ret}=-I_C$。所以从式(3)的同一演化得到

$$
\Delta p_\xi(g)
=-2g^2(\xi-1)I_C+\mathcal R_\xi(g),\qquad
I_C=-\frac{22700475}{27917287424}.
\tag{10}
$$

这不是把核振幅取模平方，也没有将期望源替换进全 Hamiltonian。独立关系 $\Delta p_1^{(2)}=0$ 仅是领先阶；本比较父没有宣称非微扰全阶 Einstein 因果性。

## 4. 真空准备、全部未知输入的严格余项

对场粒子数不超过 $n$ 的任意物质／参考向量，

$$
\|V_1\psi\|\le2L\sqrt{n+1}\|\psi\|.
$$

$W$ 不改变场粒子数。含 $l$ 次线性顶点、$k$ 次 $W$ 的 Dyson 项，对所有交错顺序求和后的范数至多

$$
\frac{u^l}{\sqrt{l!}}\frac{v^k}{k!},\qquad
u=2|g|LT,\quad v=g^2wT.
\tag{11}
$$

理由是每次线性场顶点至多增加一个粒子，给 $\sqrt{l!}$；时间单纯形与交错数给 $T^{l+k}/(l!k!)$。这使真空准备上的级数绝对收敛并与自伴演化相同；不要求 $V_1$ 在整个 Fock 空间有界。$R$ 只被张量恒等作用，界对全部声明物质输入／任意 $R$ 一致，不只对见证态成立。

截到 $g^2$（同时保 $V_1V_1$ 及一次 $W$），写输出为 $\psi_0+\psi_1+\psi_2+\epsilon$。对 $u,v<1$，可以保守取

$$
\|\psi_1\|\le u,\quad\|\psi_2\|\le d=u^2+v,\qquad
\|\epsilon\|\le r_3=
\frac{u^3}{1-u}+\frac{uv}{(1-u)(1-v)}+\frac{v^2}{1-v}.
\tag{12}
$$

三个尾项分别覆盖 $l\ge3,k=0$、$l\ge1,k\ge1$ 和 $l=0,k\ge2$。对任意效果，把二阶概率后的项逐项估计得单设置余项

$$
B_{\rm prob}=2r_3(1+u+d)+r_3^2+2ud+d^2,
\qquad |\mathcal R_\xi(g)|\le2B_{\rm prob}.
\tag{13}
$$

这里没有把全场末态投影到0或1粒子，也没有删噪声、控制关联或失配支路。场宇称还使本模型概率对 $g$ 为偶函数，但本页无需借此改善粗界。

## 5. 一份固定、极小且非器件性能宣称的数值证书

只取 $g=10^{-15}$、$L\le300$、$w\le60000$、$T=\pi<22/7$，其余物理参数不扫描。用有理数计算 $u,v,d,r_3,B_{\rm prob}$。当 $|\xi-1|\ge1/2$ 时，

$$
|\Delta p_\xi(g)|\ge g^2|I_C|-2B_{\rm prob}>0.
\tag{14}
$$

代码给精确有理下界并展示量级。差约在 $10^{-34}$，**不是已能实际使用或已经观测到的信号**，也不值得据此继续优化 g、头或仪器。它只说明该明确正常父模型的有限耦合余项不会自动吞掉二阶系数。

若后来独立证明两设置的完整实施各自与该父预测的概率误差不超过 $\eta_{\rm actual}$，式(14)再减 $2\eta_{\rm actual}$。可反推所需精度，不能把它当作已成立的材料合同，更不能将式(10)本身写成匹配公理。有限测量仍只能排候选区间，不能无限精确选择实数 $\xi$。

## 6. 当前真正完成与未完成

完成的是有限自治来源、连续正常传播场、全部约束自项及同一个非选择记录之间的一个**数学父过程**，并有任意内部输入／参考下的有限耦合余项。一般来源图范数逼近和 C¹ 编译已见623—625及1123，不重复计为新工具；这里只用它们的范围原则，实际误差直接由本父 Dyson 展开给出。

仍未完成：该固定世界管父与包含动态 COM、身体标架、真实支撑、新鲜局部设置、有限钟窗及真实末读的同一总过程之间的来源保持匹配。1171 的已知参数门族可以保存给定参数导数，但不能决定这里的物料背景耦合，亦不能预读未来的背景。当前不能由本页宣布实际无信号原则已经选择 $\xi$、认知公理已选引力核、或 GR／暗部门目标已完成。

[复算代码](response_check.py) · [保存结果](response_results.json) · [准入范围](admission_scope.md)
