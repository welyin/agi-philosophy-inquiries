# 1282 证明：有限正常态的共同幅度校准选择同一规范边比较

2026-10-17。状态：待独立审阅。新公理0。本稿只给公开候选类别内的反向分类；[报告](../research_note_1282.md)、[代码](check.py)、[结果](results.json)、[下一项](NEXT.md)。

## 1. 旧对象、独立输入与新问题

继承[557 §2](../../archive_554_584/research_note_557.md)的两个实四维双重态 $X_1,X_2$、singlet $s_1,s_2$、动态链路 $U\in SU(2)$，以及同一两端局部描述作用

$$
(X_1,X_2,U)\longmapsto
(R(g_1)X_1,R(g_2)X_2,g_1Ug_2^{-1}).
\tag{1}
$$

$R$ 是左四元数乘法对应的实四维正交表示。源空间仍为 $L^2(\mathbb R^{10}\times SU(2))$，Haar测度归一；物理态支持在局部Gauss不变部门。原正动能、正四次势 $W$ 和群动能不改变。仅开放原来直接指定的二次边势：

$$
F=\frac12X_1^T A(U)X_1-X_1^T B(U)X_2
  +\frac12X_2^T D(U)X_2\ge0.
\tag{2}
$$

$A,D$ 为实对称矩阵，$B$ 为一般实矩阵，系数光滑依赖 $U$；本候选 $F$ 不依赖旧singlet $s_i$ 或其它变量。若允许系数依赖singlet，可能留下新的函数自由，本轮不排除它。无导数、齐次二次阶、逐点非负及式(1)下不变均为本轮公开建模条件，不是已经从现行认知公理推出的全部物理规律。

认知动机是：差别分辨率相同，增大双方共同内容的幅度不应制造额外失配成本。可检验版本是下面两份**有限正常态**的 $F$ 均值相同；不要求实际量子态支撑在 $X_1=R(U)X_2$ 的零测集上。

## 2. 完整局部协变保留四元数交叉自由

由式(1)逐项比较系数，独立右作用使 $A(U)$ 为常矩阵；左作用使它与全部 $R(g)$ 对易。交换两端同理处理 $D$。实四维双重态的对易代数是右四元数代数

$$
\operatorname{Comm}(R)=\operatorname{span}_{\mathbb R}
\{I,J_1,J_2,J_3\},\qquad
J_j^T=-J_j,\quad
J_jJ_k+J_kJ_j=-2\delta_{jk}I.
\tag{3}
$$

右乘的乘法次序与左乘相反，不在证明中混用其结构常数。对易代数的实对称部分只有 $\mathbb RI$，故 $A=aI,D=dI$。交叉块不必对称：

$$
B(U)=B_0R(U),\qquad
B_0=b_0I+\sum_{j=1}^3b_jJ_j,\qquad
B_0B_0^T=(b_0^2+|\mathbf b|^2)I.
\tag{4}
$$

例如取 $g_1=I,g_2=U$ 将 $U$ 送到单位元，给 $B(U)=B(I)R(U)$；再取 $g_1=g_2=g$，得到 $B(I)$ 与 $R(g)$ 对易。实不可约并不意味着一般实对易矩阵只能是标量；若删去三个 $J_j$，会在校准前偷删合法候选。

记 $\beta_B^2=b_0^2+|\mathbf b|^2$。二次型矩阵

$$
M(U)=\begin{pmatrix}aI&-B_0R(U)\\-R(U)^TB_0^T&dI\end{pmatrix}
\tag{5}
$$

的两个本征值分别为

$$
\lambda_\pm=\frac{a+d\pm\sqrt{(a-d)^2+4\beta_B^2}}2,
\qquad
M(U)\ge0\ \Longleftrightarrow\ a,d\ge0,\ ad\ge\beta_B^2.
\tag{6}
$$

每个本征值按重数4计算，退化处允许合并。该判据包含 $a=0$ 或 $d=0$，不借不可逆Schur补排除边界。

## 3. 同一Gauss空间中的有限正常校准族

定义相对于原链路的共同与差别变量

$$
y=\frac{X_1+R(U)X_2}{\sqrt2},\qquad
z=\frac{X_1-R(U)X_2}{\sqrt2}.
\tag{7}
$$

令 $s=\sigma^2>0,t=\delta^2>0$，都为有限数。取

$$
\psi_{s,t}(X_1,X_2,U)=
\frac{1}{(2\pi s)(2\pi t)}
\exp\!\left(-\frac{|y|^2}{4s}-\frac{|z|^2}{4t}\right).
\tag{8}
$$

另乘两个singlet的归一中心Gaussian，每个概率方差1。对每个 $U$，式(7)是正交变量变换；概率分布下 $y,z$ 独立且各有协方差 $sI_4,tI_4$。积分 $U$ 的归一Haar测度后总范数仍为1。

式(1)使 $y,z$ 同作 $R(g_1)$ 变换，故该波函数严格Gauss不变，不只是密度与群作用对易。每个有限正 $s,t$ 下，它关于实变量一致Schwartz、关于紧群光滑；任意固定多项式、动量或群微分的有限次矩都有限。因此它属于原同一物理 $L^2$ 空间，原固定正多项式Hamiltonian的任意有限整数阶矩也有限。

这是一份公开的正常测试态族，**未证明原自主Hamiltonian自行制备它，也未授予免费准备或任意读口**。波函数通常有物质—链路关联，不可换成两端与链路独立的裸乘积。没有采取 $s\to\infty$、$t\to0$ 或delta准备极限。

## 4. 两个有限均值选择原边势

令

$$
q=a+d-2b_0\ge0,\qquad r=a+d+2b_0\ge0,
\qquad V=\frac{a-d}{2}I+\sum_jb_jJ_j.
\tag{9}
$$

在 $(y,z)$ 坐标中，$M$ 变为

$$
\widetilde M=\begin{pmatrix}qI/2&V\\V^T&rI/2\end{pmatrix},
\qquad
F=\frac q4|y|^2+y^TVz+\frac r4|z|^2,
\qquad
\langle F\rangle_{s,t}=sq+tr.
\tag{10}
$$

交叉项平均为零由这份独立中心Gaussian给出，并非一般输入上的删关联许可。固定同一个 $t>0$，取任意两个不同有限 $s_2>s_1>0$，则

$$
\langle F\rangle_{s_2,t}-\langle F\rangle_{s_1,t}
=(s_2-s_1)q.
\tag{11}
$$

另一方面，正性式(6)给

$$
(a-d)^2+4|\mathbf b|^2
\le(a+d)^2-4b_0^2=qr.
\tag{12}
$$

故两均值相同当且仅当 $q=0$；进而 $a=d=b_0=\kappa\ge0$、$\mathbf b=0$。反向代回也成立。因此

$$
\boxed{\quad
\langle F\rangle_{s_2,t}=\langle F\rangle_{s_1,t}
\ \Longleftrightarrow\
F=\frac\kappa2|X_1-R(U)X_2|^2,\quad\kappa\ge0.
\quad}
\tag{13}
$$

被选形式正是557原边势，系数对应为 $\kappa_{\rm here}=\kappa_{557}v_{557}$；没有消去旧体积参数 $v$。两份态上的成本都等于 $4\kappa t$，在 $\kappa,t>0$ 时严格正。校准的是**共同幅度改变不增加差别成本**，不是正常态具有精确零成本。

## 5. 非零有限读数误差与系数界

令实际均值读数为 $\widehat f_i$，在声明置信事件内满足 $|\widehat f_i-\langle F\rangle_{s_i,t}|\le\eta_i$。误差须包含被授予合同中的准备偏差、校准及采样，不由浮点残差代替。设 $\Delta s=s_2-s_1>0$，则

$$
0\le q\le\varepsilon:=
\frac{|\widehat f_2-\widehat f_1|+\eta_1+\eta_2}{\Delta s}.
\tag{14}
$$

若另有独立可核上界 $r\le R_{\max}$，式(12)给

$$
|a-d|\le\sqrt{\varepsilon R_{\max}},\qquad
|\mathbf b|\le\frac12\sqrt{\varepsilon R_{\max}}.
\tag{15}
$$

例如，由非负性和任一份有限均值的上界即可取 $r\le(\widehat f_1+\eta_1)/t$；亦可使用独立来源或系数合同，不能凭希望宣称 $r$ 有界。

取相比较的目标系数 $\kappa_*=r/4$，$M_*$ 为式(13)对应矩阵。由于 $\|V\|^2=((a-d)^2+4|\mathbf b|^2)/4$，

$$
\|M(U)-M_*(U)\|
\le\frac{q+\sqrt{qr}}2
\le\frac{\varepsilon+\sqrt{\varepsilon R_{\max}}}2.
\tag{16}
$$

此界对全部 $U$ 一致。它是有限系数矩阵范数界，**不是无界势算符的全空间范数界**。对任意物理态及任意旁挂伙伴，若 $\langle|X_1|^2+|X_2|^2\rangle\le B_2$，包括物质—链路相关输入，则点态二次型界给

$$
|\langle F-F_*\rangle|
\le\frac{B_2}{4}(q+\sqrt{qr})
\le\frac{B_2}{4}
(\varepsilon+\sqrt{\varepsilon R_{\max}}).
\tag{17}
$$

两个待比模型使用同一物理输入来评价式(17)。该式不证明有限时间通道距离、参考制备稳定性或原系统演化始终满足同一 $B_2$；需要这些结论时另列任务，不能直接从均值误差推出。

## 6. 两份准备并非同总能或免费

以 $a_{\rm kin},b_{\rm el}>0$ 表示557坐标动能及群动能系数，避免与式(4)系数混淆。式(8)的两项真实动能为

$$
\left\langle-a_{\rm kin}(\Delta_{X_1}+\Delta_{X_2})\right\rangle
=a_{\rm kin}\left(\frac1s+\frac1t\right),\qquad
\langle b_{\rm el}(-\Delta_G)\rangle
=\frac{3b_{\rm el}(s-t)^2}{16st}.
\tag{18}
$$

后一式沿用557的 $-\sum_A J_A^2=3I/4$ 归一。令 $c=X_1\cdot R(U)X_2$，则

$$
\log\psi= -\frac{|X_1|^2+|X_2|^2}{8}(s^{-1}+t^{-1})
+\frac c4(t^{-1}-s^{-1})+\mathrm{const},\qquad
\left\langle\sum_A|J_Ac|^2\right\rangle=3st.
\tag{19}
$$

因为 $\sum_A|J_Ac|^2=(|X_1|^2|X_2|^2-c^2)/4$，括号内等于 $|y|^2|z|^2-(y\cdot z)^2$，均值 $16st-4st=12st$。实波函数分部积分给式(18)，不存在被删去的自由群动能。

完整账还包括两个singlet的动能 $a_{\rm kin}/2$、原 $W$ 及式(10)的边势。记 $h=s+t$；每端

$$
\langle|X_i|^2\rangle=2h,\qquad
\langle|X_i|^4\rangle=6h^2,\qquad
\langle s_i^2\rangle=1,\quad\langle s_i^4\rangle=3.
\tag{20}
$$

按557的 $W=(v/4)\sum_i (|X_i|^2-u_h,s_i^2-u_s)^T L(|X_i|^2-u_h,s_i^2-u_s)$，同一全势均值为

$$
\langle W\rangle=\frac v2\left[
L_{11}(6h^2-4u_hh+u_h^2)
+2L_{12}(2h-u_h)(1-u_s)
+L_{22}(3-2u_s+u_s^2)\right].
\tag{21}
$$

式(18)—(21)加总为这份测试态的完整原源平均能源。它们一般随 $s$ 改变，不能从式(11)为零宣布两准备总能源相同。这里核对的是态与原源的身份，不是提供完整实验准备工厂。

## 7. 删除条件的准确反例

1. **删正性：** $F=(|X_1|^2-|X_2|^2)/2$ 保式(1)和二次阶，两测试均值皆0，但不属于式(13)。正负项抵消。
2. **删完整局部描述协变：** 令 $A_0=\operatorname{diag}(1,2,3,4)$，$F=\tfrac12(X_1-R(U)X_2)^TA_0(X_1-R(U)X_2)$。它正且两测试均值相同，却有方向依赖。独立右端描述仍相容，缺失的是左端全部旋转相容，不能含混说所有规范性都消失。
3. **删二次有效阶：** $F=\kappa|z|^2+4\lambda|z|^4$，$\kappa,\lambda>0$，保正性、完整协变及校准；均值 $4\kappa t+96\lambda t^2$。因此本轮没有选出全部非线性作用。
4. **不要求真实非零响应：** $F=0$ 合格，$\kappa$ 的非零性及大小未获选择。声明独立非零比较响应才可排除该分支。
5. **删相同差别宽度：** $a=d=1,B_0=0$ 时 $q=r=2$。$(s,t)=(1,2),(2,1)$ 的均值均为6，却有 $q=2$；两次必须固定 $t$，或把其变化及误差明确校正。

反例只删除各自标明的条件，不否定认知纲领。尤其第三项不是有限有效二次近似失败，而是说明其范围不能无条件外推。

## 8. 正式贡献与可复用旧接口

[544](../../archive_531_553/research_note_544.md)是标量运行；[548 §2](../../archive_531_553/research_note_548.md)的轨道传递工具直接复用；[557 §2](../../archive_554_584/research_note_557.md)及[559 §2](../../archive_554_584/research_note_559.md)原先直接输入边势。933的运输分类、959的Gauss父表示、1281的有限块成本和1016的全部局部作用变分均不重算。

本轮新连接是：**在557原载体与完整局部描述作用中，两份正常态的有限共同幅度校准，排除四元数交叉自由并选择原边势，同时给有限精度下的系数与有限二阶矩预测界。** 选中 $\kappa>0$ 后，557的Gauss读口、源二阶矩与注能结论，559的有限时探针结论，按原权限、原误差及原参数重标范围直接继承。它们不是本轮的新科学组。

未知候选 $F$ 的测量权限需作为有效合同明示：557已有选定 $F$ 的数学仪器不能自动成为辨识任意未知作用的免费黑箱。全部均值证书的实际误差已留在式(14)中。没有新增认知公理、改变原准备记录、生成SU(2)、空间维数、时间分量、规范场独立动力学或Yang–Mills场论。
