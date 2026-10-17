# 1256：原第二记录流的局部canonical场、输运与完整来源

2026-10-16。依据[事前准入](admission.md)。新增认知公理0；CT仍为待验证猜想。认知动机是：原参与者换一种描述后，其相位、运动、记录及资源仍应属于同一过程。这里给可逆经典表示连接，不把表示所需的技术条件宣传为新物理定律。

## 1. 原材料和待补接口

保留[1245](../1245/proof.md)的第二连续记录流：固定正标签测度 $n_r(z)d^3z$，嵌入 $X_r^i(t,z)$、每标签机械动量 $P_i(t,z)$、内部相位 $\vartheta_L(t,z)$ 和动作 $I_L(t,z)$。完整辛一形式为

$$
\Theta_r=\int d^3z\,n_r(z)
 [P_i(z)\,dX_r^i(z)+I_L(z)\,d\vartheta_L(z)] .
\tag{1}
$$

下标L表示标签资料；不表示新的物种。原质量、法向机械能与约束是

$$
M_r=\frac1{10}+\frac{1+p(X_r)}{10}I_L,\qquad
E_{\mathrm{part}}=\sqrt{h^{ij}(X_r)P_iP_j+M_r^2},
$$
$$
H_r[N]=\int d^3z\,n_rN(X_r)E_{\mathrm{part}},\qquad
D_r[\beta]=-\int d^3z\,n_r\beta^i(X_r)P_i .
\tag{2}
$$

负号沿原ADM约定 $dx^i+\beta^idt$。原世界线内部相位率为 $d\vartheta_L/d\tau_r=(1+p)/10$，不是原第一钟的单位固有时变量；两流不能同时固定共动。

[937](../../archive_935_955/research_note_937.md)的根规约方法、[1255](../1255/proof.md)的原互动钟规约及一般cotangent变换均不重计。本轮补的是(1)—(2)这份**原记录载体**的完整局部场字典。旧文已警示不可漏内部梯度，但没有给下面可直接接入1255的全部公式。

## 2. 可逆域、逆标签和密度

取原第二流无焦散、定向保持、空间映射可逆的有限条带，定义

$$
Z^A(t,x)=(X_r^{-1})^A(t,x),\qquad
F^A{}_i=\partial_iZ^A,\qquad
E^i{}_A=(F^{-1})^i{}_A,
$$
$$
\eta(x)=n_r(Z(x))\det F(x)>0,\qquad
\vartheta(t,x)=\vartheta_L(t,Z(t,x)).
\tag{3}
$$

$P_i(x),I(x)$下文表示原标签资料与 $Z$ 的复合。原标签为周期坐标时，$Z$使用相容周期图及原绕数；没有引入全局实值环面坐标。局部表达在相容图间运输。

$\eta$是原标签数测度在坐标体积 $d^3x$ 下的密度。它不是proper数密度，也不是质量密度，更不是新的独立canonical变量。所需逆变换存在于 $\det F>0$、$n_r>0$ 的图；相位沿原角变量或有限lift处理，$I$与 $M_r$ 的原合法部门保留。图外的流体交叉/焦散需要其他描述，本轮不宣称覆盖它们。

## 3. 从辛一形式导出canonical动量

在固定标签 $z$ 取任意配置变分。由 $Z(X_r(z))=z$ 及 $\vartheta_L(z)=\vartheta(X_r(z))$，有

$$
\delta X_r^i\circ Z=-E^i{}_A\delta Z^A,
\qquad
\delta\vartheta_L\circ Z
 =\delta\vartheta-(\partial_i\vartheta)E^i{}_A\delta Z^A .
\tag{4}
$$

换积分变量 $n_r(z)d^3z=\eta(x)d^3x$ 后，(1)成为

$$
\Theta_r=\int d^3x\,[\Pi_A\,\delta Z^A
                         +\Pi_\vartheta\,\delta\vartheta],
$$
$$
\boxed{\Pi_\vartheta=\eta I,\qquad
\Pi_A=-\eta(P_i+I\partial_i\vartheta)E^i{}_A.}
\tag{5}
$$

这是一形式的精确恒等，无需删边界项。虽然 $\eta$ 随配置变化，它已是(5)动量定义的系数；取外微分时它的全部变化自然保留，不得另把它固定。

定义实际机械动量密度 $m_i$ 和原动量约束密度 $D_{r,i}$：

$$
\boxed{D_{r,i}=\Pi_A\partial_iZ^A
              +\Pi_\vartheta\partial_i\vartheta
              =-\eta P_i=-m_i.}
\tag{6}
$$

逆变换明确为

$$
I=\Pi_\vartheta/\eta,\qquad P_i=m_i/\eta,
\qquad X_r=Z^{-1},\qquad
\vartheta_L(z)=\vartheta(X_r(z)),
\tag{7}
$$

其余资料亦在 $X_r(z)$ 取值。把(6)—(7)代回(5)恢复同一 $\Pi_A$；自由度未增加。于是正规片内

$$
\{Z^A(x),\Pi_B(y)\}=\delta^A{}_B\delta^3(x-y),
\quad
\{\vartheta(x),\Pi_\vartheta(y)\}=\delta^3(x-y),
\tag{8}
$$

其余新基本括号为零。这不是把原每标签 $I$ 的带 $1/n_r$ 括号直接搬为密度括号，而是由完整辛形式得到。

内部相位梯度项是表示转换的一部分，不是新增机械力。若只在(6)保留该梯度，却仍以 $\Pi_A=-\eta P_iE^i{}_A$ 定义动量，就混用了两套变量；反之，正确的 $\Pi_A$ 固定后再漏掉(6)的相位项，也不再表示原机械动量。

## 4. 同一质量、局部能源与来源

令

$$
W:=\eta M_r=\frac{\eta+(1+p)\Pi_\vartheta}{10}>0 .
\tag{9}
$$

由(2)逐项换元，原能量和动量约束完全局部化为

$$
\boxed{\mathcal E_r=\eta E_{\mathrm{part}}
 =\sqrt{h^{ij}m_i m_j+W^2}
 =\sqrt{h^{ij}D_{r,i}D_{r,j}
      +\frac{[\eta+(1+p)\Pi_\vartheta]^2}{100}},}
$$
$$
H_r[N]=\int N\mathcal E_r\,d^3x,
\qquad D_r[\beta]=\int\beta^iD_{r,i}\,d^3x .
\tag{10}
$$

没有另选每格质量或删除静质量 $\eta/10$。$\eta$、$D_r$含 $Z,\vartheta$ 的一阶空间导数，故这是同一原流的局部有限jet表达；它仍是连续有效材料，不是有限点粒子或已经量子化的空间cell。

同一作用的Hilbert应力由(10)还原为

$$
\sqrt h\,T^r_{nn}=\mathcal E_r,\qquad
\sqrt h\,T^r_{ni}=-m_i=D_{r,i},\qquad
\sqrt h\,S^r_{ij}=\frac{m_i m_j}{\mathcal E_r} .
\tag{11}
$$

尤其 $T_{ni}$ 的符号不能由能源中的平方猜测。proper数密度 $\nu_r$ 和proper质量密度是

$$
\nu_r=\frac{\eta W}{\sqrt h\mathcal E_r},
\qquad \nu_r M_r=\frac{W^2}{\sqrt h\mathcal E_r}.
\tag{12}
$$

它们与原 $T_r^{\mu\nu}=\nu_rM_rU_r^\mu U_r^\nu$ 相同。

固定canonical资料时，字段源和逆度规变分满足

$$
\frac{\partial\mathcal E_r}{\partial p}
 =\frac{W\Pi_\vartheta}{10\mathcal E_r},
\qquad
\left.\delta\mathcal E_r\right|_{\delta h^{ij}}
 =\frac{m_i m_j}{2\mathcal E_r}\delta h^{ij} .
\tag{13}
$$

第一式乘 $N$ 等于原标签推前的 $N\eta(M_r/E_{\mathrm{part}})(I/10)$；场动量方程取其负号。原作用的 $\nabla_\mu T_r^{\mu\nu}=-\nu_r(I/10)\nabla^\nu p$ 与字段交换仍成立，不额外加 $\nabla I$ 力。

## 5. 密度、载体和内部相位的共同输运

对任意 $\delta Z$，Jacobian恒等式给

$$
\delta\eta=\partial_i(\eta E^i{}_A\delta Z^A).
\tag{14}
$$

它同时包含 $n_r(Z)$ 非常数时的导数。可用 $\partial_i[(\det F)E^i{}_A]=0$ 的cofactor恒等式直接展开核对；不能只在常密度上忽略该项后宣称一般变换成立。

先对原未规约 $H_r[N]+D_r[\beta]$ 求第二流体的Hamilton方程，将 $N,\beta$ 视作原lapse/shift乘子。实际坐标速度为

$$
v_r^i=N\frac{h^{ij}m_j}{\mathcal E_r}-\beta^i
     =-N\frac{h^{ij}D_{r,j}}{\mathcal E_r}-\beta^i .
\tag{15}
$$

(8)—(10)给

$$
\dot Z^A+v_r^i\partial_iZ^A=0,\qquad
\dot\eta+\partial_i(\eta v_r^i)=0,
$$
$$
\dot\Pi_\vartheta+\partial_i(\Pi_\vartheta v_r^i)=0,
\qquad
(\partial_t+v_r^i\partial_i)I=0 .
\tag{16}
$$

后两式利用质量不依赖 $\vartheta$。相位率则为

$$
\boxed{(\partial_t+v_r^i\partial_i)\vartheta
 =N\frac{W}{\mathcal E_r}\frac{1+p}{10}
 =N\frac{M_r}{E_{\mathrm{part}}}\frac{1+p}{10}.}
\tag{17}
$$

这正是原第二记录器固有时率乘 $d\tau_r/dt=NM_r/E_{\mathrm{part}}$。它通常既不等于1，也不等于在第一钟材料图中未经运动修正的静止率。

可逆变换同时保原机械方程；写为Eulerian资料是

$$
\begin{aligned}
(\partial_t+v_r^j\partial_j)P_i={}&-E_{\mathrm{part}}\partial_iN
 -\frac{N}{2E_{\mathrm{part}}}P_jP_k\partial_i h^{jk}
 -N\frac{M_r}{E_{\mathrm{part}}}\frac I{10}\partial_i p
 +P_j\partial_i\beta^j .
\end{aligned}
\tag{18}
$$

右端是在原位置变分固定标签内部 $I$ 的公式。由(5)的梯度修正、(14)及(16)，内部准备的Eulerian梯度不会成为额外力；(18)也可直接从(2)经(7)运输得到。若改 $n_r(Z)$ 或相位分布，必须使用同一变换，不能删除这些相消所需的项。

## 6. 接入原1255的材料时间Hamiltonian

1255已给原第一互动钟的材料图及

$$
\mathcal R=\sqrt{C_0^2-h^{ij}D_{0i}D_{0j}}>0,
\quad H_{\rm mat}=\int[n_cU-\mathcal R]d^3x,
$$
$$
N=-C_0/\mathcal R,\qquad
\beta^i=h^{ij}D_{0j}/\mathcal R .
\tag{19}
$$

现在同一第二流体项精确写成

$$
C_0=C_g+C_{p,X,\phi}+\mathcal E_r,
\qquad D_{0i}=D_{g,i}+D_{p,X,\phi,i}+D_{r,i},
\tag{20}
$$

没有新增参考流，也没有把第二流换成第一钟的共动副本。其他字段、初始约束、正根及原 $J=\mathcal R/n_c-U>0$ 的范围沿1255保留。

对于第二流体变量的变化，原根的完整微分包含

$$
\left.\delta H_{\rm mat}\right|_r
 =\int[N\,\delta\mathcal E_r+\beta^i\delta D_{r,i}]d^3x .
\tag{21}
$$

因此(15)—(18)可代入(19)的同一 $N,\beta$。**(21)是原根微分的相应部分，不是另把 $H_{\rm mat}$ 换为 $N\mathcal E_r+\beta D_r$ 后再变分。** 后一种做法会产生不相符的乘子导数，并漏掉第一钟的 $n_c\delta U$ 及根的显式逆度规项。字段/度规变化时必须使用1255完整变分式；(13)只是其中原第二流体的源贡献。

此时母Hamiltonian全部使用材料坐标上的局部字段及有限阶空间导数：引力曲率、原字段梯度、$Z,\vartheta$及其canonical密度动量，还有原第一钟的标签内部变量。两个流仍能相对运动；原1255初片上第二流 $P_r=0$ 时，(15)仍给 $v_r=-\beta\ne0$。

这个结果解决的是**同一连续母作用的局部表示**，没有证明有限格点上的约束闭合、量子CCR表示、正源域尾或实际cell仪器。它提供后继有限任务承接的明确起点，不能替代那些证明。

## 7. 固定删除见证：同canonical资料下漏相位项

使用原第二流密度 $n_r=1/10$。在一处初片取

$$
Z=x,\quad \eta=\frac1{10},\quad I=1,\quad
p=0,\quad h^{ij}=\delta^{ij},\quad P_i=0,
\qquad \partial_1\vartheta=1,\quad\partial_{2,3}\vartheta=0 .
\tag{22}
$$

这可由周期为 $L=10$ 的光滑准备

$$
\vartheta(x)=\frac L{2\pi}\sin\frac{2\pi x^1}{L}
$$

在 $x^1=0$ 实现，无须不合法的全局非周期实相位。质量不含该相位，故此准备不增加原机械动量、质量或Hilbert应力；相对于原1255初片只是更改内部相位资料，不能称完整初态没变。

正确canonical资料为

$$
\Pi_\vartheta=\frac1{10},\qquad
\Pi_A=(-\tfrac1{10},0,0),\qquad
D_{r,i}=0,\qquad m_i=0,\qquad
W=\mathcal E_r=\frac1{50} .
\tag{23}
$$

现在**保持(23)的canonical资料不变**，只在机械动量恢复式中漏掉相位项：

$$
m_i^{\rm bad}=-\Pi_A\partial_iZ^A,
\quad D_{r,i}^{\rm bad}=-m_i^{\rm bad}.
$$

于(22)的点得到

$$
m^{\rm bad}=(\tfrac1{10},0,0),\qquad
D_r^{\rm bad}=(-\tfrac1{10},0,0),\qquad
\mathcal E_r^{\rm bad}=\frac{\sqrt{26}}{50},
$$
$$
\boxed{\mathcal E_r^{\rm bad}-\mathcal E_r
 =\frac{\sqrt{26}-1}{50}>\frac2{25}.}
\tag{24}
$$

还产生原本不存在的 $\sqrt hS_{11}^{\rm bad}=1/(2\sqrt{26})$ 及非零 $T_{n1}$。因此错误只在动量定义中出现，也会污染原几何和字段来源。

这里不另设 $\Pi_A^{\rm bad}=0$：那是同时改变canonical变换/初值的另一种错误，不能与本固定资料删除混算。(24)仅拒绝漏项字典，不反对相位梯度准备、Eulerian变量或另一份自洽材料模型。

相位不进入原 $M_r$，所以正确(22)—(23)可与1255原初始约束和来源相容。本轮没有据此领取新有限未来几何预测、原末读精度或量子仪器资格。

## 8. 证据分层与下一步边界

作者代码计划分三组：canonical与逆变换；输运/源；固定资料删除。相应的解析量词由(3)—(21)承担，有限精确核验只检查身份与实现。实际代码、结果和独审由主任务另行保存；其是否通过以正式验收文件为准，不以本计划代替运行结果。

得到的是：同一原第二记录流的机械运动、内部相位、全部源及与1255的接续，可以完整写成一份局部canonical场表示。这为局部母Hamiltonian提供明确变量，但不增加物理原则、物种或自由能源，也不是从认知原则新推出一条物理定律。

下一项应先判断一份小固定空间任务的离散作用、局部效果、约束误差和完整源能否共同取得非空证书。量子准备、CCR/有限类型、实际材料时钟及Jacobson热—几何角色未由本变换获得。若候选仍无明确母对象，就停在具体缺项，不把连续流的所有技术问题变成整个研究纲领的前置要求。
