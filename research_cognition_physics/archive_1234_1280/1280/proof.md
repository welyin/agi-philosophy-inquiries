# 1280证明：自然共形参考、绝对来源与有限模—几何余额

2026-10-17。规格见[fixed_spec.md](fixed_spec.md)，准入见[admission.md](admission.md)。本轮复用1279的经典解，改用明确bulk参考；不重新推导Einstein方程，不把成熟共形场论结果计为原创。

## 1. 同一几何与不同于1279的参考合同

固定签名 $(-,+,+,+)$、$4\pi G=1/36000$、$\Lambda=0$。取最小无质量实标量作用及

$$
ds^2=-2du\,dv+r(u)^2(dx^2+dy^2),\qquad
\varphi_{\rm cl}=f(u),\qquad
r''=-q r,\quad q=f'^2/36000,\quad r(0)=1,\ r'(0)=0.
\tag{1}
$$

$f$是1279已固定的脉冲：令 $x=u+2$，在 $0<x<1$ 上 $f=\exp[1-1/(4x(1-x))]$，其余为零。继承其解析证书

$$
4\le E:=\int f'^2du<36,\quad
0<A_q:=E/36000\le1/1000,\quad
499/500\le r\le1\quad(-2\le u\le0),\quad 0\le r'\le A_q.
\tag{2}
$$

这里明确取横向 $\mathbb R^2$，不沿用周期环面。固定有限横向测试域不会加入墙或周期边界；一旦紧化，参考须另核，不能领取零Casimir来源。

在 $u\le-2$，$r(u)=r_-+b(u+2)$，其中

$$
b=r'(-2)=\int_{-2}^{-1}q r\,du>0,\qquad
u_c=-2-r_-/b<-2.
\tag{3}
$$

故最大正Rosen支为 $(u_c,\infty)$。定义

$$
\eta(u)=\int_0^u\frac{dt}{r(t)^2},\qquad
\eta(0)=0,\qquad
g=r^2[-2d\eta\,dv+dx^2+dy^2].
\tag{4}
$$

过去 $r=b(u-u_c)$ 给 $\eta\to-\infty$；未来 $u\ge-1$ 时 $r=1$、$\eta=u\to\infty$。所以(4)将该支共形映到整个四维Minkowski空间；它不是跨越 $r=0$ 的坐标延拓。该支因而全局双曲，所用正共形因子在每个紧域上光滑且有正下界。

## 2. 正Hadamard参考和处方

令 $X=(\eta(u),v,x,y)$，$W_M$为四维Minkowski真空二点分布，定义

$$
W_c(p,p')=\frac{W_M(X,X')}{r(u)r(u')}.
\tag{5}
$$

对测试函数相应作体积及共形权运输，(5)保持正性、CCR和KG方程，定义准自由参考 $\omega_c$。共形映射保持Hadamard微局部谱条件；此处使用 [Pinamonti 0806.0803v2](https://arxiv.org/html/0806.0803) §2 Lemma 2.4。由于本背景 $R=0$，质量零最小与共形KG算符相同，因此它也是原最小材料的合法Hadamard态。这一步不声称两种Hilbert应力相同。

处方在计算前固定：平直Minkowski真空的应力和Wick平方为零；Wick平方取共形协变选择；共形应力取标准局部异常变换，最小应力以同一兼容Wick处方定义。Pinamonti §3.3式(9)给Wick平方的局部自由度为 $\alpha R$，共形协变选择存在；两端标量曲率都为零，故

$$
\omega_c(\varphi^2)_{\rm ren}=0.
\tag{6}
$$

不是针对工作态再拟合减除常数。重整化 $G,\Lambda$ 按(1)保持；未登记的 $cG_{ab}$ 或额外维度参数不属于本轮可任意添加的自由度。

## 3. 绝对参考应力及最小工作来源

共形真空应力的成熟局部变换使用 [Brown–Cassidy, Phys. Rev. D16,1712](https://doi.org/10.1103/PhysRevD.16.1712)；可读公式复核见 [Awad–Johnson hep-th/0008211v3](https://arxiv.org/pdf/hep-th/0008211) §3.1式(17)、§3.2式(24)。后者的文献37是Brown–Cassidy D15,2810，不能误写为D16,1712。此处只用：变换结果由平直参考应力及 $H^{(1)},H^{(3)}$ 局部张量组成。

具体地，令 $\ell_a=(du)_a$。它协变常数，曲率满足

$$
R_{ab}=2q\ell_a\ell_b,\qquad R=0,\qquad C_{abcd}=0,
\qquad R_a{}^cR_{cb}=R_{ab}R^{ab}=0,
\qquad \Box R_{ab}=0.
\tag{7}
$$

最后一式不是“标量不变量零”的猜测：$\nabla_cR_{ab}=2q'\ell_c\ell_a\ell_b$，再求导后收缩 $\ell_c\ell^c=0$。因此由 $R,\nabla\nabla R,\Box R,R R_{ab}$ 组成的 $H^{(1)}$ 全零；由 $R_a{}^cR_{cb},RR_{ab},g_{ab}R_{cd}R^{cd},g_{ab}R^2$ 组成的 $H^{(3)}$ 全零。标准四维曲率平方重整化歧义也消失：$R^2$变分为零，Weyl平方的Bach张量因 $C=0$ 为零，Euler项无bulk变分，Ricci平方可由这三者组合。这里没有声称任意新增高阶作用或非局部修改均受本结论约束。

平直参考应力为零，于是 $\omega_c(T^{\rm conf}_{ab})_{\rm ren}=0$。在兼容处方中，原最小与辅助共形Hilbert张量之差为改进项（可能的标准局部几何歧义在(7)上也为零）：

$$
T^{\rm conf}_{ab}-T^{\rm min}_{ab}
=\frac16\bigl(G_{ab}+g_{ab}\Box-\nabla_a\nabla_b\bigr)\varphi^2.
\tag{8}
$$

由(6)得

$$
\omega_c(T^{\rm min}_{ab})_{\rm ren}=0.
\tag{9}
$$

这是选定参考和已固定兼容处方下的绝对来源结论，不是仅有“无粒子产生”，也不要求对任意参考成立。

以光滑经典解 $f$ 作代数相干位移 $\varphi\mapsto\varphi+f$，得 $\omega_f$。它保持Hadamard性；在任一相对紧区域，可将对应Cauchy数据作紧支延拓而获得相同局部位移，故局部正常。整个平面位移不声称有限全局Fock范数或有限总能。因参考一次均值零，原**最小**工作来源为

$$
\omega_f(T^{\rm min}_{ab})_{\rm ren}
=\omega_c(T^{\rm min}_{ab})_{\rm ren}+T^{\rm min}_{ab}[f]
=f'^2\ell_a\ell_b.
\tag{10}
$$

工作态的改进项一般不为零，不能把 $T^{\rm conf}[f]$ 代替(10)。式(1)、(9)—(10)给本限定材料的一份自洽**平均源**半经典几何；没有另加或漏掉本模型内的参考平均应力。它不证明涨落小、稳定性、量子引力动力学或加入仪器后的解。

## 4. 同一bulk参考的类光限制及有限横向解释

在 $v=0$，规范场 $\psi=r\varphi$ 经(4)就是Minkowski标量的类光数据。采用去除零频问题的lightfront导数/CCR代数，其参考由(5)诱导，而非再次任意选频率。复用 [Halpern–Neiman 1502.04106v2](https://arxiv.org/html/1502.04106) §II、IV的半面真空：负半面 $\eta<0$ 的未来向模参数为 $s_m=-\log(-\eta)$。

下文首先定义**每单位横向坐标面积的响应密度**。在有限可测横向集合 $A_\perp$ 上，也可限定到它的lightfront子代数：参考二点核的横向 $\delta^2(x-y)$ 给正频单粒子空间按 $L^2(A_\perp)\oplus L^2(A_\perp^c)$ 分解，Fock真空及半面模作用相应因子化。位移在该集合内的横向因子 $1_{A_\perp}$ 属 $L^2$，纵向 $rf$ 光滑紧支，其一粒子位移范数有限。于是该子代数内相干响应正常且等于 $|A_\perp|$ 乘下列密度。此陈述不是在bulk制造一条不光滑边缘脉冲，也没有加入真实墙或授予新的读口。

同一参考的半面模能响应是

$$
\Delta K_{\eta}/{|A_\perp|}=2\pi M_\eta,\qquad
M_\eta=\int_{-2}^{-1}\xi_\eta(rf)'^2du,
\qquad\xi_\eta=-\frac{\eta}{\eta'}=-\eta r^2>0.
\tag{11}
$$

若只使用密度，式左记作 $\Delta k_\eta$。积分区间缩到脉冲是因位移紧支；它不是任意有限纵向小区间的模Hamiltonian，也不是整个无限平面的有限总模能。

## 5. 新的有限模权界与物理能流

在脉冲上 $r$向未来不减。由 $\eta(u)=-\int_u^0r(t)^{-2}dt$ 及 $r(u)\le r(t)\le1$，

$$
r(u)^2\le\frac{\xi_\eta(u)}{-u}\le1,\qquad
0\le1-\frac{\xi_\eta}{-u}\le
\delta_\eta:=1-(499/500)^2=\frac{999}{250000}.
\tag{12}
$$

物理未来向boost取 $\chi=-\kappa u\partial_u$，由绝对工作来源(10)独立定义

$$
Q_\chi/(\kappa|A_\perp|)=B:=\int_{-2}^{-1}(-u)r^2f'^2du>0.
\tag{13}
$$

令 $J=\int(-u)(rf)'^2du$。沿用1279已证明的紧支分部积分，而不沿用其参考：

$$
C:=J-B=\int[rr'+(-u)q r^2]f^2du\ge0,\qquad
\frac CB\le\delta_C=\frac{375}{498002}.
\tag{14}
$$

式(12)给 $0\le J-M_\eta\le\delta_\eta J$。故本轮连续严格界为

$$
\frac{|M_\eta-B|}{B}
\le\delta_C+\delta_\eta(1+\delta_C)
=\frac{591628623}{124500500000}
<\frac1{200}.
\tag{15}
$$

严格余量为 $30873877/124500500000$。不同参考使匹配预算变宽，但本轮新增了明确bulk态及绝对参考源证书；没有把参考变化隐去以声称原误差不变。

## 6. 独立几何余额与固有钟

取 $a=r^2$、$\theta=2r'/r$、剪切零、$\theta(0)=0$，在 $[-2,0]$ 上独立计算

$$
E_g=\int_{-2}^{0}(2u+2)aR_{uu}du
-\int_{-2}^{0}(u+2)a\theta^2/2\,du,
\qquad \Delta a-E_g=8\pi G B.
\tag{16}
$$

这是旧有限Raychaudhuri恒等式，不是新拟合。已采用 $s=1/(4G)$，由(15)得到

$$
|\Delta k_\eta-s(\Delta a-E_g)|<\frac1{200}(2\pi B).
\tag{17}
$$

面积由(1)确定，模量由(5)、(11)确定，来源由(10)确定；三者并未相互定义。此式验证已选Einstein有效解上的局部共同接口，不是第二次引力生成。

纵向 $x=y=0$ 沿1279曲线 $u=-e^{-\kappa\tau}/(\sqrt2\kappa)$、$v=e^{\kappa\tau}/(\sqrt2\kappa)$，$\tau$仍为精确固有时、几何加速度为 $\kappa$。比较 $\tau_m=-\log(-\eta)/\kappa$，

$$
1\le\frac{d\tau_m}{d\tau}=\frac{-u}{\xi_\eta}
\le\frac{250000}{249001}.
\tag{18}
$$

这仅是脉冲对应窗的参数率合同；不证明沿该时样世界线的二点函数精确proper-time KMS、不赋予真实探针响应或内部钟来源。HG的温标角色已采用，仍须区分本具体参考已证的null模热性与尚未实证的观察者仪器。

## 7. 检查与结论边界

[check.py](check.py)的严格部分核(12)、(15)、(18)及局部null张量收缩的精确分数例；连续场论结论由上述证明及已指明的成熟定理承担。固定RK/Gauss仅诊断 $r,\eta,B,M_\eta,C,E_g$，网格差不充当严格误差。

本轮得到：同一明确bulk参考、原最小物质平均源、独立经典几何与null模响应之间的有限共同证书。绝对零参考源只在已固定处方及非紧横向覆盖成立。未得到全局有限资源、同一真实仪器的温标、稳定参考的形成恢复、所有记录/控制源或全目标完成。具体停止线见[NEXT.md](NEXT.md)。
