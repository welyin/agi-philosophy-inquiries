# 1287：恒定局部混合产生空间量子记录，并保全部物质来源

2026-10-18。正式状态以独立审阅和验收记录为准。本证明只实现[既定准入](../_shared/source_scope_and_round1287_admission.md)中的给定背景量子过程、因果接续与几何来源；新公理为零。

## 1. 固定对象、参数与旧结果差分

采用四维Minkowski背景、签名$(-+++)$及$\hbar=c=1$。质量单位固定后，全部长度与时间使用其倒数单位；不把质量1设为Planck质量，也不从数值大小领取弱引力误差。两个实标量组成$\varphi=(\phi_1,\phi_2)^T$，另有两个简并、规范中性的Dirac味$q=0,1$。同一个局部协变作用为

$$
S=\int\sqrt{-g}\left[-\frac12\nabla_a\varphi^T\nabla^a\varphi-\frac12\varphi^TM^2\varphi+\sum_{q=0}^1\left(\frac i2\bar\psi_q\gamma^a\overleftrightarrow\nabla_a\psi_q-3\bar\psi_q\psi_q\right)\right],\quad
M^2=\begin{pmatrix}5/2&-3/2\\-3/2&5/2\end{pmatrix}.
\tag{1}
$$

Dirac项复用[942 §§3—5](../../archive_935_955/research_note_942.md)的正能等距，但显式转换签名：在本$(-+++)$正交标架取$\{\gamma^a,\gamma^b\}=-2\eta^{ab}I$、$\bar\psi=\psi^\dagger\gamma^0$，沿用942的矩阵，故$\beta_D=\gamma^0$、$\alpha_i=\gamma^0\gamma^i$及$H_D=\boldsymbol\alpha\cdot p+3\beta_D$不变。以本度规下降$\gamma_a$后，在壳Hilbert式为

$$
T^D_{ab}=-\frac i4\bar\psi\left(\gamma_a\overleftrightarrow\nabla_b+\gamma_b\overleftrightarrow\nabla_a\right)\psi,\qquad
T^D_{00}=\frac i2\psi^\dagger\overleftrightarrow\partial_t\psi.
$$

因此正能包给正总能源；不将942的$(+---)$下降指标应力号直接搬到本式。两个KG本征场$\phi_+=(\phi_1+\phi_2)/\sqrt2$、$\phi_-=(\phi_1-\phi_2)/\sqrt2$的质量分别为1、2。质量混合恒定，没有外加时变开关、远处无来源控制或丢弃的相互作用项。

固定

$$
a=1/32,\quad L=1/4,\quad I=[7/16,1/2],\quad\zeta=1/256,
$$

$$
b_a(x)=\begin{cases}\exp[-1/(1-|x|^2/a^2)],&|x|<a,\\0,&|x|\ge a,\end{cases}
\quad Z_a=\int b_a,\quad f_A=b_a/Z_a,\quad h_B(x)=f_A(x-Le_1).
\tag{2}
$$

这是预先固定的正常空间准备及接收测试，无参数或积分扫描。$A$、$B$是这两个不交球。下文“接收者”专指第二KG场的指定代数，不是B中所有物种。Dirac正能包可以在B有尾，不能因此声称B全部物理资料起初均不含输入。

[524](../../archive_467_530/research_note_524.md)已给局部线性场散射；942已给含能的内部味载体；[1124后查重 §3](../../archive_1120_1158/_shared/1120_1156/after1124_execution_audit.md)明确禁止只换味标签重复。新增合同是：同一完整正常码在指定接收代数起初与输入及伙伴独立，随后在固定有限空间窗形成严格可分记录，同时全部局域来源对未知码相同。它不沿用[1155](../../archive_1120_1158/research_note_1155.md)的外开关或[1160](../../archive_1159_1160/_shared/field_bridge_after_1160.md)的机械读取权限。

## 2. 正常参考、含能标签及未知输入

令$\Omega$是质量1、2的Minkowski自由真空经上述固定正交变换得到的耦合KG真空。它为正常Gaussian Hadamard参考，通常不在原$\phi_1/\phi_2$划分下因子化。Dirac采用同一平直真空、固定正能spin与归一动量包

$$
g(p)=\pi^{-3/4}e^{-|p|^2/2},\quad
|q\rangle_{\rm tag}=b_q^\dagger(g)|\Omega_D\rangle,\quad
\langle p^2\rangle=3/2,\quad E_{\rm tag}\le3+\langle p^2\rangle/6=13/4.
\tag{3}
$$

正能旋量通过942的等距归一，故单粒子$L^2$范数为1；两个味严格正交、相同质量和波包，全部标签应力为味单态。它不是无能源全局比特。正能包并非紧空间支撑；此处不领取标签的紧区自主制备。

取$W_\pm=\exp[\pm i\phi_1(0,f_A)]$。约定使$W_+^\dagger\pi_1(x)W_+=\pi_1(x)+f_A(x)$，初经典位移为零，第一场速度为$\pm f_A$，第二场速度为零。定义

$$
J|0\rangle=|0\rangle_{\rm tag}\otimes W_+\Omega,\qquad
J|1\rangle=|1\rangle_{\rm tag}\otimes W_-\Omega.
\tag{4}
$$

标签正交给$J^\dagger J=I_2$，不要求两个相干态正交。紧支光滑Cauchy数据使位移正常且有限能；各固定能量矩有限。对任意未知$\rho_{QR}$，完整过程是$(U_tJ\otimes I_R)\rho_{QR}(J^\dagger U_t^\dagger\otimes I_R)$，保所有场、标签及伙伴，不做后选择。$J$是公开的含能初态输入，未证明它由本次自治等待制造。

对第二场**完整初片代数**$\mathcal A_{2,0}$，$[W_\pm,B]=0$。对任意$B\in\mathcal A_{2,0}$及有界$C_R$，标签内积杀去输入矩阵的非对角项，而两对角项均给$\omega(B)$，所以

$$
\omega_{\rm in}(BC_R)=\omega(B)\operatorname{Tr}(\rho_R C_R).
\tag{5}
$$

这是接收代数与R的联合乘积身份，不借真空因子化。未来接收场只需记录指定逻辑基的区别，不声称接收者复制完整未知量子态。

## 3. 有限空间与时间窗内的精确信号

$K_m(t)=\sin(t\sqrt{-\Delta+m^2})/\sqrt{-\Delta+m^2}$是速度数据传播算符。其三维核为

$$
K_m(t,r)=\frac{\delta(t-r)}{4\pi r}-\mathbf1_{t>r}\frac{mJ_1(m\sqrt{t^2-r^2})}{4\pi\sqrt{t^2-r^2}},\quad t>0.
\tag{6}
$$

该表达直接由KG基本解或$J_1$级数代入径向方程及初值得到。负混合约定给$\phi_{2,+}=(K_1-K_2)f_A/2$。两质量的光锥delta项抵消；内部尾核为

$$
D(t,r)=\frac{2J_1(2s)-J_1(s)}{8\pi s},\qquad s=\sqrt{t^2-r^2}.
\tag{7}
$$

当$x\in B,y\in A,t\in I$，$3/16\le|x-y|\le5/16<t$，且$3/32\le s^2\le55/256<1/4$。由

$$
J_1(z)=\sum_{n\ge0}\frac{(-1)^n(z/2)^{2n+1}}{n!(n+1)!}
$$

的递减交错项（[NIST DLMF 10.2.2](https://dlmf.nist.gov/10.2.E2)），有$2J_1(2s)/s\ge2-s^2$、$0\le J_1(s)/s\le1/2$、$2J_1(2s)/s\le2$。因此整个空间乘积与时间窗上

$$
\frac5{32\pi}\le D(t,|x-y|)\le\frac1{4\pi},\qquad
m(t):=\int h_B(x)\phi_{2,+}(t,x)\,d^3x\in\left[\frac5{32\pi},\frac1{4\pi}\right].
\tag{8}
$$

负输入的均值为$-m(t)$。这是传播后的局域记录，未借直接读取标签。混合KG的主符号仍是度规光锥，Peierls/CCR传播支撑在因果锥内；当$0<t<3/16$时整个接收球与初始准备支持类空，其接收态保持参考限制。有限传播是本线性场的精确身份，不是齐次截断。

## 4. 有界效果及全窗严格概率界

取$F_B(t)=\phi_2(t,h_B)$、$E_B(t)=[I+\sin(\zeta F_B(t))]/2$。这是接收球局部von Neumann代数中的效果。参考的等时方差与t无关，按Parseval归一为

$$
v=\frac14\int\frac{d^3k}{(2\pi)^3}|\widehat h_B(k)|^2\left(\frac1{\omega_1(k)}+\frac1{\omega_2(k)}\right)\le\frac38\|h_B\|_2^2.
\tag{9}
$$

内半球上$b_a>1/4$，全域$\sup b_a<1/2$，故$Z_a>\pi a^3/24$、$\|h_B\|_2^2<12/(\pi a^3)<131072$，于是$v<49152$。Gaussian Weyl特征函数给

$$
p_\pm(t)=\frac12[1\pm e^{-\zeta^2v/2}\sin(\zeta m(t))].
\tag{10}
$$

由于$\zeta^2v/2<3/8$且$e^{3/8}\le1/(1-3/8)=8/5<2$，衰减因子$>1/2$；$0<\zeta m<1$时$\sin(\zeta m)>\zeta m/2$。再用$\pi<4$得

$$
\boxed{p_+(t)-p_-(t)>\frac5{131072}\quad\text{对所有 }t\in[7/16,1/2].}
\tag{11}
$$

任何事前给定非负归一$\chi\in C_c^\infty((7/16,1/2))$定义$E_B^\chi=\int\chi(t)E_B(t)dt$，其强／弱算符积分仍是同一有限worldtube代数中的有界效果，保同一下界。它不要求理想瞬时末读。这里没有制造相应读取器，也没有授予新的控制或读取源。

若调用的两个预测各有独立事件误差$\epsilon_\pm$，差下界减去$\epsilon_++\epsilon_-$；例如各不超过$10^{-6}$仍有正余量。此为公开误差合同，不是本轮已校准设备。迹距离控制有界效果，不能代替下节无界应力的独立来源预算。

## 5. 同一局域完整来源、有限账与未知伙伴

标量Hilbert来源取式(1)的最小耦合变分：

$$
T^{\rm KG}_{ab}=\nabla_a\varphi^T\nabla_b\varphi-\frac12g_{ab}\left[\nabla_c\varphi^T\nabla^c\varphi+\varphi^TM^2\varphi\right],\quad
T^{\rm full}_{ab}=T^{\rm KG}_{ab}+\sum_qT^{D,q}_{ab}.
\tag{12}
$$

量子场使用一个固定、局域协变、保守的Hadamard/Wick处方；平直共同真空应力为零是该处方公开的归一选择，而不是省略真空。标量可在固定质量本征基逐场采用[Moretti，gr-qc/0109048，摘要及正文的保守应力构造](https://arxiv.org/abs/gr-qc/0109048)，再正交转回；Dirac沿942的同一平直正常序来源。曲率反项、宇宙常数及非最小耦合不是由本记录实验选定；本轮固定最小耦合及平直真空归一，不把加改进项后的另一来源冒认同一式(12)。

参考一次均值为零，两个相干位移相反，二次应力的参考相减均值相同。场应力的码间非对角元由标签正交杀掉；标签应力的非对角元由味单态杀掉，而两个同包对角元相同。故对任意实光滑紧支四维张量测试$w^{ab}$，有二次型身份

$$
J^\dagger U_t^\dagger T^{\rm full}(w)U_tJ=t_w(t)I_2.
\tag{13}
$$

等价地，可将时间演化并入w的Heisenberg来源。因而任意未知Q/R，包括其相干项，都给相同无条件平均局域来源；与R的任意有界测试联合也分离。这里不声称所有应力高阶涨落、选定结果分支来源或测量后几何均相同。

KG壳上$\nabla^aT^{\rm KG}_{ab}=(\Box\varphi-M^2\varphi)^T\nabla_b\varphi=0$，Dirac来源同样守恒。混合项是此身份的一部分；只留对角质量时一般产生非零缺陷。全部源守恒不是只检查总H。

初位移为零使混合势初均值为零，但其后不得删除。相干能量为

$$
E_{\rm KG}=\tfrac12\|f_A\|_2^2<65536,\qquad E_{\rm full}<65536+13/4=262157/4.
\tag{14}
$$

两个符号同能，任意未知码同源。Hadamard参考加光滑位移、正常Schwartz标签使本合同的有限涂抹来源存在；各固定能量矩有限。

更具体地，正质量矩阵保证经典每个$|T^{\rm KG}_{\mu\nu}|\le T^{\rm KG}_{00}$，其空间积分不超过$E_{\rm KG}$。归一Dirac正能波函数满足$\|\partial_0\psi\|_2^2=\langle9+p^2\rangle=21/2$，各空间导数不大于此值。由Dirac对称双线性、Cauchy–Schwarz和固定正交标架，$\int|t^{D}_{\mu\nu}|\,d^3x\le\sqrt{21/2}<4$。所以

$$
|t_w|<65540\int dt\sum_{\mu,\nu}\|w^{\mu\nu}(t,\cdot)\|_\infty<\infty.
\tag{15}
$$

这是独立的完整平均源投影界；没有将记录概率误差误用来控制应力，也未将量子局部能量普遍假设为正。公开定态基线如另选非零，则须显式加其测试贡献。

## 6. 旧比较角色在同一作用中的限定落点

本例$M^2=I+B$，$B=(3/2)\left(\begin{smallmatrix}1&-1\\-1&1\end{smallmatrix}\right)$；势能包含独立裸质量及正比较势$3(\phi_1-\phi_2)^2/4$。只对事前单列的比较能源，设一般$B=\left(\begin{smallmatrix}a&-b\\-b&d\end{smallmatrix}\right)\succeq0$。固定同一Hadamard协方差和零均值参考，施加共同配置位移$\theta u(x)(1,1)$，则

$$
\Delta E_{\rm comp}=\frac{\theta^2}{2}(a+d-2b)\int u^2.
\tag{16}
$$

比较 $\theta=0$ 与任一 $\theta\ne0$ 的校准，比较成本增量为零，给$(1,1)B(1,1)^T=0$；PSD遂给$B(1,1)^T=0$，即$a=b=d=\kappa\ge0$。这是[1282](../1282/proof.md)已有正二次共同校准选择的直接例化，不另计定理，也不称完整成本为零。裸质量1、$\kappa=3/2$、场表、背景及实际校准权限仍是输入；本轮不借此领取SU(2)群、E-PHASE动能或所有参数的认知生成。

## 7. 删除、范围与验收责任

删去混合后，第二场动力学自治。对任意合格、起初接收场与Q/R独立的准备，接收演化不能产生输入信息。若展示具体删除模型，须使用该模型自身的Hadamard参考；不在原耦合真空上突然换质量再隐去quench来源或UV问题。此删除只验证本局域机制，不是同一次实验的免费开关。

本轮的新证据是正常未知码、指定空间因果记录、整个有限窗概率间距与完整含能来源的联合交付；不是重新证明一般Weyl公式、Schur引理或Noether定理。1286仍保持其齐次有效范围，本例是另一个公开局部场实现。

给定背景上可耦合几何的来源，不等于Minkowski满足此非零来源的自引力方程；本轮没有给新的半经典响应解、噪声反作用预算、实际末仪器、稳定经典记录、热温标或整个阶段完成。参考是比较对象，未调用的制备工厂不是隐含免费装置；初态的有限材料资源已显式计入。

[作者复算](check.py)只做一组精确有理证书：质量/投影、Bessel交错尾与几何窗、归一与噪声/能量界、有限jet的完整守恒及漏混合缺陷、同源味矩阵选择与校准零空间。连续量词由以上解析证明承担；不以有限样本替代Hadamard、因果或未知伙伴定理。
