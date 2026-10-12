# 1138资格补全：同一壳层的全 Fock 图界与读取域

2026-10-11。接续 [未完成半稿](local1138_kernel_preflight.md)，保留原稿。本文件是候选1138的资格引理，**不登记1138完成，也不独立签来源参考迁移、实际信号或区域动力学验收**。固定原 $M=10^{80}$、$A=10^{34}$、$\epsilon=10^{-5}$、正常Y、完整Dirac海及AB/rail副本、标量Fock、D/CP/Z/O/Q/R、$b=0$ 和 $0\le t\le4$。

认知动机是同一过程的区域资源账；区域、平滑核和有效读取是额外操作建模输入。以下为固定模型的解析引理；代码只复算常数及有限诊断。本引理未生成新认知公理、逐点应力张量、完整场论或自治探测器。

## 1. 算符与表象

采用 [动力学预审](local1138_dynamics_preflight.md) 和 [范围预审](local1138_scope_preflight.md) 的同一裸量。LLP只酉变换完整正常输入，不重置Y/K/源的关联。在此表象内区相对实际Y定位，写成中心零的自由区域密度；可酉拉回实际 $E(Y)$。

$$
L=3M+H_D+1+H_F+H_\phi+H_Z+H_O+V,\quad
A_0=H_F+H_\phi+H_Z+H_O+|K|+1.
$$

内区取自由正常序密度的空间/时间平滑，外区用 $1-\chi$。区域总能源及去恒等项量为

$$
E_w=3M+\mathcal E_w,\qquad
\mathcal E_w=E_{F,w}+E_{\phi,w}+H_Z+H_O+H_D+1+V.
$$

内外求和为L；实际H的 $R_{rec}+g$ 另列。物理能源未删除3M，只在差值和图估计中单列。裸量沿L作短时标定，实际态沿H；不是H历史平均。$D(H)=D(A_0)$及图等价复用动力学预审。

## 2. 固定空间和时间核

令 $P(u)=35u^4-84u^5+70u^6-20u^7$ 于 $[0,1]$，两端延拓为0、1。先取 $g(r)=1-P((r-6.01)/0.98)$，再在实变量r上与非负归一 $C_c^\infty(-0.01,0.01)$ bump卷积，最后令 $\chi(x)=g_{sm}(|x|)$。它在 $B_6$ 为1，$B_7^c$为0；原点邻域常数，故是三维径向光滑紧支函数。

$$
P'=140u^3(1-u)^3,\quad
\|P'\|_\infty=35/16,\quad
\|P''\|_\infty=84/(5\sqrt5),\quad
\|P'''\|_\infty=105/2.
$$

宽度49/50缩放及正卷积保持

$$|\chi'|\le3,\quad |\chi''|\le8,\quad |\chi'''|\le84.$$

傅里叶变换为酉归一。记 $M_j=\int|q|^j|\widehat\chi(q)|^2d^3q$。球/壳体积、Plancherel和Cauchy给

$$M_0<1437,\quad M_2<4797,\quad M_1\le\sqrt{M_0M_2}<3117.$$

令 $V_{sh}=4\pi(7^3-6^3)/3$。径向拉普拉斯给

$$
\|\Delta\chi\|_2^2\le81V_{sh},\quad
\|\nabla\Delta\chi\|_2^2\le(84+16/6+6/36)^2V_{sh}.
$$

定义 $F_j=(2\pi)^{-3/2}\int|q|^j|\widehat\chi(q)|d^3q$。加权Cauchy与积分 $\int(1+q^2)^{-2}=\pi^2$、$\int q^2(1+q^2)^{-3}=3\pi^2/4$ 给

$$F_0^2\le\frac{M_0+2M_2+81V_{sh}}{8\pi}<50^2,$$

$$F_1^2\le\frac{3[M_0+3M_2+243V_{sh}+(84+16/6+6/36)^2V_{sh}]}{32\pi}<420^2.$$

较紧数值为46.405和352.176，证书只用50和420。更高Fj有限，因χ光滑紧支。

取 $u_a=(2a)^{-1}1_{[-a,a]}$、$a=1/32$ 及同支撑非负归一光滑bump $b_a$，固定

$$w=u_a*u_a*u_a*b_a,\quad \operatorname{supp}w\subset[-1/8,1/8],\quad \int w=1.$$

w非负且光滑，特征傅里叶变换满足

$$|\widetilde w(s)|\le\min\{1,(32/|s|)^3\}.$$

令 $q=k+l,p=(k-l)/2$，Jacobian为1，质量零或质量1都有 $\omega_k+\omega_l\ge2|p|$，所以

$$
I_n=\int |p|^n\min\{1,(16/|p|)^6\}d^3p
=4\pi16^{3+n}\left(\frac1{3+n}+\frac1{3-n}\right),\quad n=0,1,2.
$$

三次sinc尾只用于数值预算；强交换子的额外频率权使用完整光滑bump的快速衰减。

## 3. 从单粒子到全 Fock 的图引理

若正h及对称A满足图界 $\|A\psi\|\le c_1\|h\psi\|$ 和形式界 $-c_0h\le A\le c_0h$，则对称或反对称Fock芯上

$$\|d\Gamma(A)\Psi\|\le\max(c_0,c_1)\|d\Gamma(h)\Psi\|.$$

证明：对角项 $A_i^2\le c_1^2h_i^2$。异粒子项通过 $h^{-1/2}Ah^{-1/2}$ 的有界自伴形式闭包作张量积，得 $A_iA_j\le c_0^2h_ih_j$。展开 $(\sum_iA_i)^2$ 后取每个粒子数直和并闭包。零能处可先正则化h+δ。这里不要求A与h对易；仅单粒子形式界不足以替代该图证明。

## 4. 标量 CCR 的所有数项与配对项

正常序密度 $\tfrac12:(\pi^2+|\nabla\phi|^2):$ 给

$$E_{\phi,w}=d\Gamma(A_\chi)+Q+Q^\dagger,\quad Q=\tfrac12\int B(k,l)a_ka_l\,dk\,dl.$$

数核系数为 $\sqrt{\omega_k\omega_l}(1+\widehat k\cdot\widehat l)/2$，配对B同绝对系数（Q定义已含1/2）。数/配对核分别乘 $(2\pi)^{-3/2}\widehat\chi(k-l)\widetilde w(\omega_k-\omega_l)$ 和 $(2\pi)^{-3/2}\widehat\chi(k+l)\widetilde w(\omega_k+\omega_l)$。

数项形式界为F0。对 $A_\chi h^{-1}$ 用Schur权 $|k|^{-1/2}$，列界F0，行界

$$|k|\int\frac{f(k-l)}{|l|}dl\le\int f(q)dq=F_0,\qquad f=(2\pi)^{-3/2}|\widehat\chi|.$$

因f径向非负，球壳上 $1/|k-q|$ 平均为 $1/\max(|k|,|q|)$。因此处理了质量零红外，未把1/|l|当有界。第3节给全Fock数项系数50。

定义 $\kappa=\|B/\sqrt{\omega_k\omega_l}\|_{HS}$、$\eta=\|B/\sqrt{\omega_k}\|_{HS}$、$\beta=\|B\|_{HS}$。p/q变换与 $|l|\le|p|+|q|/2$ 给

$$\kappa^2\le I_0M_0/(2\pi)^3<447^2,$$

$$\eta^2\le(I_1M_0+I_0M_1/2)/(2\pi)^3<4\cdot10^6,$$

$$\beta^2\le(I_2M_0+I_1M_1+I_0M_2/4)/(2\pi)^3<10^8.$$

带权向量Cauchy给 $\|Q\Psi\|\le(\kappa/2)\|H_\phi\Psi\|$。CCR完整重排为

$$[Q,Q^\dagger]=d\Gamma(B^\dagger B)+\beta^2/2,\qquad d\Gamma(B^\dagger B)\le\eta^2H_\phi.$$

这里B对称，行/列加权Hilbert–Schmidt范数相同；公式按 $Q=\tfrac12\sum B_{ij}a_i a_j$ 的复核约定写成B†B。因而

$$\|(Q+Q^\dagger)\Psi\|\le\kappa\|H_\phi\Psi\|+\eta\|\sqrt{H_\phi}\Psi\|+\beta\|\Psi\|/\sqrt2.$$

用谱不等式 $\eta\sqrt{x}\le50x+\eta^2/200$，有

$$\boxed{\|E_{\phi,w}\Psi\|\le547\|H_\phi\Psi\|+30000\|\Psi\|.}$$

从全粒子芯按图连续性延拓至D(Hφ)，不限制玻色粒子数，不丢真空涨落。

## 5. 完整 Dirac 粒子/空穴及四个副本

原空间 $L^2(\mathbb R^3,\mathbb C^4)\otimes\mathbb C^2_{AB}\otimes\mathbb C^2_{rail}$ 含四个Dirac副本。每副本正负谱各2维，全部正/负谱各8维。$h_+=\sqrt{k^2+1}$，$H_F\ge N$。

对称密度 $\tfrac12:[\psi^\dagger h\psi+(h\psi)^\dagger\psi]:$ 的同号谱数核系数为 $(\omega_k+\omega_l)/2$，异号配对系数为 $(\omega_k-\omega_l)/2$。投影/自旋矩阵算符范数≤1。$|\omega_k-\omega_l|\le|k-l|$、ω≥1给数核图Schur界及形式界同时不大于 $F_0+F_1/2<260$。全部粒子空穴数项是正能源直和上的单粒子算符，第3节提升后仍为260，不按副本数重复乘算符范数。

配对q=k+l时 $|\omega_k-\omega_l|\le|q|$。四个副本各有一个2×2粒洞块，共16个非零矩阵条目，每条模长≤1；配对系数1/2的平方给

$$\beta_F^2=\|B_F\|_{HS}^2\le4I_0M_2/(2\pi)^3<3\cdot10^6.$$

这明确计入原16维Cauchy标签，未把一般8×8矩阵无理由只算16条目。

SVD给 $Q_F=\sum_j s_jd_jb_j$。有限截断的向量Cauchy和CAR给

$$\|Q_F\Psi\|\le\beta_F\|\sqrt N\Psi\|,\quad [Q_F,Q_F^\dagger]=\sum_js_j^2(1-n_{b,j}-n_{d,j})\le\beta_F^2I.$$

所以 $\|(Q_F+Q_F^\dagger)\Psi\|\le2\beta_F\|\sqrt{N+1}\Psi\|$。用 $2\beta_F\sqrt{x+1}\le x+1+\beta_F^2$ 及N≤HF；SVD尾的HS范数趋零，允许在D(HF)取极限。海配对本身未被误称为全Fock有界。

$$\boxed{\|E_{F,w}\Psi\|\le261\|H_F\Psi\|+3000001\|\Psi\|.}$$

## 6. 系数650与平均壳流

令 $B_+=H_F+H_\phi+H_Z+H_O$。正强对易自由项和Cauchy给

$$547\|H_\phi\Psi\|+261\|H_F\Psi\|+\|(H_Z+H_O)\Psi\|\le\sqrt{367331}\|B_+\Psi\|<607\|B_+\Psi\|.$$

留裕量取650；冻结 $\|H_D\|<27$、$v=16(11/7)C_*+558/7$、$C_*=33536000000000000000$，全部有限余项小于 $10^{60}$，得

$$\boxed{\|\mathcal E_w\Psi\|\le650\|B_+\Psi\|+10^{60}\|\Psi\|.}$$

状态资源账必须控制此B+或能严格控制它的正自由账，不可把带符号HD混入后直接假定正性。

平均壳流由独立自由连续方程定义。$a=25/8,b=29/8,h=1/2,t_*=27/8$，$\overline J=h^{-1}\int_{-h/2}^{h/2}\alpha_L^s(J_w)ds$。外区光锥与V分离时，Jbar等于两个**自由演化的区域自由能源**之差除以h；HZ/HO/HD/V不在该自由差内。自由传播保持HF+Hφ，第4、5节给

$$\boxed{\|\overline J\Psi\|\le2600\|B_+\Psi\|+4\cdot10^{60}\|\Psi\|.}$$

因子4为2/h，不是未经证明的 $\|J_w\Psi\|\le\|w'\|_1\|E_w\Psi\|$。未经平均的流会改变频率尾，须另证。

## 7. 强域、闭包与有限后资源

共同芯C取有限标量/Dirac粒子空穴、光滑紧动量且避开质量零原点的波包，K光滑紧支且避开K=0、Z光滑紧支、有限球谐转子向量及全部有限D/CP。它是A0²的芯且A0C⊂C，任意被动参考可张量。设O为ℰw或Jbar。

自由数核的 $i[H_{free},O]$ 多出 $i(\omega_k-\omega_l)$，由|k−l|控制；第3节及同Schur证明给有限图常数（Dirac再用有限F2）。配对核多出ωk+ωl；完整光滑时核对任意固定n满足 $|\widetilde w(s)|\le C_n(1+|s|)^{-n}$，故所需加权HS积分仍有限。再多一个能源权证OC⊂D(A0)。有限时平移/平均不改变这些有限性。

真实V的CAR因子由原光滑模式及固定ρε平滑，$[H_F,c(f)]=-c(hf)$和乘积规则给有界导数。对实平滑f和实际门B，强域公式为

$$[H_\phi,B(\phi(f)/A)]=-\frac{i}{2A}\{\pi(f),B'(\phi(f)/A)\}.$$

Segal界、B′/B″有界及固定f给 $\sqrt{H_\phi+1}$ 图界。可先用Weyl共轭Hφ的精确一次场加常数公式，再积分B的傅里叶表示验证域；阶跃门S用其光滑紧支导数表示 $S(X)=S(0)+\int\widehat{S'}(s)(e^{isX}-1)/(is)ds$（统一傅里叶归一），无需S本身可积。更高有限阶Weyl公式保证芯向量上所需乘积有定义。空间积分在固定纤维作图Bochner积分。$[H_O,V]$ 是有界一阶角系数乘Lj，加有界二阶角系数，1137冻结角导数给A0图界。b=0时HZ与V强对易；LLP纤维K与V强对易。

因此存在有限常数c′O，芯上有**强**身份

$$A_0O\psi=OA_0\psi+C_O\psi,\qquad\|C_O\psi\|\le c'_O\|A_0\psi\|.$$

以A0²图芯逼近，右侧收敛，A0闭性将身份延拓至D(A0²)。O的A0图界及强身份给

$$\|O\psi\|\le c\|A_0^2\psi\|,\qquad|\langle O\psi,A_0^2\psi\rangle-\langle A_0^2\psi,O\psi\rangle|\le2c'_O\|A_0\psi\|^2.$$

Nelson交换子准则取比较算符A0²，得O在C上本质自伴。仅A0的弱形式/半域结论不够。

再取 $R_\varepsilon=(1+\varepsilon A_0^2)^{-1}$、$O_\varepsilon=R_\varepsilon O R_\varepsilon$，后者有界自伴，强交换子为 $R_\varepsilon C_OR_\varepsilon$。图范数Grönwall给

$$\|A_0e^{isO_\varepsilon}\psi\|\le e^{c'_O|s|}\|A_0\psi\|.$$

Oε在O芯上收敛，故强预解式和酉群强收敛；A0图球弱闭性传递一致界。因此 $e^{isO}$ 保持D(A0)=D(H)，有限个这类酉群的线性组合也保持。

例如 $K_\pm=[\cos(\tau\overline J/2)\pm\sin(\tau\overline J/2)]/\sqrt2$ 未归一后态具有有限H图资源。不授予所有有界谱函数相同域性质，不给稀有归一后选择支路统一能源界，也不把有限后资源当成自治探测器做功已闭合。

## 7.1 外侧量与短时强传播身份

非紧支外侧量直接定义为

$$E_{out,w}=H_F+H_\phi-E_{F,w}-E_{\phi,w}.$$

它在D(A0)有上述图定义。这里没有假定无限空间截断的局域能源以强算符拓扑自动收敛，也没有使用未作时间平滑的Wick平方作为处处有定义的算符。令α0为自由传播。以点分裂Wick双线性和全部正/负Dirac CAR先在共同测试芯上检验：V的所有场支撑在 $B_{5+\epsilon}$，自由外区测试点经 $|s|+1/8\le3/8$ 回溯仍在半径 $6-3/8=45/8$ 之外。因此外区的标量/偶CAR表达式与V对易；正常序减项为常数，不改变交换子。对全能源项与区域项先作此交换子相减，外侧非紧支部分只承担这个局域恒等式，而不要求截断能源本身强收敛。

先作真实C∞时平滑，其核为前节的图有界核；点分裂极限在该图估计/测试矩阵元中收敛。利用V有界且保持D(A0)，由芯逼近得到

$$[V,\alpha_0^s(E_{out,w})]\psi=0,\quad\psi\in D(A_0),\quad |s|\le1/4.$$

此处强交换子中的两个乘积有定义，不能只用几何支撑代替域论证。为比较L和自由传播，仅在 $D(A_0)$ 上采用双端矩阵元的弱Duhamel论证，不要求L保持 $D(A_0^2)$。写 $U(t)=e^{-itL}$、$O_u=\alpha_0^u(E_{out,w})$。L群保持 $D(A_0)$ 且在A0图范数中连续；自由核及其一阶时间导数均为A0相对有界，并在固定短窗一致。对 $\psi,\varphi\in D(A_0)$ 考察

$$f(u)=\langle U(s-u)\psi,O_uU(s-u)\varphi\rangle.$$

两条L轨道在Hilbert范数可微、在A0图范数连续，因此差商的线性项有定义。双轨道交叉项除以步长δ，其绝对值由 $\|\Delta U\psi\|/|\delta|$ 的有界性与 $\|O_u\Delta U\varphi\|\to0$ 控制而趋零；涉及 $\Delta O_u$ 的交叉项同样由一致图界和图连续性控制。由此无需写 $LO_u\varphi$ 这种未获域许可的乘积，即可逐项求矩阵元导数。自由导数与L的对易型先在共同芯相消，再由A0图连续性延拓到两端的 $D(A_0)$，剩余

$$f'(u)=-i\langle U(s-u)\psi,[V,O_u]U(s-u)\varphi\rangle=0.$$

积分得到所有此类双端矩阵元的弱身份。两侧都已是共同域 $D(A_0)$ 上的算符，利用该域稠密即可提升为强身份。因此

$$\alpha_L^s(E_{out,w})=\alpha_0^s(E_{out,w}),\qquad |s|\le1/4.$$

再由 $E_w=L-E_{out,w}$ 及自由连续方程，获得第6节差商与独立自由壳流的一致。空间/时间平滑、完整正负谱与域控制在这里共同使用；没有只用正频场作虚假的精确支撑推断。

## 8. 复算与尚待验收

[代码](local1138_kernel_completion.py) 默认只读核 [结果](local1138_kernel_completion.json)，显式 --write-results 才写结果。Machin有理交错余项围住π，精确核算整数界和冻结来源SHA；非对易数核张量、保留全部目标输出的CCR输入子空间及完整有限CAR提供独立诊断。连续域证明由本文承担。

资格结论为cE650/dE1e60、平均壳流2600/4e60、明确裸量读取的强资源域。实际1138仍须同一正常Y/源/全材料参考图界、正信号余量、三点有限平衡误差、读取精度与完整扰动预算及独审。不能闭合就停在资格引理，不调整M/A/窗口掩盖失败。
