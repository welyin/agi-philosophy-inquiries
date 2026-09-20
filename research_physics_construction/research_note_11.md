# 构造路线第11轮：径向线性稳定性、移动表面与动态钟记录

日期：2026-09-17。C11。只论证充分性，必要性只登记。

**结果：保持C10同一个EOS和中心参数，得到径向线性本征问题，构造与本征数值无关的正能量下界，并用精确有理比较证明所有允许径向线性模满足omega²>0.015。前3个数值模与独立射击相容；最低模的质量、粒子数与钟率扰动保留移动表面项，并接回时变量子探针。**

这是给定Einstein–流体合同下的条件性径向稳定结果，不是非线性稳定或完整量子引力。范围见[C11约定](C11_CONTRACT.md)。

## 1 保持原来的物理问题

不换源或EOS：rho_s=.01、s=1/3、p_c=.0005、G=.1，R约4.13486431、M约3.13038623，背景来自C10。使用 $ds^2=e^{2\nu}dt^2-e^{2\lambda}dr^2-r^2d\Omega^2$，B=e^(2lambda)。时间仍按无穷远静止钟归一化，频率不是已标定的Hz。

扰动为球对称、绝热、无耗散。位移xi对应面积半径变化；Euler变化delta在固定r定义，Lagrange变化 $\Delta f=\delta f+\xi f'$。模式取 $\xi(r,t)=\xi(r)e^{i\omega t}$，实际小振幅在后面的读数任务中另行声明。

## 2 从约束与流体运动得到本征方程

线性质量与径向度规约束给

$$
\delta m=-4\pi r^2(\rho+p)\xi,\qquad
\delta\lambda=\frac{GB}{r}\delta m.
$$

定义 $\zeta=r^2e^{-\nu}\xi$。粒子数守恒和绝热EOS给

$$
\frac{\Delta n}{n}=-\frac{e^\nu}{r^2}\zeta',\quad
\Delta\rho=(\rho+p)\frac{\Delta n}{n},\quad
\Delta p=s\Delta\rho=-s(\rho+p)\frac{e^\nu}{r^2}\zeta'.
$$

径向Einstein约束线性化为

$$
\delta\nu'=\frac{G(\delta m+4\pi r^3\delta p)}{r(r-2Gm)}
+\frac{2G\delta m}{r-2Gm}\nu'.
$$

流体Euler方程的径向线性部分为

$$
-\omega^2(\rho+p)Be^{-2\nu}\xi
+\delta p'+(\delta\rho+\delta p)\nu'
+(\rho+p)\delta\nu'=0.
$$

代入 $\delta p=\Delta p-\xi p'$、$\delta\rho=\delta p/s$及上述约束，使用背景TOV关系，得到自伴方程

$$
\boxed{(P\zeta')'+(Q+\omega^2W)\zeta=0,}
$$

$$
P=\frac{s(\rho+p)e^{\lambda+3\nu}}{r^2},\qquad
W=\frac{(\rho+p)e^{3\lambda+\nu}}{r^2},
$$

$$
Q=\frac{e^{\lambda+3\nu}}{r^2}
\left[\frac{p'^2}{\rho+p}-\frac{4p'}r-8\pi G e^{2\lambda}p(\rho+p)\right].
$$

这是既有相对论径向扰动结构在当前EOS上的应用，不宣称原创定理，也不重新研究引力作用量的必要性。代码还从重建场独立差分核对Euler方程和 $\delta m'=4\pi r^2\delta\rho$，不是只检验矩阵自身的本征方程。

## 3 自束缚表面不是稀薄表面

中心位移xi=O(r)，所以zeta=O(r³)。表面为随流体移动的自由边界，要求 $\Delta p(R)=0$。因为 $s(\rho+p)\to s\rho_s>0$，得到

$$
\zeta'(R)=0.
$$

不能因p(R)=0就把Gamma p当成0，其中 $\Gamma p=s(\rho+p)$。C10的径向接口在r=R返回外侧真空，本轮显式覆盖为内侧rho_s和对应p'，不更改C10文件。

中心条件与自由表面使分部积分边界项消失，频率平方满足Rayleigh比

$$
\omega^2=\frac{\int_0^R(P\zeta'^2-Q\zeta^2)dr}
{\int_0^R W\zeta^2dr}.
$$

惯性W正。有限维Ritz本征值给上方近似，故不能只看它们为正就宣布全部连续径向模稳定。

## 4 只用输入参数构造充分下界

对任意本合同的正则TOV平衡，p随r下降，rho_s≤rho≤rho_c。因为

$$
-p'\ge\frac{4\pi G\rho_s^2}{3}r,
$$

可得

$$
R^2\le D:=\frac{3p_c}{2\pi G\rho_s^2}=\frac{75}{\pi},\qquad
\frac{2Gm(r)}r\le C:=\frac{8\pi G\rho_cD}{3}=\frac{23}{100}.
$$

因此B≤(1-C)^(-1)。由C10的焓第一积分，整个内部ell有参数下界

$$
\ell\ge L:=\sqrt{1-C}\left[1+\frac{(1+s)p_c}{s\rho_s}\right]^{-s/(1+s)}
\approx0.83839765,
$$

并且ell≤1。设

$$
P_*:=s\rho_sL^3,\quad
W^*:=\frac{\rho_c+p_c}{(1-C)^{3/2}},\quad
A:=\frac{4\pi G(\rho_c/3+p_c)}{1-C}.
$$

则 $P\ge P_*/r^2$、$W\le W^*/r^2$、$\nu'/r\le A$，以及

$$
\frac QW=e^{2\nu-2\lambda}\left[\nu'^2+\frac{4\nu'}r-8\pi GBp\right]
\le A^2D+4A.
$$

还需一个控制径向梯度的独立不等式。令 $f(r)=\sin(\pi r/R)-(\pi r/R)\cos(\pi r/R)$，它在(0,R]为正，中心为O(r³)，且f'(R)=0。由基态变换得

$$
\int_0^R\frac{\zeta'^2}{r^2}dr-
\frac{\pi^2}{R^2}\int_0^R\frac{\zeta^2}{r^2}dr
=\int_0^R\frac{f^2}{r^2}\left[(\zeta/f)'\right]^2dr\ge0.
$$

于是对整个允许的径向扰动空间，

$$
\boxed{\omega^2\ge\frac{\pi^2P_*}{DW^*}-(A^2D+4A).}
$$

代入精确输入后的较紧表达式浮点评估为0.01624503019。这个推导不依赖数值半径、数值本征值或节点扫描；条件是C10合同所述正则平衡存在，并且扰动遵循上述绝热径向方程。

## 5 用有理数固定稳定性符号

为了不把正负号交给浮点近似，使用已知 $157/50<\pi<22/7$并作保守界：

$$
D<24,\quad L>419/500,\quad W^*<9/500,\quad
A<71/10000,\quad \pi^2>197/20.
$$

例如L界可通过四次方比较 $(419/500)^4<(77/100)^2/(6/5)$核对；W界也只需平方比较。所有中间比较均由Python标准库Fraction用整数比完成。所得严格下界为

$$
\omega^2>
\frac{(197/20)(1/3)(1/100)(419/500)^3}{24(9/500)}
-24(71/10000)^2-4(71/10000)
=\frac{4897743463}{324000000000}>\frac3{200}>0.
$$

这给当前径向线性问题的充分稳定性结论，区别于小残差或正Ritz值的数值证据。它不是对C10平衡存在唯一性的独立证明，也不证明非径向、有限振幅或量子稳定性。

相应自伴系统的二次能量 $\tfrac12\int[W\dot\zeta^2+P\zeta'^2-Q\zeta^2]dr$为正并保持。这里采用的Sturm–Liouville能量规范化只用于线性分析；其与完整引力–物质二阶作用量的整体正则系数及ADM二阶能量，需要在量子化前明确，不能直接把数字当作已校准的物理能量。

## 6 两种独立数值路径

Ritz采用 $\zeta=(r/R)^3T_j(2(r/R)^2-1)$，积分用192点Gauss–Legendre。中心条件在基底中满足；表面导数条件是变分的自然边界，不通过硬加多余Dirichlet条件改变问题。

另令y=xi/r、z=Delta p，直接积分一阶方程

$$
y'=-\frac{3y+z/[s(\rho+p)]}{r}-\frac{p'}{\rho+p}y,
$$

$$
z'=y\left[\omega^2(\rho+p)Be^{-2\nu}r-4p'
+\frac{rp'^2}{\rho+p}-8\pi GB(\rho+p)pr\right]
+z\left[\frac{p'}{\rho+p}-4\pi GBr(\rho+p)\right].
$$

中心起点r=1e-5R，取y=1、z=-3s(rho_c+p_c)，再射击至z(R)=0。Brent求根的区间由Ritz谱中点给出；DOP853积分、背景容差和中心起点分别精化。中心常数初值只到领先正则阶，较小起点对照用于数值诊断。

| 模编号 | omega²（射击） | omega | 内部位移节点数 |
|---:|---:|---:|---:|
| 0 | 0.125641672822 | 0.354459691393 | 0 |
| 1 | 0.562006665347 | 0.749671038087 | 1 |
| 2 | 1.289430045495 | 1.135530732959 | 2 |

4维Ritz对第三模给约1.52111，不能假装已经收敛；6、8、12维结果趋于1.28943004550。最低模在6维后已接近射击值。更严格背景、1e-12射击容差和起点减半后，最低omega²变化约2.50e-11。

节点数和模间比较仍是有限数值检查，不承担第5节正能量证明。基频按无穷远坐标时间计，不混成表面固有时间频率。

## 7 移动表面恢复守恒

最低模归一化为xi(R)=R。局部Euler密度变化为 $\delta\rho=(\Delta p-\xi p')/s$，粒子数变化为 $\delta n=n\delta\rho/(\rho+p)$。固定半径积分并非整个星体的守恒量变化，还需边界：

$$
\delta N=\int_0^R4\pi r^2\sqrt B(\delta n+n\delta\lambda)dr
+4\pi R^2\sqrt{B(R)}n_s\xi(R),
$$

$$
\delta M=\int_0^R4\pi r^2\delta\rho\,dr+4\pi R^2\rho_s\xi(R).
$$

数值中，粒子数体积项约-9.64374627383、表面项+9.64374627363；质量体积项约-8.88369523376、表面项+8.88369523355。各总和约-2.04e-10、-2.17e-10，属于单位模归一化下的求积残差。

由质量约束，$\delta m(R)+4\pi R^2\rho_s\xi(R)=0$解析成立。遗漏表面项会错误地认为该径向模改变了一阶ADM质量。表面移动不被重新解释为额外薄壳物质。

## 8 重建钟率而不引入虚假外部信号

在无限远时间规范不变、外部一阶delta M=0的条件下，真空外部delta nu=0。内外背景nu'在p(R)=0时匹配，所以内部边界delta nu(R)=0，再向中心积分第2节的delta nu'。

最低模单位归一化下，r/R=.2、.5、.8的delta nu约为.39015989、.25018807、.07927975。外部为零。这是球对称约束与真空匹配，不是发射或探测到张量引力波；二阶模能量对背景的修正未计算。

重建后的五点差分检验给Euler残差约1.76e-13、质量约束残差约1.08e-9。两种微分形式、积分守恒和约束重建分别核对，但这些小数仍不是全域误差认证。

## 9 同一时变背景接回量子链

声明最低模相位cos(omega t)，表面相对振幅epsilon=1e-3。C10同一环上固定节点的度规一阶变化为

$$
\delta\ln\ell=\delta\nu,\qquad
\delta\ln b=\frac{B\,\delta\lambda\,(dr/dx)^2}{b_0^2},\qquad
\delta\ln v=\delta\nu-\delta\ln b.
$$

按C5同一中点规则及交错质量线性化，形成有限矩阵

$$
H_D(t)=H_0+\epsilon\cos(\omega t)H_1.
$$

这是本轮实际演化的定义；没有把系数公式中的高阶项偷偷当作已求解的非线性扰动。节点与边中点离静态表面的距离必须大于epsilon R，否则拒绝使用固定内外侧系数。环仍有固定导向和理想度规耦合，碰撞与探针对流体的反作用未计。

钟位置固定，其相位为

$$
\theta_r(T)=\omega_r^0\ell_{0,r}\left[T+
\epsilon\delta\nu_r\frac{\sin(\omega T)}{\omega}\right].
$$

不能只把T时刻的钟率乘T。外侧A钟delta nu=0，内侧C钟系数约0.3908959958，实际钟率相对调制约3.91e-4。

N=32、质量.35、T=2的完整数据–双钟记录相对静态C10任务的TV约0.00019715249。A钟X+仍约.3332295261，C钟Y+从.9274646071变为.9272674546。新H随时间变，数据能量由.3186934147变为.3186829435；差约-1.047122418e-5，与

$$
\int_0^T\langle\partial_tH_D\rangle dt
=-\epsilon\omega\int_0^T\sin(\omega t)\langle H_1\rangle dt
$$

一致。不能因数值演化酉就要求此能量保持，也不能据做功恒等式宣称“流体模＋探针”总能量已闭合：本轮未把这份交换反馈为模耗尽。

完整准备混合一致性、共同J实表示、零振幅恢复C10及相位积分均检查。时间演化采用DOP853，U的酉性残差约4.46e-13；更严格背景与模、探针积分的记录TV差约1.21e-11。后二者不包含遗漏的epsilon²引力误差，不能用来认证实际观测可分辨性。

## 10 记录和运输没有被省略

保留32位置比特、两钟读数，共34原始位，平滑输出用12个独立公平位和416位阈值表。后处理不丢弃原始记录，不重置源或钟。

C10的静态零线行时证明不自动适用于这次时变协议。本轮明确使用每跳1单位时隙的独立经典消息服务，272bit-hop，T=2读出、参数274收齐后发布。这个服务是接口输入，未构造时变度规中的实际通信装置或其零线误差界。

模的经典振幅和相位准备、流体微观资源、探针测试极限、同步、导向和测量能源仍未闭合；没有新采样或实验数据。

## 11 结论范围与下一轮12

本轮对同一个因果EOS球取得了明确的径向线性稳定性充分证明，而非只观察数值不爆炸。真正新增的是参数正能量界、正确的自束缚表面守恒项，以及与时变探针读数的连接。

但引力作用量仍是输入，非径向与非线性稳定、量子源约束和完整认知推导尚未完成。径向模也不是C9自由TT波的替代名称。

下一轮从该正能径向模出发，**先核定物理二阶作用量与正则幅度规范化，再构造模、数据和钟共同演化的有限量子H**。能量交换、准备混合及弱扰动有效域必须同账；不再把经典驱动做功当作已计入模耗尽。量子化单模仍不等于完整量子引力，必要性继续暂停。

## 12 复算

[代码](radial_fluid_modes.py)、[结果](radial_fluid_modes_results.json)、[实际验证记录](research_round_11_checks.json)。17项新检查、全路线180项通过，旧224篇编号笔记与705份源证据快照一致；公式、本地链接和编辑器诊断检查通过。参数正号证书用精确Fraction算术；ODE、本征频率、守恒残差与记录使用浮点，分别标注。

```powershell
python -B -X utf8 research_physics_construction/radial_fluid_modes.py --write-results
python -B -X utf8 -m unittest discover -s research_physics_construction -p "*.py"
python -B -X utf8 research_physics_construction/inheritance.py --verify
```

使用已有NumPy和SciPy，无新安装、Git、后台、现实实验或必要性研究。C0—C10材料不覆盖。