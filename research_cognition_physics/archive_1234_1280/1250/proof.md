# 1250：同一完整动态历史的有限 null 面积—boost 余额

日期：2026-10-15。两组作者科学检查；CT 新采用数为 0。最终范围见[验收记录](acceptance.json)及[独立审阅](independent_review.md)。

## 1. 命题、继承与净增

本轮在 [1247 的同一完整材料历史](../1247/proof.md)上，连接独立度规面积、完整物质 boost 流与 CT 已条件选出的 Einstein–Hilbert（EH）作用的 Wald 面积泛函。它是选后理论的有限一致性桥，不是第二次独立生成 Einstein 方程。

[316](../../archive_301_341/research_note_316.md)已有光学积分；[1243](../1243/proof.md)已有有限面积/热的不同权重及误差运输；[1248](../1248/proof.md)已有保持历史的时样来源见证。本轮新增是：将这些工具落实到 1247 原两支实际解，给出同对象的严格有限 null 余额，并判断一个具体记忆来源遗漏在何种预算下可辨。没有新选择物质、窗口或控制脉冲。

结果是：完整账本符合**扣除几何余项后的余额**；只漏记忆项会产生严格正的有限差。未经校正的 Clausius 容差在此窗口不能排除漏账，甚至漏账的原始残差更小。这一失败与校正后成功同时保留。

## 2. 固定历史与完整来源

取旧原三环面周期 $L=10$，两支仍为 1247 的程序 $alpha=0,1$。仅从冻结的 $t=2$ 端点严格区间继续到 $t_*=9/4$；不重跑旧 $[0,2]$。保持期 $w=f=0$，所以 $J,K,q^2+r^2$ 为常量，原完整钟质量为

$$
M_c=\frac{71}{100}+J+\frac{q^2+r^2+K^2}{2},\qquad
m_K=\frac1{100}+\frac{K^2}{2}.
\tag{1}
$$

这里 $m_K$ 包括记忆静质量与保留记录的动能；$J$ 已包含此前自治控制与写入代价。它不能重新求解或补偿遗漏。记 $u=1/a$，原变分材料给

$$
R_{\rm rod}=\frac3{20}\left(1+\frac p5\right)u^2,\quad
D_{\rm dust}=\left(M_c+\frac1{50}+\frac p{100}\right)u^3,
$$
$$
\rho=\frac{v^2+p^2+V^2+X^2}{2}+R_{\rm rod}+D_{\rm dust},\quad
P=\frac{v^2-p^2+V^2-X^2}{2}-\frac13R_{\rm rod}.
\tag{2}
$$

两标量 $v=\dot p,V=\dot X$、原三杆、原记录器、钟、程序、振子和记忆的来源均保留。演化沿用

$$
\dot u=-Hu,\quad \dot H=-\frac32(\rho+P),\quad
\dot p=v,\quad \dot v=-3Hv-p-\frac3{100}u^2-\frac1{100}u^3,
$$
$$
\dot X=V,\quad \dot V=-3HV-X,\quad \dot a=Ha.
\tag{3}
$$

由完整材料链式求导，$\dot\rho+3H(\rho+P)=0$，故 $d(H^2-\rho)/dt=0$。旧精确初约束因此延续；计算始终独立计算式(2)，不以 $H^2$ 定义物质密度。

## 3. 由独立度规固定局部平衡片

几何是原杆坐标中的 $ds^2=-dt^2+a^2(t)(dr^2+r^2d\Omega^2)$。事前固定

$$
t_* =\frac94,\quad w_t=\frac1{64},\quad R_* =H_*^{-1},\quad
\Omega_{\rm cap}=1,
\tag{4}
$$

取径向向内、未来指向的 null 束。它满足 $\dot r=-u$、$\dot R=HR-1$。固定 $k^t(t_*)=1$，令 $\lambda(t_*)=0$，则

$$
k^t=\frac{u}{u_*},\quad \frac{d\lambda}{dt}=\frac{u_*}{u},\quad
\theta=2k^t\left(H-\frac1R\right),\quad \sigma=0.
\tag{5}
$$

因此 $\theta_*=0$ 是由独立几何所选瞬时平衡片的精确性质，并非拟合热量。定义 $s=t_*-t$、$x=-\lambda$、$y=R/R_*$；则 $x(0)=0,y(0)=1$，

$$
x_s=u_*a,\qquad y_s=H_*-Hy,\qquad
\ell=x(w_t)>0,\qquad \mathcal F=A/A_*=y^2.
\tag{6}
$$

严格管给出 $y,u,H,a>0$。两支整片共动半径均小于 $1.430<5$，故无环面回绕；局部 cap 无焦散。该 cap 只是开放 null 片，不要求它是全局事件视界或封闭热循环。

由度规的 Ricci 张量与式(3)分别得到

$$
R_{kk}=-2\dot H(k^t)^2,\qquad
T_{kk}=(\rho+P)(k^t)^2,\qquad R_{kk}=3T_{kk}.
\tag{7}
$$

式(7)前两项的字典不同：一项是度规导数，一项是式(2)的变分来源；末项才是选后动力学的联系。程序的几何核由 $-2\dot H$ 构造，热核从完整 $\rho+P$ 单独构造。这里是同模型的几何计算，不宣称已独立实测曲率；它也不把热积分直接定义为面积或几何修正。

## 4. Wald 系数与两个有限权重

[1245](../1245/proof.md)的条件选择给 $16\pi G=6$、$8\pi G=3$。EH 的 Wald 局部面积泛函为

$$
S_W=s_W A,\qquad s_W=\frac1{4G\hbar}=\frac{2\pi}{3\hbar}.
\tag{8}
$$

这消去了本对象中另选面积系数的自由；并不证明 $S_W$ 等于裸区域 von Neumann 熵。沿用 [CT→Jacobson 审计](../_admission/ct_to_jacobson_bridge_after1249.md)中的条件参考温标 $T_U=\hbar\kappa_B/(2\pi)$。定义 boost 核 $\chi^a=-\kappa_B\lambda k^a$，归一热量

$$
\mathcal B=\frac{\delta Q/T_U}{s_W A_*}
=3\int_{-\ell}^0(-\lambda)\mathcal F T_{kk}\,d\lambda.
\tag{9}
$$

没有把这一定义提升为实际 KMS/Unruh 仪器资格；有限窗中的 $\chi$ 是声明的 boost 权，不声称为精确 Killing 场。

由 Raychaudhuri 方程 $d\theta/d\lambda=-\theta^2/2-R_{kk}$ 及 $\theta_*=0$，精确面积余额是

$$
\delta_A:=\frac{A_*-A_-}{A_*}
=\int_{-\ell}^0(\lambda+\ell)\mathcal F
       \left(R_{kk}-\frac{\theta^2}{2}\right)d\lambda.
\tag{10}
$$

面积权 $\lambda+\ell$ 与热权 $-\lambda$ 不可混同。定义纯几何修正

$$
E_{\rm weight}=\int_{-\ell}^0(2\lambda+\ell)\mathcal F R_{kk}\,d\lambda,
\quad E_{\rm exp}=-\int_{-\ell}^0(\lambda+\ell)\mathcal F\frac{\theta^2}{2}\,d\lambda,
\quad E_g=E_{\rm weight}+E_{\rm exp}.
\tag{11}
$$

故在完整来源下，$\delta_A-\mathcal B-E_g=0$。程序直接以 $1-y^2$ 计算 $\delta_A$，另积式(11)，从未以 $\delta_A-\mathcal B$ 定义 $E_g$。若实际末端只有近似 $\theta_*=0$，式(10)须加 $\ell\theta_*$；相应误差必须进入几何菜单预算，不能免费忽略。

## 5. 同几何记忆漏账

只在热来源账本删去 $m_Ku^3$、保留真实 $J$、几何与其他材料。保持期该项的压力为零，所以

$$
\mathcal B_- =\mathcal B-\mathcal D,
\qquad
\mathcal D=3\int_{-\ell}^0(-\lambda)\mathcal F m_Ku^3(k^t)^2\,d\lambda>0.
\tag{12}
$$

相同几何余额变成 $\delta_A-\mathcal B_- -E_g=\mathcal D$。这检验同一历史的完整来源身份；漏账不是另一套自洽 Einstein 解，也不排除改变材料字典的其他理论。

## 6. 严格复算与固定误差合同

[check.py](check.py)先校验冻结 [1247 结果](../1247/results.json)与 [1242 整数区间运行时](../1242/check.py)的 SHA256，只读导入后者的 $10^{-36}$ 外向舍入。先以 64 步续算 $[2,9/4]$，再以 16 步反向覆盖固定 $1/64$ null 窗。每步验证 Picard 包含，以八阶 Taylor 加整管九阶余项封闭；旧端点相关性被放宽，只会扩大包络。保持期向量场是解析多项式，不跨原 $C^2$ 脉冲接点求高阶导数。

反向累计 $G_0,G_1,Z_0,Z_1,\mathcal B,\mathcal D$，其中

$$
(G_0)_s=-2\dot H\,y^2k^t,\quad (G_1)_s=x(G_0)_s,
\quad (Z_0)_s=2k^t(Hy-H_*)^2,\quad (Z_1)_s=x(Z_0)_s.
\tag{13}
$$

因此面积光学式是 $\ell G_0-G_1-\ell Z_0+Z_1$，几何修正是 $\ell G_0-2G_1-\ell Z_0+Z_1$；$\mathcal B$ 另由式(2)积分。检查覆盖约束、仿射光学、正面积、两核与严格来源遗漏。

额外声明：归一面积、热量、几何修正及遗漏积分四个菜单，每项总误差不超过固定 $\eta=10^{-9}$。这是本轮的有限读口/模型匹配合同，包含采样、来源读数、端点定位、反作用与数值误差；实际仪器和该精度未从认知原则推导。完整校正余额误差至多 $3\eta$；由完整热量减遗漏计算的漏账余额误差至多 $4\eta$。两个允许区间分离的充分条件是 $\mathcal D>7\eta$。

## 7. 结果：原始容差不排漏账，校正余额排除

下表只作易读近似，最终严格端点见 [results.json](results.json)。

|归一量|$\alpha=0$|$\alpha=1$|
|---|---:|---:|
|独立面积 $\delta_A$|$2.0336325484\times10^{-5}$|$2.0325278752\times10^{-5}$|
|完整热核 $\mathcal B$|$2.0461841412\times10^{-5}$|$2.0450611466\times10^{-5}$|
|几何修正 $E_g$|$-1.2551592779\times10^{-7}$|$-1.2533271357\times10^{-7}$|
|其中膨胀修正|$-6.8843566193\times10^{-11}$|$-6.8768871667\times10^{-11}$|
|漏账的未经校正余额|$3.4413548936\times10^{-8}$|$3.4603321773\times10^{-8}$|

两支遗漏的严格包络分别为

$$
\mathcal D_0\in[1.59929476726864406412378282023,
1.59929476726864416861684526085]\,10^{-7},
$$
$$
\mathcal D_1\in[1.59936035345326980356464234168,
1.59936035345326998290716404774]\,10^{-7}.
\tag{14}
$$

完整校正残差的严格区间均包含零，绝对包络小于 $2.76\times10^{-19}$。两账按上述观测合同分离后的最小余量大于 $1.52929476726\times10^{-7}$。

相反，原始 Clausius 采用共同容差 $\Gamma=\max_\alpha|E_g|+2\eta\simeq1.2751593\times10^{-7}$ 时，两支漏账的全部误差区间仍在接受范围内。**不得以更小 raw 残差声称来源更完整；本窗口只能在显式保留几何修正后排除这项遗漏。** 本轮没有为 raw 分离重新选窗。

## 8. 结论与尚缺接口

在同一 CT 条件选后的有效 EH 历史中，Wald 系数、独立 null 面积、完整物质 boost 流与有限余项得到共同可复算余额。两组作者检查通过；独立审阅与最终验收分别保存。本轮没有新采用认知公理，没有交付实际 Unruh 温度、指定局域量子参考的模哈密顿量、几何熵的统计身份或这些四菜单的自治读出。

下一步若推进热学角色，应固定**同一历史、区域、观察者和完整参考来源**，再检验实际热/模能源到本轮热核的误差连接；不能把本轮有限余额重新计为独立 Einstein 推导，也不因缺少该参考而要求修复全部 UV 理论。参见 [NEXT](NEXT.md)。
