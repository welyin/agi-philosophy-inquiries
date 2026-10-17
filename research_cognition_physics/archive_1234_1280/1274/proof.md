# 1274：同一自治边界的几何充分性、结构分账与有限熵角色检验

2026-10-17。本轮作者核验1组。依据[准入](admission.md)，只复用[1237原证明](../1237/proof.md)及其完整有限载体；不增加材料、门或中途仪器。本文是旧固定过程的新轨迹推论，不是新的普遍熵定理。CT及新增几何充分性均未采用为公理。

## 1. 同一对象、参考和完整过程

原物理边界为
$$
\mathcal H_B=\mathcal H_4\oplus\mathcal H_6,\qquad
\dim\mathcal H_4=2,\quad\dim\mathcal H_6=5.
$$
它来自原144维运动学载体的完整SU(2)不变子空间。独立面积与边界能源为
$$
\widehat A=4aP_4+6aP_6,\qquad H_B=4eP_4+6eP_6,\qquad a,e>0.
\tag{1}
$$
面积来自原Casimir字典，不由熵或热流反定义；但该字典尚不是实际Lorentz钟尺或视界面积的证书。

指定原准备为$\sigma_4=I_2/2$，工作单元处于能源$2e$的激发态，相干旗标为720个置换标签的均匀纯态，程序为三站链的起点。全程保工作$W$、旗标$F$、程序$C$及原输入的任意伙伴$R$。只在推导实际边界约化态时取偏迹，不实际删除任何输出。

原等距$\mathcal W:\mathcal H_4\to\mathcal H_6$添加两腿singlet，令
$$
\tau=\mathcal W\sigma_4\mathcal W^\dagger,\qquad
\sigma_6=I_5/5.
$$
原$U_1$以工作能源激活六腿；原$U_2$按全部720个置换相干控制。其边界通道满足
$$
\operatorname{Tr}_{WF}\!\left[
U_2U_1(X\otimes|1,+\rangle\langle1,+|)U_1^\dagger U_2^\dagger
\right]=\frac{\operatorname{Tr}X}{5}I_5
\tag{2}
$$
对全部$X\in M_2$成立；完整输出等距保留未知输入及$R$。

实际生成元仍为
$$
H_{\rm tot}=I_C\otimes H_0+
g\bigl(|1\rangle\langle0|_C\otimes U_1+
|2\rangle\langle1|_C\otimes U_2+\mathrm{h.c.}\bigr),\qquad
H_0=H_B+H_W.
\tag{3}
$$
$4e>\sqrt2g$保证原完整正谱下界。路径装饰酉
$\mathcal D=\operatorname{diag}(I,U_1,U_2U_1)$给
$\mathcal D^\dagger H_{\rm tot}\mathcal D=I_C\otimes H_0+gA_3\otimes I$。
这些均直接继承1237，不把静态耦合当成免费的外部门。

## 2. 从完整路径等距到整个边界轨迹

令
$$
\vartheta=\sqrt2gt,\quad c=\cos\vartheta,\quad t_*=\frac{\pi}{\sqrt2g}.
$$
三站链从起点演化的振幅为
$$
a_0=\frac{1+c}{2},\qquad
a_1=-\frac{i\sin\vartheta}{\sqrt2},\qquad
a_2=\frac{c-1}{2},
$$
故
$$
p_0=\frac{(1+c)^2}{4},\qquad
p_1=\frac{1-c^2}{2},\qquad
p_2=\frac{(1-c)^2}{4},\qquad \sum_jp_j=1.
\tag{4}
$$
原输入码的自由总能恒为$6e$。完整等距是
$$
V(t)=e^{-6eit}\sum_{j=0}^2a_j|j\rangle_C\otimes V_j,\qquad
(V_0,V_1,V_2)=(V_{\rm in},U_1V_{\rm in},U_2U_1V_{\rm in}).
\tag{5}
$$
程序三站正交，所以偏迹$CWF$后没有路径交叉项。由式(2)，对全部输入矩阵单位及其线性组合，
$$
\Phi_t(X)=p_0X\ \oplus\
\left[p_1\mathcal WX\mathcal W^\dagger+
p_2\frac{\operatorname{Tr}X}{5}I_5\right].
\tag{6}
$$
这不是在路径中实际测量程序，也不是先选择程序分支。全输出仍按式(5)相干演化；式(6)张量恒等通道可处理原任意$R$。

对本轮指定的$\sigma_4$，
$$
\rho_B(t)=p_0\sigma_4\ \oplus\ (p_1\tau+p_2\sigma_6).
\tag{7}
$$
$\tau$在二维$\mathcal W$像空间的本征值为$1/2$，其余三维为0，故完整七维谱为
$$
\boxed{\operatorname{spec}\rho_B(t)=
\left\{\frac{p_0}{2}\ [\times2],\
\frac{p_1}{2}+\frac{p_2}{5}\ [\times2],\
\frac{p_2}{5}\ [\times3]\right\}.}
\tag{8}
$$
约定$0\log0=0$。端点为$(1/2,1/2,0,\ldots,0)$及$(0,0,1/5,\ldots,1/5)$，没有丢弃零谱。

### 2.1 限定的几何充分性正结果

原独立面积统计给
$$
\overline A(t)=6a-2ap_0,\qquad
p_0=\frac{6a-\overline A}{2a},\qquad c=2\sqrt{p_0}-1.
\tag{9}
$$
因此在**同一$\sigma_4$、同一纯旗标、工作、程序、父$H_{\rm tot}$的原轨迹**上，$\overline A$确定式(4)、完整$\rho_B$及$S_B=-\operatorname{Tr}\rho_B\log\rho_B$。由于$\widehat A$只有两个本征值，其期望也确定完整面积概率分布；并非一次面积测量的单个结果就确定熵。

这是无拟合的$S_B=\mathcal S(\overline A)$限定正结果，不只是初末两个点的比值。它不声称所有准备、所有边界码、所有普通物质资料或所有过程中几何都充分；工作和记录沿轨迹改变，故它也不是旧“固定全部普通资料$m$”纤维判据的全族证明。式(6)对未知输入成立，不把该全码量词偷偷升级为式(9)的全输入熵结论。

## 3. 独立参考计数的结构项与完整熵分解

本模型已经分别给出两扇区的实际条件参考$\sigma_4,\sigma_6$。由它们定义
$$
L=(\log2)P_4+(\log5)P_6.
\tag{10}
$$
这是原参考的扇区信息算符；不是先用面积定义一个熵。将它与原Casimir算符比较，才得到
$$
L=\log2\,I+\frac{\log(5/2)}{2a}(\widehat A-4aI),\qquad
\Delta\langle L\rangle=\frac{\log(5/2)}{2a}\Delta\overline A.
\tag{11}
$$
这里只有两个固定扇区；这种仿射关系不是一般面积律或普适系数。

令$p=p_0$并以实际面积扇区权重构造比较参考
$$
\omega_p=p\sigma_4\oplus(1-p)\sigma_6,\qquad
h_2(p)=-p\log p-(1-p)\log(1-p).
$$
这一步是声明的参考运输规则，不是新增Gibbs准备。由$\log\omega_p$在各扇区为标量，直接有
$$
\boxed{
S_B=h_2(p)+\langle L\rangle-D(\rho_B\Vert\omega_p).
}
\tag{12}
$$
$p=0,1$时按共同支持及连续极限理解。相对熵有限；实际态在比较参考的支持内。

式(12)把三个独立可算角色分开：面积扇区混合、条件参考结构计数、对该参考的非均匀偏离。普通记录与关联可以承担其中的微观来源，无须另设非物质熵库；也不得重复计入。它给用户“结构熵”的一个更精确候选$\langle L\rangle$，并说明为何不应未经检验把全部$S_B$代替它。

原数值给候选边际比例
$$
\frac{d\langle H_B\rangle}{d\langle L\rangle}
=\frac{2e}{\log(5/2)}.
\tag{13}
$$
此比例尚未认证为同一实际观察者的热温标。条件扇区参考本身因能源简并不能选择$\beta$；以下用完整父生成元另外确定一个参考身份。原$a,e,g$及$e/a$继续是输入。

### 3.1 同一完整相互作用的Gibbs参考与有限修正

不改变父Hamiltonian，对任意$\beta>0$定义其正常有限维参考$r_\beta=e^{-\beta H_{\rm tot}}/Z_{\rm tot}$。由原路径装饰，
$$
r_\beta=\mathcal D\bigl(r_{C,\beta}\otimes r_{H_0F,\beta}\bigr)\mathcal D^\dagger,
\quad
r_{C,\beta}=\frac{e^{-\beta gA_3}}{Z_C},\quad
r_{H_0F,\beta}=\frac{e^{-\beta H_0}}{Z_0}\otimes\frac{I_F}{720}.
\tag{13a}
$$
各$\mathcal D_j$与$H_0$对易，因而对$C$取偏迹后，每个对角项都给同一$r_{H_0F,\beta}$。所以
$$
\operatorname{Tr}_{CWF}r_\beta
=r_{B,\beta}
=\frac{e^{-\beta H_B}}{2e^{-4\beta e}+5e^{-6\beta e}}.
\tag{13b}
$$
这是完整相互作用参考的精确边际，不是弱耦合近似，也不是把裸参考冒充全参考；全$r_\beta$通常仍含程序与数据相关。

由原能源与原条件参考计数独立取
$$
\beta_*=\frac{\log(5/2)}{2e}>0,\qquad
r_{B,\beta_*}=\frac{P_4}{4}+\frac{P_6}{10},\qquad
-\log r_{B,\beta_*}=L+\log2\,I.
\tag{13c}
$$
故参考边界模算符与$L$只差中央常数。对本轮每个真实时刻，令
$D_B(t)=D(\rho_B(t)\Vert r_{B,\beta_*})$，则
$$
D_B=\langle L\rangle+\log2-S_B,\qquad
\boxed{\Delta S_B=\beta_*\Delta\langle H_B\rangle-\Delta D_B,\quad
\Delta\langle L\rangle=\beta_*\Delta\langle H_B\rangle.}
\tag{13d}
$$
有限修正$-\Delta D_B$由同一全参考和真实后态独立定义，不是根据失配任意补入。全$BE$实际态相对于同一全$r_{\beta_*}$的相对熵则因同$H_{\rm tot}$酉演化保持，局部修正不代表全账资源凭空产生。

因此本轮不仅删除一种错误识别，也给出一个正向共同身份：原结构计数、原几何算符、全父参考的边界模量及完整边界能源增量在这份有限模型中相符。该身份未证明参考的实际准备、原过程初态热化、实际观察者时钟与$\beta_*$归一、物理视界热流或GR。$\omega_p$只是式(12)的面积权重比较参考；不能把所有$p$都称为同一固定温度的Gibbs态。

## 4. 同一实际过程的有限区间：能量增加、边界熵净减少

比较同一自治演化的$t_m=t_*/2$与$t_f=t_*$，不在$t_m$插入操作。中点$c=0$时
$$
(p_0,p_1,p_2)=\left(\frac14,\frac12,\frac14\right),\qquad
\operatorname{spec}\rho_B(t_m)=
\left\{\frac18[\times2],\frac3{10}[\times2],\frac1{20}[\times3]\right\}.
\tag{14}
$$
终点$S_B(t_f)=\log5$。所以
$$
D_*:=S_B(t_m)-S_B(t_f)
=\frac{33\log2-5\log5-12\log3}{20}
=\frac1{20}\log\frac{2^{33}}{5^5\,3^{12}}> \frac1{15}.
\tag{15}
$$
严格界只用整数与正项级数：
$2^{33}>4\,5^5\,3^{12}$，以及
$\log2=2\operatorname{arctanh}(1/3)>2/3$。
数值显示仅作诊断，严格证据不依赖浮点对数。

同一区间有
$$
\Delta\overline A=\frac a2,\qquad
\Delta\langle H_B\rangle=\frac e2,\qquad
\Delta\langle H_W\rangle=-\frac e2,\qquad
\Delta S_B=-D_*<-\frac1{15}.
\tag{16}
$$
不是整个后半段逐点降熵：在$c=0$，
$$
\left.\frac{dS_B}{dc}\right|_0
=\frac1{10}\log\frac{1152}{3125}<0,
$$
而$c$随时间减小，故熵最初仍上升，之后才下降。式(16)只使用严格有限净差。

结构项则给$\Delta\langle L\rangle=\tfrac14\log(5/2)>0$。这与式(16)相容：式(12)中的混合与偏离项同时变化，不能丢掉后把总区域熵当作结构项。

### 4.1 只排除一种字面联合识别

$t_m\le t\le t_f$上$\langle H_B\rangle$单调增加。若额外声称：

1. 本过程全部$d\langle H_B\rangle$均为所选边界的热量；
2. 全部$S_B$就是所提结构热熵；
3. 具有正温$T(t)$，且没有其他已推导修正；

则$\int d\langle H_B\rangle/T(t)>0$，与式(16)的负熵变化冲突。此结论不需要拟合温度，但它**仅否决上述联合识别**。本模型实际供能来自原工作单元，程序和相干旗标生成相关；此前并未认证这些能源全是热流。用户的新猜想允许结构功、记录、关联和有限修正，未断言第1—3项，所以本轮不是用户猜想或完整认知纲领的反例。

## 5. 全资源、旁挂伙伴与数学稳定性

全自由能源始终为$6e$，相互作用$K=H_{\rm tot}-H_0$满足
$$
\langle K\rangle=0,\qquad\langle K^2\rangle=g^2.
\tag{17}
$$
零均值没有删除相互作用。全父Hamiltonian、工作、旗标、程序均未在中途重置或切断。至$2t_*$原自治演化返回原态；未领取永久末态或免费循环参考。

记$E=WFC$。对指定$\sigma_4$，全$BE$熵恒为$\log2$，故每时刻精确有
$$
S_B-\log2=I(B:E)-S_E.
\tag{18}
$$
若用二维$R$纯化初始$\sigma_4$，完整$BER$纯态给$S_E=S_{BR}$。核验器通过实际路径态的$14\times14$约化矩阵计算此谱，未构造巨大密度矩阵，也未丢失环境。原未知输入及$R$的全输出可恢复性继承式(5)；本轮的熵数值仍只针对指定参考准备。

为说明有限结论不是零误差尖点，另给**纯数学**稳健预算。若两端实际边界态分别与式(7)半迹距至多$\delta=10^{-3}$，七维熵连续界给每端
$$
|\delta S_B|\le\delta\log6+h_2(\delta)<\frac1{100}.
$$
这里$\log6<2,\log1000<7$及$-(1-\delta)\log(1-\delta)\le\delta$已足够。因此净熵下降仍严格大于
$$
\frac1{15}-\frac1{50}=\frac7{150}.
\tag{19}
$$
面积和边界能源的谱幅分别为$2a,2e$，两端差的误差至多$4a\delta,4e\delta$，故两项净增分别至少$62a/125,62e/125$。此为对已指定两态的稳定性条件，不是新增中途读口、时钟、停止门或仪器的误差认证；原1237实际权限不扩大。 两端误差只保负的净熵差与正的净能源差：稳健的零修正热身份删除可直接用于恒定正温。若温度随时间改变，还须另证实际全段能流非负或相应加权能流的正下界；两端预算本身不授予此结论。理想轨迹的全段能源单调性已在§4.1证明。

## 6. 查重、执行和停止线

- 1237 §2.4及原检查器已经核中途面积、工作与全能源；§3和信息核验只算初末、停止在$\tau$的熵及末旗标账。它没有给式(7)—(12)的整条轨迹或式(15)的中途净差。
- [热机制准入旧稿](../_admission/thermal_mechanism_admission_after1270/geometric_thermal_identity_candidate_after1270.md)已说明1237条件最大混合参考代入317不能选择面积系数。本轮不重开该选择。
- [1238后范围修正](../_admission/thermal_admission_independence_audit_after1238.md)允许普通关联成为候选几何熵来源；本轮依此保留结构项，不要求扣去全部物质之后另找熵库。
- 与1178、1270—1272的公共关系、实际钟尺及热部门不拼接。它们的源、时钟和热仪器权限不是本边界载体已有的权限。

[check.py](check.py)只读导入1237的实际载体构造，取得144维Gauss来源的等距与全部720置换；重新构造原$U_1,U_2$、三站路径和完整输出，先偏迹后比较式(6)—(8)，再核独立面积、能源、结构分解、环境账与严格有理余量。对原720个旗标块中的实际42维完整Hamiltonian分别指数并取偏迹，另核式(13b)—(13d)，不只把模身份重复代入自身。它不只把七谱公式重复代入自身。整条连续轨迹的证明在正文，有限采样仅作复算诊断。本轮1组；不把继承载体重复计作新组。默认只读比对[results.json](results.json)，仅首次显式--write排他创建。

本轮正结果是原指定参考轨迹的几何充分性及独立结构分账；删除结果仅限§4.1的字面联合识别。尚未得到同一实际钟归一的温标、真实因果热流、Lorentz边界面积、几何动力学或三项共同交付。完成本次检验后停止该载体扩建；不得靠另造温度、旗标或事后修正把失配拟合掉。

