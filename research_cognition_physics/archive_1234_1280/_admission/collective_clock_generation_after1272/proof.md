# 精确受约束初值、集体钟窗口及框架选择的范围

本项是成熟候选中的有限集体解析筛查，不是新的认知公理或完整量子引力证明。本文$T_c$专指Lee最终版的$\tilde t_c$，不是其式(6)—(10)中的另一组同名变量。$\alpha$指该文的$\tilde\alpha$。

## 1. 接受的候选输入

在[Lee最终版式(11)—(13)、(22)、(29)](https://arxiv.org/html/2212.14011v5#S3)的经典集体层，采用

$$
Q=q^Tq+p_1+p_2,\qquad
U=ss^T+\sum_{c=1}^2\left(4T_cp_cT_c+\frac1{16}p_c^{-1}\right),
$$

$$
\mathcal H=-U+\alpha UQU=0,\qquad
\mathcal G=sq+2\sum_cT_cp_c=0.
$$

$1/(16p_c)$来自原复变量移位，不能丢去。本文不作有限$N$物理态或测量概率断言。有限两站仅检验集体方程的代数子问题，未证明是原文大站数空间近似的有效截断。

原演化式为

$$
\dot T_c=-4T_cvT_c-\alpha UvU+\frac1{16}p_c^{-1}vp_c^{-1}-yT_c-T_cy^T,
$$

$$
\dot p_c=4p_cT_cv+4vT_cp_c+p_cy+y^Tp_c,
\quad \dot s=-2\alpha UvUq^T-ys,
\quad \dot q=2s^Tv+qy.
$$

这是额外候选规则，并未由四原则、五操作公理或用户的生成要求唯一选出。论文的GR代数选择条件、局域态、框架及长波近似分别保留；没有“少纠缠钟单独推出全部GR”的主张。

## 2. 精确初值族

令$X$为两站交换矩阵，$X^2=I$。取

$$
q=p_1=p_2=I,\quad T_1=rX,\quad
T_2=bI-\frac r2X,\quad s=-2bI-rX,
$$

$$
b^2+\frac34r^2=1,\qquad \alpha=\frac8{195}.
$$

直接展开，所有$X$项抵消：

$$
U=\left(8b^2+6r^2+\frac18\right)I=\frac{65}8I=:uI,
\quad Q=3I,\quad \alpha u=\frac13.
$$

于是两项约束精确为零。本文取

$$
b=\frac{9997}{10003},\qquad r=\frac{400}{10003},
$$

来自有理参数$b=(1-3a^2)/(1+3a^2)$、$r=4a/(1+3a^2)$、$a=1/100$。此处$a$是初值参数，不是物理截止尺度。

## 3. 初片：框架旋转与钟矩阵切向量不同

任取对角lapse $v=\theta$，原线性shift类别在本初片为

$$
Y_\kappa=-2\kappa\theta s+(1-\kappa)(-4s\theta+2\theta s).
$$

定义$R=4T_1\theta+Y_\kappa$。逐项相加得

$$
R=4b\theta+(8-4\kappa)rX\theta+(4\kappa-2)r\theta X,
$$

$$
R-R^T=2(5-4\kappa)r[X,\theta],
\qquad
\dot p_1=R+R^T=8b\theta+6r\{X,\theta\}.
$$

这两条是该**初片**的精确式，无小$r$遗漏项。对所有非交换对角lapse要求前式为零，在此类别选出$\kappa=5/4$；取$\theta=I$则不能选择它。后式对$\kappa$完全独立，例如$\theta=\operatorname{diag}(1,2)$时非对角切向量为$18r\ne0$。

因此“无反对称框架旋转”与“$p_1$不生非对角项”不是同一条件。这不反驳作者明确限定的框架条件，也没有证明不同$\kappa$存在可辨物理概率。在$p_1=I$处反对称生成元本来就不改变它的一阶值。这里非对角项是集体关联诊断，未单独认证量子纠缠熵。

在此初片，站点对角钟读数对对角lapse的速率Jacobian为$8bI_2$，行列式$(8b)^2>49$。这是初片的横截性，不是全窗口多局部钟规范唯一性的证明。

## 4. 同一原方程中的非零钟窗

固定$v=I$，同一类别给$y=-2s$，与$\kappa$无关。所有矩阵在$I,X$代数内，且$q=I$保持。对$X$特征值$\sigma=\pm1$分别处理，令

$$
s_{0,\sigma}=-2b-\sigma r,\quad d_\sigma=s_{0,\sigma}^2+\alpha u^2,
\quad D_\sigma=1-4s_{0,\sigma}\tau+4d_\sigma\tau^2,
$$

$$
U_\sigma=\frac u{D_\sigma},\qquad
s_\sigma=\frac{s_{0,\sigma}-2d_\sigma\tau}{D_\sigma}
=-\frac{D_\sigma'}{4D_\sigma}.
$$

分母满足$D_\sigma=(1-2s_{0,\sigma}\tau)^2+4\alpha u^2\tau^2>0$。令

$$
B_{1,\sigma}=8b+12\sigma r,\qquad B_{2,\sigma}=16b,
$$

$$
p_{c,\sigma}=1+B_{c,\sigma}\tau+\frac{B_{c,\sigma}^2+1}{4}\tau^2,
\qquad T_{c,\sigma}=\frac{s_\sigma}{2}+\frac{p_{c,\sigma}'}{8p_{c,\sigma}}.
$$

这些表达式精确回到第2节初值。定义$z_c=T_c-s/2$，原方程在每分支化为

$$
\dot U=4sU,\quad \dot s=-2\alpha U^2+2s^2,
\quad \dot p_c=8p_cz_c,\quad
\dot z_c=-4z_c^2+\frac1{16p_c^2}.
$$

上述二次多项式有恒等式$2p_cp_c''-(p_c')^2=1$，故满足后两式；$D,s$的导数直接满足前两式。再验证**定义与约束**而非仅假定它们传播：

$$
Q=1+p_1+p_2=3D,\qquad \alpha UQ=1,
$$

$$
\mathcal G=s+2\sum_cT_cp_c=sQ+\frac{Q'}4=0,
$$

$$
s^2+\sum_c\left(4p_cT_c^2+\frac1{16p_c}\right)
=s^2Q+\frac{sQ'}2+\frac{Q''}8
=\frac{Q''}8-s^2Q=\frac uD.
$$

因而原$U$定义、两项约束和全部原ODE共同成立。矩阵解由两个特征分支固定变换回原站点基底即可。

### 有限范围证书

在$0\le\tau\le1/100$，$8b-12r>7$、$16b>15$，所以$p_c\ge1$且各分支$p_c'>7$。原站点的两个$p_1$对角读数相等，记为

$$
C_1(\tau)=1+8b\tau+\left(\frac{65}4+24r^2\right)\tau^2,
\qquad C_1'(\tau)\ge8b>7.
$$

它可在此单参数轨道上用作集体计时变量。$p_2(\tau)I$另给严格单调且站点非对角项为零的集体钟。两者是可比较的不同读数函数，并非要求所有工作钟同速。

与此同时

$$
(p_1)_{12}=12r\tau(1+4b\tau),\qquad
0\le\frac{|(p_1)_{12}|}{\sqrt{(p_1)_{11}(p_1)_{22}}}
\le\frac{78}{15625}<\frac1{200}.
$$

关联有严格有限上界，却在任何正$\tau$上非零。不能写“所有站点钟始终独立”。上述有限窗只是集体方程的精确性质，尚不含被忽略量子涨落或高导数误差界。

## 5. 本次没有获得的身份

本解不含三维空间的长波结构，不将$U$、$p$或任意其函数自行定义为已经测出的度规。它也没有证明独立边界读口、完整物质热流、Unruh温度、面积熵或Jacobson交换关系。现行TR共同热参考继续继承，但未自动在这一新候选中获得其特定几何身份。

同一模型原文在额外框架与近似下给出的传播结果，是后继可核的接口；完整正常量子概率及跨观察者共同实现仍须分别证明。先检验一个同对象有限几何预测，再判断是否投入更多细节。不得因本筛查通过便宣布强生成猜想或广义相对论完成，也不得因其尚缺这些接口便宣称全部认知公理不足。

## 6. 可复算性

[check.py](check.py)使用现有Python、标准库精确分数及NumPy对象矩阵；默认只读比对[results.json](results.json)。18项多项式系数恒等式检查覆盖全过程的显示代数；12项矩阵代入检查初片lapse/框架式；22个有理数分支时点的全部原ODE与约束替代另作实现诊断。解析全窗结论来自上述证明，不以有限采样代替。
