# 1225后准入：由原热参考推出规模一致的边界资源与面积上界

2026-10-15。短准入；不占新编号，不替换1225的面积读口工作，不修改冻结文件。以下是下一轮可证明的候选链及核查结果，尚未交付其逐项常数与正式验收。

## 1. 结论：值得进入，但条件应限定准确

**可以优先尝试由原完整Gauss Gibbs态本身推出单点/单边资源界，不另加每条边的热能预算。** 有限范围、有界度数、原束缚势和有限CAR次数允许把旧固定图估计改成节点数一致的估计；平移不变随后把总热能密度界转换成边界附近的局部能源界。这能补[1221](../research_note_1221.md)明确保留的“常数尚未证明随图规模一致”缺口，并接正在完成的1225独立面积字典。

声明范围必须是：固定格距/单位、固定有限代数与局域耦合；周期立方图族且局域作用范围、每点项数有共同上界；均匀正几何在事前紧集内，$0<w_-\le w\le w_+$、$g_-I\le G\le g_+I$、$\det G=1$；固定$\beta>0$或紧正温度区间；原**完整有限体积Gauss Gibbs态**，未额外挑一个排除试态的全局扇区。所有来源仍使用原H，不删除相互作用或费米质量。这些是所选旧物理候选的范围，不采用为新认知公理。

## 2. 去重结果

- [598](../../archive_585_628/research_note_598.md)有原完整CAR质量的束缚势控制，但其全图写法允许常数随N增长。[603](../../archive_585_628/research_note_603.md)给所有正β的完整Gauss热迹和能源矩；明确是固定图，没有规模一致密度结论。
- [637](../../archive_629_652/research_note_637.md)已给原单点H⁵和商群Casimir的热迹因子、区域熵有限性。其区域熵上界仍含区域体积及全部能源，不是本候选边界面积界。
- [639](../../archive_629_652/research_note_639.md)的体积加权读口具有节点数无关注能条件；它不是完整热态的单点束缚势/Casimir预算。
- [529](../../archive_467_530/research_note_529.md)已有另一个统一凸关系模型的量子BL协方差界；不能把那份平直构形和正路径测度移植成当前非线性Gauss/CAR定理。
- [1221证明§3](../1221/proof.md)给双向形式控制和保匹配比较态，但明确未证随图规模/格距一致。其双相对熵身份、中心/配对分账直接复用，不重新计新定理。
- 635—645、707—710、837—843的区域/参考/来源共同身份继续按[前份筛查](after1224_joint_geometry_thermal_screen.md)复用。本次没有找到其中对**原完整无界Gauss材料、原匹配比较态、热参考和真实面积字典**给出以下规模一致论证的已验收版本。

全库检索范围包括编号报告中的“规模一致热”“局部热能”“平移热态”“自由能密度”“边界资源面积”等，并回读上述原文件。该去重结论不是对全部外部文献的原创性宣称。

## 3. 第一个关键：逐节点控制，不能沿用可能超线性的全图Young常数

令N为站点数，每点原CAR次数$f=32n_g$固定。原完整质量在单点满足

$$
\|B_v(\phi_v)\|\le aU_v^{1/4}+b
\le\tfrac12w_-U_v+c_{\rm mass}.
$$

系数来自598原$F^{-2}\le AU+B_0$及固定有限矩阵，不随N变。先逐节点Young、再求和；**不能先把整图改写成$C+cW^{1/4}$，再把N依赖藏进常数。** 有限范围有界跳跃至多产生$c_{\rm hop}N$；正电动能、梯度、磁势保留。因此应逐项核得

$$
L_N:=H_N+c_*N\ \ge\ K_N
=c_T\Bigl(\sum_vL_v+\sum_eL_e\Bigr)+u\sum_vU_v\ \ge0,
\tag{1}
$$

$c_*,c_T,u>0$对全图族和紧几何域一致。$L_v=-\Delta_{H^5,v}$，$L_e$为原商群参考Casimir。非对角电动能由每个固定三边星的最小特征值控制；不能丢掉电动能后对无限维链路空间取迹。$c_*N$仅作分析移位，不从物理源或真空账中扣除。

## 4. 同一Gauss热迹夹界推出每点资源

603/637的局域热迹直接给

$$
Z_{L_N}(s)\le
\left[2^f z_\phi(s)\,z_{\mathsf G}(sc_T)^3\right]^N,
\quad z_\phi(s)=\operatorname{Tr}e^{-s(c_TL_v+uU_v)}<\infty.
\tag{2}
$$

这里用min–max本征值比较及Gauss子空间迹不超过全空间迹；不用错误的算符指数单调性。$z_{\mathsf G}$保原$\mathbb Z_6$商，不更换群。定义$a(s)=\max\{0,\log[2^fz_\phi(s)z_{\mathsf G}(sc_T)^3]\}$。

下界可由显式Gauss试向量取得：各节点同一规范不变径向紧支光滑包、各链路Haar常函数、全CAR空真空。它不是H的不变空间或基态，但有

$$
\langle L_N\rangle_{\rm trial}\le e_*N.
$$

有限支撑与固定范围控制各梯度星/磁项；全部质量、配对和跳跃按真实期望或局域范数入账。试态不受演化保持不妨碍此变分下界。谱测度上的Jensen给$\log Z_{L_N}(\beta)\ge-\beta e_*N$。

由$F_N(\beta)=\log Z_{L_N}(\beta)$凸，

$$
\langle K_N\rangle_\beta\le\langle L_N\rangle_\beta
\le\frac2\beta\left[F_N(\beta/2)-F_N(\beta)\right]
\le N e_\beta,
\qquad e_\beta=2e_*+\frac{2a(\beta/2)}\beta.
\tag{3}
$$

这使用同一完整Gauss热态；无需算全部热谱。平移在原CAR上由合法第二量子化置换实现，保持H和Gauss，所以其有限体积Gibbs态严格平移不变。把$K_N=\sum_vk_v$，每点含其三条出边，得

$$
\langle U_v\rangle_\beta\le e_\beta/u,
\qquad \langle L_{(v,i)}\rangle_\beta\le e_\beta/c_T.
\tag{4}
$$

无需旋转各向同性，不同方向可有不同实际均值；式(4)逐方向由非负和控制。任意非平衡有限总能态不能继承此单点结论；不能省去平移不变参考的范围。

## 5. 从同一边界能源到独立面积

固定区域A为原格胞并集，用1221同一个切分J及保匹配比较态$\tau[\sigma_\beta]$。它保持每点$U_v$及区域内比较能源；原切边Casimir也由同一标签概率保持。梯度簇用

$$
|\langle\operatorname{Log}_v(w),\operatorname{Log}_v(z)\rangle|
\le\tfrac12(d_{vw}^2+d_{vz}^2),\qquad r_v^2\le c_r(1+U_v)
$$

控制；磁簇和原CAR跳跃逐项有界。若某合法边界分账保留电交叉形式，亦由邻近Casimir的二次型Cauchy界控制，不得暗删。

每个跨区作用簇的支持在固定半径内且连通，至少经过一条切边；每条切边承载的簇数、各簇站点/链路数有统一上界。于是可证明存在不随N或区域A增长的$C_\partial(\beta)$，使

$$
0\le\Delta E_\partial
=\operatorname{Tr}[(\tau-\sigma_\beta)H_\partial]
\le C_\partial(\beta)\,b_A,
\tag{5}
$$

$b_A$为切边数。右侧使用式(4)和比较态的同一局域边缘，不需另加每边预算。真正的新增工作是核对各原作用簇及其N一致常数，不是再写一次三角不等式。

1225的独立字典如按其范围成立，给同一均匀格胞Gram矩阵$h=w^{2/3}G$的面面积$a_i=\sqrt{\operatorname{cof}(h)_{ii}}$。原格胞并集的边界面积

$$
A_\partial=\sum_{e\in\partial A}a_{i(e)},\quad
a_i\ge a_-:=w_-^{2/3}/\sqrt{g_+}>0,
\quad b_A\le A_\partial/a_-.
$$

与1221精确身份合并，候选主命题为

$$
\boxed{D(\sigma_\beta\Vert\tau)+D(\tau\Vert\sigma_\beta)
=\beta\Delta E_\partial
\le\frac{\beta C_\partial(\beta)}{a_-}\,A_\partial.}
\tag{6}
$$

面积来自独立几何记录/字典，不能由熵倒定义。若面积实际估计满足$|\widehat A-A_\partial|\le\delta_A$，安全概率/仪器范围内的上界使用$\widehat A+\delta_A$，不得把估计值无误差代入。

## 6. 可同轮保留的完整互信息扩展

1221给扩展区域互信息

$$
I(A:B)=H(p_{\rm joint})+2\sum_{e\in\partial A}\mathbb E\log d_{\lambda_e}
+D(\sigma_\beta\Vert\tau).
$$

不能免费删除前两项。对每条原商群切边，以
$d_\lambda^2e^{-tC_\lambda}/z_{\mathsf G}(t)$为比较概率，KL非负给

$$
H(p_e)+2\mathbb E\log d_{\lambda_e}
\le t\mathbb E C_{\lambda_e}+\log z_{\mathsf G}(t),\qquad t>0.
$$

再用$H(p_{\rm joint})\le\sum_eH(p_e)$及式(4)，得到同样N一致的完整互信息上界

$$
I(A:B)\le
\frac{\beta C_\partial(\beta)+t e_\beta/c_T+\log z_{\mathsf G}(t)}{a_-}\,A_\partial.
\tag{7}
$$

有限边界熵与无穷表示和的合法性直接复用637；若正式证明使用截断，也须共同保留原商群及其热迹极限。不要求标签独立，更不将$p_{\rm joint}$换成直积。

## 7. 成功标准与停止线

成功应交付一份对原局域项逐类核验的N一致常数链、同一真实热态的局域资源推论，以及(6)或连同(7)的有限区域面积上界；可复算局部常数和图簇计数即可，不需要全热谱数值，也不应另造少量qubit替代原Gauss系统。

这是**有限温度、固定格距、明确参考态族的关联/互信息面积上界**。不声称普通热区域熵$S(A)$满足同一面积律，不声称$S=sA$、$\delta S=s\delta A$、普适UV系数、Newton常数、热力学极限唯一态、真实热化或Einstein动力学。对$\beta\to\infty$、$w\to0$、节点非均匀、非平衡高能集中态，常数可失效或发散，应按其范围另审，不能升级成整个研究纲领失败。

如果逐点质量/动能常数、平移对称、切边能源保持或跨簇计数任何一项不成立，应记录具体失败项；不得用另加每边预算掩盖本轮欲消去的输入。1225实际面积读口未验收前，式(6)—(7)只保为带几何字典条件的候选，不抢先宣布TG3/TG4完成。
