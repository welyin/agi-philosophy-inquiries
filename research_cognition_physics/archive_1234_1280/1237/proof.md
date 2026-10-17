# 1237证明：同一有限载体的指定参考扩展与全输出身份

2026-10-15。与[准入](admission.md)、[复算](check.py)、[结果](results.json)一起使用。解析命题在公开输入下成立；NumPy只核验这一个固定有限实现，不替代全码证明。

## 0. 工具来源、认知动机及适用域

认知动机是扩大共同参考时，既不把旧未知信息销毁，也不把新参考当作免费现成资源。新增候选角色只要求已接入同类面片的**指定参考**在实际许可交换下不变；它不是所有认知态都须均匀，也不是被动坐标协变。

几何字典复用[Bianchi–Donà–Speziale](https://arxiv.org/html/1009.3402)的交织子／量子多面体与Casimir面积工具；singlet计数及指定最大混合参考熵直接复用[Livine–Terno视界模型](https://arxiv.org/html/gr-qc/0508085)和[零角动量混合态](https://arxiv.org/html/quant-ph/0502043)。这里没有证明任意singlet是平滑多面体、Lorentz视界或同一原钟的距离记录。

[1181](../../archive_1161_1200/research_note_1181.md)的原EPRL全码配对与[1182](../../archive_1161_1200/research_note_1182.md)的纯标签旋量障碍均不改动；本例不是在原纯EPRL码内添加免费旋量。本文也不重证[1042](../../archive_1009_1043/research_note_1042.md)的中心熵自由、[317](../../archive_301_341/research_note_317.md)的平衡选择恒等式、[1226](../../archive_1219_1233/research_note_1226.md)的面积上界或[1230](../../archive_1219_1233/research_note_1230.md)的面片响应。净新增是下面三项接口在同一有限过程中的共同实现和失败边界。

## 1. 独立几何与指定参考的选择

### 1.1 同一载体、两个实际支持

取

$$
\mathcal H_{\rm kin}=V_{1/2}^{\otimes4}\otimes(V_0\oplus V_{1/2})^{\otimes2},\qquad
\mathcal H_B=\operatorname{Inv}_{SU(2)}\mathcal H_{\rm kin}
=\mathcal H_4\oplus\mathcal H_6.
$$

运动学维数为144，Gauss不变空间维数为$2+5$；下文“物理空间”仅指这个约束空间，不表示已认证实际时空身份。后两槽位同为$V_0$时得到四半自旋singlet；同为$V_{1/2}$时得到六半自旋singlet；只激活一腿时，中心元$-I\in SU(2)$在五半自旋空间上作用为$-1$，故没有不变态。维数2和5是既有Catalan计数。

独立声明

$$
\widehat A=a_*\sum_{i=1}^6\sqrt{J_i^2}=4aP_4+6aP_6,
\qquad a=\frac{\sqrt3}{2}a_*>0,
\qquad Q_{56}=J_5\cdot J_6.
$$

两条额外槽位一直在原载体中，激活改变实际Casimir支持，不是比较时重新划分系统。$a_*$仍自由，不代入含$G$的Planck面积。面积谱排序、SU(2)字典及其可观测含义均为新有效模型输入。

令$|s\rangle_{56}$为归一化二半自旋singlet。等距

$$
W:\mathcal H_4\longrightarrow\mathcal H_6,\qquad W|\psi\rangle=|\psi\rangle\otimes|s\rangle_{56}
$$

在同一物理空间两扇区之间定义，$W^\dagger W=I_2$。

### 1.2 真正的选择力及删弱

六半自旋singlet的$S_6$作用是不可约$(3,3)$表示。故如果**指定参考本身**在全部许可主动交换下不变，Schur引理给

$$
\sigma_6=I_5/5.
$$

SU(2)在这个重数空间上作用平凡，不能选出此态；仅说重标号时描述协变也不够。实际执行全部交换及指定参考不变性是组合／参考层新增角色输入，尚未采用为现行公理。

取$\sigma_4=I_2/2$、$\tau=W\sigma_4W^\dagger$。在原成员／新成员分别交换的$S_4\times S_2$下，$\tau$仍不变。其面积和下面定义的能源均与$\sigma_6$相同，但秩2、熵$\log2$，且

$$
\frac12\|\tau-\sigma_6\|_1=\frac35.
$$

一个原形状读口就足以分辨，无须把全层析塞入几何资料。$\tau$的末两腿是singlet，故$\langle Q_{56}\rangle_\tau=-3/4$。在$\sigma_6$中，全部15对相等，而总自旋为零给

$$
\sum_{i<j}J_i\cdot J_j=-\frac12\sum_iJ_i^2=-\frac94,
\quad\Longrightarrow\quad
\langle Q_{56}\rangle_{\sigma_6}=-\frac3{20}.
$$

因此$P^{(0)}_{56}=I/4-Q_{56}$的结果概率由1变为$2/5$。全交换角色删除了同面积、同能源下的一份参考自由；它还可能抹去指定经典形状，不能因此称为平滑视界。

## 2. 一份自治保能总过程

### 2.1 所有准备和能源输入

取$e,g>0$，声明

$$
H_B=4eP_4+6eP_6,\qquad H_W=2e|1\rangle\langle1|_W,\qquad H_0=H_B+H_W,
\qquad 4e>\sqrt2g.
$$

形状简并和所取能源差均是模型输入；尤其$e/a$不是研究选择出的响应系数。对$H_B$，条件扇区Gibbs态对任意$\beta$都为$\sigma_N$，其全7维Gibbs态的两个权重则为$2e^{-4\beta e}$与$5e^{-6\beta e}$的归一化值。它们并不是交互开启时完整$H_{\rm total}$的Gibbs参考，也未取得TR完整材料热参考资格。跨扇区权重、温标与Clausius归一未选出。

另有720维零自由能源旗标$F$和3维零自由能源程序$C$。初态登记为

$$
\rho_{\mathcal H_4R}\otimes|1\rangle\langle1|_W
\otimes|+\rangle\langle+|_F\otimes|0\rangle\langle0|_C,
\qquad |+\rangle_F=\frac1{\sqrt{720}}\sum_{\pi\in S_6}|\pi\rangle_F.
$$

$\rho_{\mathcal H_4R}$任意，$R$任意有限维。工作激发、纯相干旗标、纯程序、固定耦合及初始载体都是公开准备资源，不声称由TR免费产生。指定参考转换另取$\rho_{\mathcal H_4}=\sigma_4$，没有把一般未知输入直接假定为参考。

### 2.2 两个保能算符

在$\mathcal H_B\otimes\mathcal H_W$中，令$S$为嵌入$|\psi\rangle\mapsto|\psi\rangle\otimes|1\rangle_W$，$T$为$|\psi\rangle\mapsto W|\psi\rangle\otimes|0\rangle_W$。二者正交等距且都在能源$6e$支持。定义

$$
U_1=I-SS^\dagger-TT^\dagger+TS^\dagger+ST^\dagger.
$$

于是$U_1=U_1^\dagger$、$U_1^2=I$、$U_1S=T$、$[U_1,H_0]=0$。它对全部旧码确定性激活末两腿，工作下降$2e$、边界能源上升$2e$、面积上升$2a$。

令$U_\pi$为六腿物理singlet中的交换表示，在$\mathcal H_4$块取恒等，并定义

$$
U_2=\sum_{\pi\in S_6}(I_{\mathcal H_4}\oplus U_\pi)
\otimes I_W\otimes|\pi\rangle\langle\pi|_F.
$$

这在同一$2+5$空间中是全酉且保能。只有六腿已激活时才执行其内部置换，不要求把真空槽位换入不允许的原四腿表示。全720旗标保留，不测量后选或免费取迹。对初始旗标$|+\rangle$，边界约化通道是

$$
\mathcal E(\rho)=\frac1{720}\sum_\pi U_\pi W\rho W^\dagger U_\pi^\dagger
=\frac{I_5}{5}\operatorname{Tr}\rho.
$$

最后一式是不可约表示平均；只是约化描述。完整演化保留全部未知输入，逆$U_2$与$U_1$精确恢复。恢复时须访问全部旗标和工作，不能声称边界单独保住了原码。

### 2.3 静态三站程序，消去免费定时门

以上两门不按外部脉冲执行。声明一份始终不变的Hamiltonian

$$
H_{\rm total}=I_C\otimes H_0+K,\qquad
K=g\bigl(|1\rangle\langle0|_C\otimes U_1
+|2\rangle\langle1|_C\otimes U_2+\mathrm{h.c.}\bigr).
$$

这是有限有效内部演化规则。静态耦合本身是输入，不是已从现行公理派生的设备许可。置换受控和激活不能脱离这一共同算符分别宣称免费。

设$D=\operatorname{diag}(I,U_1,U_2U_1)$，$A_3$为三站链邻接矩阵。则

$$
D^\dagger H_{\rm total}D=I_C\otimes H_0+gA_3\otimes I,
\qquad \operatorname{spec}(A_3)=\{-\sqrt2,0,\sqrt2\}.
$$

这个恒等式对整个物理数据／工作／旗标空间成立，不限初始参考。故总Hamiltonian下界为$4e-\sqrt2g>0$，无需删除相互作用来得到正能源。起点$|0\rangle_C$在

$$
t_*=\frac{\pi}{\sqrt2g}
$$

精确到达$|2\rangle_C$，并对任意旧码执行$-e^{-6eit_*}U_2U_1$。全码恒等式张量$I_R$即保任意参考；不是只核验某一参考态保真度。有限有效时间及读出时刻仍是输入，未接入实际Lorentz钟或观测误差装置。

合同只针对有限终端。程序没有在$t_*$停机，至$2t_*$同一静态Hamiltonian的前向演化会把旧码、工作、旗标和程序全部送回初态，差一个共同相位。因此既有真实前向恢复，也明确没有永久稳定的末态或免费循环参考。若中途加入实际读写、固定记录或停机动作，其源账须另加，不能由本例授予。

有限时刻误差不要求先建新装置。初始码的自由总能均为$6e$，去共同相位后，若终端误差为$\delta t$，算符范数界给全部输入及$R$的等距误差

$$
\|V(t_*+\delta t)-V(t_*)\|\le\sqrt2g|\delta t|.
$$

同一输入的全输出半迹距离至多该界（再与1取小）；实际读时精度尚未实现。它给有限误差合同，而不自动赋予额外时钟权限。

### 2.4 动态资源账

记$c=\cos(\sqrt2gt)$。起点概率$p_0(t)=[(1+c)/2]^2$；其余两站均已激活六腿。因此任意初始旧码下

$$
\langle A\rangle_t=4a+2a(1-p_0),\qquad
\langle H_B\rangle_t=4e+2e(1-p_0),\qquad
\langle H_W\rangle_t=2ep_0.
$$

全自由能源恒为$6e$。同时$\langle K\rangle_t=0$、$\langle K^2\rangle_t=g^2$；相互作用的零均值并不意味着它不存在或没有资源资格。完整$H_{\rm total}$保能，工作、程序、相干旗标、相互作用以及全部输出关联都在账内。$\Delta A$与工作变化的比例由输入$e/a$控制，未选择Newton响应。

## 3. 实际指定参考熵、全R身份与删除失败

### 3.1 不以容量代替实际态

初始取指定实际密度矩阵$\sigma_4=I_2/2$。只有激活时为$\tau$，实际熵仍$\log2$。完成全受控交换并约化边界后，实际态为$\sigma_6=I_5/5$，故

$$
\Delta A=2a,\qquad \Delta S_B=\log(5/2)>0.
$$

熵数字来自这两份实际态，不是把可用码维数直接叫已发生熵。两端面积本来由Casimir定义；但是这并不反过来证明熵变化属于几何热余额。

### 3.2 全旗标谱与未知参考的去向

记旧码$Q=\mathcal H_4$，设$V:Q\to B\otimes F$为末态等距，工作/程序已到固定末态。约化边界通道为$\mathcal E(X)=I_5\operatorname{Tr}X/5$。取旧码的最大纠缠纯化$|\Phi_2\rangle_{QR}$，则

$$
\rho_{BR}=\frac{I_5}{5}\otimes\frac{I_2}{2}=\frac{I_{10}}{10}.
$$

纯态$BFR$的旗标与$BR$非零谱相同，故旗标有十个$1/10$；其余710个特征值为零。因此，在仅$BF$而不计初始纯化$R$的账中

$$
S(B)=\log5,\quad S(F)=\log10,\quad S(BF)=\log2,
\quad I(B:F)=2\log5,
\quad\Delta S_B=I(B:F)-S(F).
$$

输入未知信息并未消失；它在全输出中，可由逆演化恢复，但局部边界与旧$R$完全脱钩。旗标既不是无变化催化剂，也不能在仍要求恢复时免费重置。若最初加入$R$，全$BFR$纯度和纠缠身份一并保留。本组明确采用纯相干旗标，不能与经典随机旗标的熵公式混用。

此增熵全部有普通记录／关联来源。本轮尚未找到扣清这部分后独立的几何热余额，故不能写为$\Delta E/T$、$s\Delta A$或声称满足317的平衡条件。所选扇区能源简并也无法由内扇区参考指定$T$。720维纯旗标及$\log10$旗标熵属于本声明实现，并非所有扩容协议的最低资源成本；不排除混合旗标或其它有限实现，也不另开压缩支线。

### 3.3 删除旗标恢复权限的精确边界

要求一个只输出5维边界的CPTP通道$\mathcal N:M_2\to M_5$，存在CPTP左逆对**全部**旧码成立，并且$\mathcal N(I_2/2)=I_5/5$。这是不可能的；如果仅要求保一个参考态而不保全码，则本结论不适用。

复用全码精确可纠正通道的Kraus正规形：误差矩阵$K_i^\dagger K_j=\lambda_{ij}I_2$经对角化成为$K_i=\sqrt{p_i}V_i$，其中$V_i^\dagger V_j=\delta_{ij}I_2$。于是

$$
\mathcal N(I_2/2)=\sum_i\frac{p_i}{2}V_iV_i^\dagger
$$

的非零支持秩为$2r\le5$，所以至多4。以支持投影为二值读口立即有

$$
\frac12\|\mathcal N(I_2/2)-I_5/5\|_1\ge1-\frac45=\frac15.
$$

界可达：先映射$\rho\mapsto\rho\otimes I_2/2$到四维支持，再嵌入五维；对四维支持取部分迹，剩余一维任意定态解码，给全CPTP恢复及参考谱$(1/4,1/4,1/4,1/4,0)$。这是一条已知恢复结构的有限推论，不作为新恢复定理领取。完整旗标输出有足够支持而逃脱此限制，同时必须付出已列的资料代价。

## 4. 核验与明确剩余自由

[check.py](check.py)用既有Python/NumPy完成三组：从144维原算符核查Gauss支持和独立几何；核查全部720置换、全码程序链及能源；穷尽2维输入的四个矩阵单位、用$10\times10$互补Gram核查旗标谱与$R$，并核查尖锐$1/5$删除边界。默认只读比较[results.json](results.json)与源码绑定，`--write`仅首次创建，不覆盖已有结果。巨大稠密总矩阵由已证明的酉装饰分解替代，未丢弃任何物理输出。

数值固定$a=e=g=1$，总Hamiltonian严格下界$4-\sqrt2$。检验不扫描腿数或参数。源文件绑定、浮点容差与每项残差保存在结果中；全码、任意$R$和正性的数学依据是上文算符恒等式。

已得到的仅是这份有限有效边界的参考扩展合同。SU(2)字典、面积尺度、能源简并、$e/a$、全交换参考角色、准备和静态耦合仍为公开输入；参考角色未采用。未取得实际时空面积、公共热温标、几何热余额、boost流、同一物质源反作用或Einstein方程。三项GR硬交付仍未完成。停止线与下一方向见[NEXT](NEXT.md)。
