# 第706轮：保边界局部近似与真实记录的共同来源响应

日期：2026-10-02。接[705](../../research_note_705.md)及[706双向域入口](local_domain_entry.md)，复用[623](../../../archive_585_628/research_note_623.md)—[625](../../../archive_585_628/research_note_625.md)、[637](../../../archive_629_652/research_note_637.md)、[704](../../research_note_704.md)。[代码](../joint_local_history_source_limit.py)、[结果](../joint_local_history_source_limit_results.json)、[核验](../research_round_706_checks.json)、[全账](../unified_physics_condition_ledger_706.md)。

## 1. 本轮问题与结论

705已经给出保全部边界表示、只压缩重数的严格局部通道。本轮证明：将这些通道插入原完整固定图的真实记录过程，可以共同保留准备、记录及总阶数二的来源响应。新增连接是局部双向能量域与真实历史的接合，不是重复625的全谱有限系统定理。

|层次|地位|
|---|---|
|认知动机|区域、记录和物质来源须来自同一过程|
|继承模型|完整原群、全部CAR、玻色变量、相互作用及Gauss限制；固定图和正几何基点|
|近似选择|来源无关、保边界载体且保CAR宇称的局部重数通道，在声明的记录后插入|
|解析增量|零阶cq记录后态及准备／动力总阶数二标量响应共同收敛|
|数值诊断|两节点、256维、冻结链路与Higgs方向的径向／中性CAR子模型，保留跨边耦合|
|未领取|整体有限输出、新有效H、自身有效Gibbs态、自动装置或空间连续极限|

每次等待仍用原完整H，整个过程仍无限维。初态是原H的Gibbs态，局部通道不是重新热化。正常CP构造不证明宇宙会自行实施这些通道；认知装置设计继续后置。

## 2. 同一局部族与宇称补充

采用617同一切口、637固定区域比较能源。按[705宇称补充](parity_completion_705.md)，回填向量逐局部CAR宇称选择：

$$
M_{A,\lambda}=\bigoplus_{\sigma=\pm}M_{A,\lambda,\sigma},\qquad
K_Ae_{\lambda\sigma j}=a_{\lambda\sigma j}e_{\lambda\sigma j},\qquad
K_Ag_{\lambda\sigma}=\mu_{\lambda\sigma}g_{\lambda\sigma}.
\tag{1}
$$

σ不是新增边界电荷，不要求两侧宇称匹配。原跨边跳跃可改变局部宇称；局部Kraus为偶算符即可。705冻结稿保留，补充另存，不声称其原夹具已经验证混合宇称。

P_{A,N}保完整谱簇。主Kraus只有一个，不先测量全部λ：

$$
\Phi_{A,N}(X)=P_{A,N}XP_{A,N}
+\sum_{a_{\lambda\sigma j}>N}K_{\lambda\sigma j}XK_{\lambda\sigma j}^{\dagger},
\qquad
K_{\lambda\sigma j}=|g_{\lambda\sigma}\rangle\langle e_{\lambda\sigma j}|\otimes I_{V_\lambda}.
\tag{2}
$$

正常性、保迹、Gauss支持、远端／被动参考保持及比较能源矩控制继承705和宇称补充。N是数学精度参数，不是新认知阈值；每扇区有限，整体仍保所有边界载体。

## 3. 双向图范数与原共同域

直接继承706入口。每次通道用一份新的、来源和N无关的辅助空间，具有ok及全部尾标签：

$$
S_N\psi=P_{A,N}\psi\otimes|ok\rangle+
\sum_{a>N}K_{\lambda\sigma j}\psi\otimes|\lambda,\sigma,j\rangle,\quad
S_\infty\psi=\psi\otimes|ok\rangle,\quad
D_E|\lambda,\sigma,j\rangle=(a-\mu)|\lambda,\sigma,j\rangle.
\tag{3}
$$

D_E非负，ok能量零，群和宇称作用均平凡。令R=K_A+K_B+c≥1，每次把辅助能量加入R_j=R_system+Σ_{l≤j}D_{E_l}，则

$$
R_{j+1}S_N=S_NR_j,\qquad
R_jS_N^\dagger=S_N^\dagger R_{j+1},\qquad
R_{j+1}(S_N-S_\infty)\psi\longrightarrow0.
\tag{4}
$$

伴随图范数收敛也已证明：先在固定有限辅助标签验证伴随强收敛，再用范数≤1和交织延拓。仅有等距强收敛不足以推出此步。辅助能量是数学域控制，未加进物理H，也不代表免费环境。

623共同算符域经637的J输送到原匹配空间。固定充分大c_H，令A=J(H_F+c_H)J†≥1，则

$$
D(A)=D(R),\qquad
c_1\|R\psi\|\le\|A\psi\|\le c_2\|R\psi\|,\qquad
\|G_a\psi\|+\|C_{ab}\psi\|\le c_{ab}\|R\psi\|.
\tag{5}
$$

R包含全部节点束缚势和群Casimir；只控一个电能矩不够。此算符域结论不能从637单边形式下界单独推出。

基点U₀与A对易，通常**不与R对易**。623的原C¹含时演化定理及图范数等价，给固定时间、小来源邻域内U、U†、L_r、L_r†的一致图界。辅助能量B=ΣD_E与系统R对易，且R、B各不超过R+B，所以

$$
\|R_j(T\otimes I)\psi\|
\le(C_T+\|T\|)\|R_j\psi\|,\quad T=U,U^\dagger,L_r,L_r^\dagger;
\qquad
\|(G_a\otimes I)\psi\|\le c_a\|R_j\psi\|.
\tag{6}
$$

这里系统算符与B对易；没有假装T与R对易。来源Hessian也有同类界。

## 4. 同一真实记录过程

固定有限等待时间、原仪器及区域菜单。第j次记录后插入所选A或B通道，统一记S_{j,N}。定义

$$
V_{\mathbf r,N}[\epsilon\lambda]
=S_{q,N}L_{r_q}U_q[\epsilon\lambda]\cdots
S_{1,N}L_{r_1}U_1[\epsilon\lambda],\qquad
\rho(\xi)=\frac{e^{-\beta H(x(\xi))}}{\operatorname{Tr}_{\rm phys}e^{-\beta H(x(\xi))}}.
\tag{7}
$$

所有N共用原归一Gauss准备态。记录后不自动热化；原跨区域相互作用保留，两边不被改为独立演化。两支共用同一辅助表示：

$$
Z_{\mathbf r,N}(\xi,\epsilon)
=\operatorname{Tr}\!\left[
V_{\mathbf r,N}[\epsilon\lambda_-]^\dagger
V_{\mathbf r,N}[\epsilon\lambda_+]\rho(\xi)\right].
\tag{8}
$$

相同两支给真实记录概率，不同两支可为复数。N=∞的S仅添ok，回到原过程。本族不同于625全谱Galerkin族和703正转移对数族，不混用其有限H。

## 5. 实际弱二阶来源

由式(4)—(6)，有限零来源前缀及伴随在相应图范数强收敛且一致有界。一次G插入只把域向量送入Hilbert空间，之后使用有界强收敛。因此

$$
V_{\mathbf r,N}(0)\psi\longrightarrow V_{\mathbf r,\infty}(0)\psi
\quad\text{于图范数},\qquad
V_{\mathbf r,N}^{(1)}\psi\longrightarrow V_{\mathbf r,\infty}^{(1)}\psi
\quad\text{于Hilbert范数}.
\tag{9}
$$

二阶有Hessian一次插入、两个有序G插入、双支各一次插入。沿625证明把较晚G移到bra一侧；两端是域向量，中间只有有界传播、记录和等距。来源插入后的向量无需再在域内。固定时间单纯形上支配积分得

$$
|b_{\mathbf r,k,N}(\chi,\psi)|\le C_k\|R\chi\|\|R\psi\|,\qquad
b_{\mathbf r,k,N}(\chi,\psi)\longrightarrow b_{\mathbf r,k,\infty}(\chi,\psi),
\quad k=0,1,2.
\tag{10}
$$

b是式(8)动力部分的标量形式；二阶记号按623—625的弱二阶Peano系数乘相应阶乘解释。不领取传播子强二阶可微、完整cq动态二阶迹范数可微，不增设共同D(H²)或来源后闭域。

以R⁻¹拉回：

$$
\langle u,C_{\mathbf r,k,N}v\rangle
=b_{\mathbf r,k,N}(R^{-1}u,R^{-1}v),\qquad
\sup_N\|C_{\mathbf r,k,N}\|<\infty,\qquad
C_{\mathbf r,k,N}\xrightarrow{\rm WOT}C_{\mathbf r,k,\infty}.
\tag{11}
$$

## 6. 同时接准备来源

704已有原热准备态的A夹权0／1／2阶导数。令B=RA⁻¹有界，直接继承

$$
T_j=R\rho^{(j)}(0)R
=B[A\rho^{(j)}(0)A]B^\dagger\in\mathfrak S_1,\qquad j=0,1,2.
\tag{12}
$$

这是夹权积的迹类延拓，不需ρ′为正。本轮不另近似准备态。式(10)的共同界允许配对准备导数和动力弱系数：

$$
\partial_\xi^j\partial_\epsilon^k Z_{\mathbf r,N}(0,0)
=\operatorname{Tr}(C_{\mathbf r,k,N}T_j)
\longrightarrow\operatorname{Tr}(C_{\mathbf r,k,\infty}T_j),
\qquad j+k\le2.
\tag{13}
$$

最后一步仅用T_j有限秩逼近、一致算符界和弱算符收敛。总阶数二包括值、两个一阶、两个纯二阶及混合项；原C²路径的曲率项依链式法则保留。Hessian由原H决定，不另拟合。

相同两支中，原效果≥I/4，局部通道保迹，q次记录有

$$
p_{\mathbf r,N}\ge4^{-q},\qquad
\partial_\xi\partial_\epsilon\log p_{\mathbf r,N}
=\frac{p_{\mathbf r,N,\xi\epsilon}}{p_{\mathbf r,N}}
-\frac{p_{\mathbf r,N,\xi}p_{\mathbf r,N,\epsilon}}{p_{\mathbf r,N}^{\,2}},
\qquad
\sum_{\mathbf r}p_{\mathbf r,N}^{(j,k)}=0\quad(j+k>0).
\tag{14}
$$

故归一得分、混合项和有限菜单Fisher量也共同匹配。不同CTP支的复系数不可套概率解释。

零阶cq记录／后态直接复用592／705：Φ_N在迹类上趋恒等，有限历史望远镜、CP收缩和部分迹收缩给

$$
\|\Omega_N-\Omega_\infty\|_1\longrightarrow0.
\tag{15}
$$

这一旧推论不单算新定理。来源还涉及G插入后的尾；初态概率尾小不保证二阶系数准。622／625量词保留，不承诺仅凭初始平均能源的普适二阶截止率，不增加图细化一致性。

## 7. 保跨边耦合的可复算诊断

两节点各有h、s各两个内部网格点，测度及逆度量沿625原径向二次型。每节点另保两中性CAR模，单节点16维，总256维。原M、束缚势、Majorana块及测地边势保留；链路和Higgs方向固定对齐，仅为条件诊断，不是全Gauss空间的正常固定配置态或已认证全谱近似。

采用655／667允许的中性自旋对角跳跃，声明系数.23。B节点CAR含A宇称串，核全CAR及总宇称；实际跳跃确实改变局部宇称。w=.8、ℏ=.7、β=1.2，共形来源按574／589／667原缩放：

$$
H(u)=e^{-6u}T+e^{6u}W+e^{2u}V+B+e^{-2u}J,\qquad
G=-6T+6W+2V-2J,\qquad C=36T+36W+4V+4J.
\tag{16}
$$

V为测地边势，B为Majorana项，J为实际跨边跳跃。先等待.12并读A，再等待.19并读B，每次读后作局部通道。全过程使用完整诊断H及原L±(s)，不重新热化。

局部参照采用单节点T+W固定谱基，按宇称分别回填最低向量；边界载体在本诊断取平凡部门。非Abel载体身份已由705及706入口验证，未借小模型替代。保存导数顺序为

$$
(j,k)=(0,0),(1,0),(0,1),(2,0),(1,1),(0,2).
\tag{17}
$$

两组检查通过：

- CAR、保迹、偶Kraus及独立Kraus求和为零残差，远端边缘差≤5.38×10⁻¹⁷。比较二阶矩4.61913932降为4.30774683，不称完整H能量下降。
- 跨边跳跃范数.46、测地边势范数.19234032均非零。原Hessian独立差分残差4.96×10⁻⁵。
- 请求1、2、3、4径向谱带，经完整简并簇处理得到局部主秩8、8、16、16。实际只有两个不同近似，不呈现为四个独立收敛数据点。
- 主秩8相对16的六类最大分支误差约2.66×10⁻⁹、6.84×10⁻⁸、.02267655、1.29×10⁻⁶、.60022232、.55662092；cq态误差.04335386。概率很准但动力来源仍未准。
- λ₊=1、λ₋=−.4时，混合系数总和由局部近似约5.61504002i变为完整诊断8.00590103i，复值不是概率。
- 独立本征指数四点混合差分，步长.0008／.0004的误差2.65×10⁻⁵／6.62×10⁻⁶。相同两支全部六类归一身份残差≤3.27×10⁻¹⁴。
- 删局部通道会令二阶来源变化.55662092，不能直接抄无通道来源。恢复完整16维后零差只是恢复诊断矩阵，不是无限维证书。

无穷维结论由式(3)—(15)承担；无空间连续数值模拟、图像检验或高精度扫描。

## 8. 旧空间接口逐项复用

|旧轮次|直接复用|本轮边界|
|---|---|---|
|382|完整反向方向接口三维上界及条件性三维|不把任意Bloch球当位置壳|
|383|连续自由对合可替代标准对径|不恢复标准对径要求|
|384|定性下界已删额外Lipschitz，传递可放宽为极小性|不恢复已删条件|
|386|渐近位移接真实邻域；有限尺度证书有覆盖／变化模预算|实际端点仍须与本过程相接|
|425|一致半幅／成本收缩给局部光滑坐标|与386为替代路线，不叠加必要性|
|522—523|热绝对差下界与实际方向仪器有共同条件实现|需识别与当前h、s、CAR／Gauss相同的参考及仪器|
|524|局域紧支撑探针及UV有限噪声已完成|不将523旧硬UV问题恢复为空白|

新增连接为

$$
\{\text{705局部保边界／偶通道},\ \text{706入口双向图域},\
\text{623—625原响应},\ \text{704夹权热导数}\}
\Longrightarrow
\text{同一局部近似的真实记录与二阶来源共同极限}.
\tag{18}
$$

补的是C02／C03／C04／C19／C20／C22固定图接口，没有补空间维数、物理体积、统一光锥或Einstein作用。内部λ和辅助E不当坐标，来源变化不等于动态量子几何反馈。

## 9. 下一项

“新局部族的来源是否存在”已接通，代价是无限边界载体及原完整H。705有限总容量权衡不撤销；相对熵仍沿705／638，不增加面积律。

返回实际共同空间尺度：区分区域合并、图细化和固定图谱精度。先回查637—646、663—667及旧空间接口，寻找已有局部代数、参考和传播合同可接合之处；不重复热迹、局部矩或容量反例，不另起认知设计。见[707入口](../../707/drafts/STATUS.md)。四分支和699指定候选的反例范围保持，目标不变。
