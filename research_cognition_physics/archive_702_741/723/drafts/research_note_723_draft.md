# 第723轮：完整物质参考记录与同一有限量子过程

日期：2026-10-04。接[722](../../research_note_722.md)、[723入口](joint_reference_entry.md)。[代码](../joint_reference_process_transport.py)、[结果](../joint_reference_process_transport_results.json)、[核验](../research_round_723_checks.json)、[条件账](../unified_physics_condition_ledger_723.md)。三组复算；主代理审查，无新增独立代理复核。

## 1. 问题及关闭的接口

722已给原Q_f的保能源读口；723入口显示，匹配四个参考各自的完整分布仍可能遗漏实际顺序记录。本轮将反例提升到原未截断Gauss模型，并把新仪器接入同一个有限CP过程。输出同时保记录、后态与能源，而不只拟合四条概率曲线。

正面结果分两种范围：

- **固定背景的有限历史：** 使用原H的真实等待，全部参考读口共用一组有限投影／CP尾项；完整记录后态在能源夹权迹范数中收敛，原能源、参考均值与末端来源形式共同收敛。
- **变化背景的立即联合读取：** 新Kraus及其一阶导数共用原能源域；原自身Gibbs准备、联合后态及几何一阶导数沿同一有限族收敛。

没有把第二项扩为含变化等待的二阶动态来源。625／704对其原仪器已成立的高阶结果保持；新仪器需要的额外域不能凭名称继承。

|层次|地位|
|---|---|
|认知动机|同一参考须支持后续操作及其共同代价，单次边缘不足以代表过程|
|继承|598／623原完整固定图H与Gauss；652四参考；722有理Kraus；625／704有限CP与热准备|
|实际任务输入|有限顺序菜单及等待；字段读口与722平方读口的实际相位|
|解析增量|原Gauss边缘相同／历史不同见证；新菜单的能源夹权共同过程及立即来源|
|数值|原一节点径向Gauss波函数积分；旧64维诊断的完整16记录和来源|
|仍开放|自主执行、空间局域有限近似、尺度一致、关系坐标逆及动态引力|

## 2. 保同一量子过程，而非只保四个量的分布

原菜单是T=h²/2、s、Q_T、Q_s；指定字段仪器与平方仪器为

$$
L_{f,\pm}=\sqrt{\tfrac12\pm\tfrac14\sin f},\qquad
K_{f,0}=\lambda_f(\lambda_f+iR_f)^{-1},\quad
K_{f,1}=I-K_{f,0},\qquad R_f=b_f-Q_f\ge0 .
\tag{1}
$$

f=T、s。字段函数有界、梯度有界，其乘法Kraus在原能源域有界；平方读口的同一性质及几何可微性由722复用，λ_f选在其充分范围。顺序可以明确指定，不声称它们构成保持四个原边缘的联合锐测量。

这与成熟过程张量方法的操作观点一致：[Jørgensen—Pollock，式(1)及附录A](https://arxiv.org/html/1902.00315v2)把中途操作与真实演化共同列入输出。这里只复用这份对象选择，不移植其Gaussian环境或张量网络效率结论。554的有限canonical联合测量、652平方压缩以及708矩层级均不重新证明。

## 3. 原Gauss模型中的四边缘不足定理

取原允许的一节点无边分支，W_T=W_s=0，保全部CAR、Dirac和Majorana。考虑径向Gauss波函数乘空CAR真空的子空间；该空间在四个参考下不变，但不要求在原H下不变。径向复共轭C给

$$
CTC=T,\quad CsC=s,\quad CV_fC=-V_f,\quad CQ_fC=Q_f .
\tag{2}
$$

身份先在原光滑核上成立，再由652的自伴闭包得到谱身份。因此ψ与Cψ对这四个参考的**所有谱集合**概率相同，不仅是若干矩相同。

令R=b+V_T²、F=(2+sin T)/4，F是式(1)的字段效果。原速度场Y_T非零，且在开集上Y_TF≠0。原核上

$$
[R,F]=[V_T^2,F]
=-i\hbar\{V_T M_{Y_TF}+M_{Y_TF}V_T\}\ne0 .
\tag{3}
$$

其一阶微分主项非零，可用该开集中的实径向紧支撑函数检测，故不是纯形式上的非对易。设K=λ(λ+iR)⁻¹，K=US，其中S=λ(λ²+R²)⁻¹/²为单射、值域稠密，U=(λ−iR)(λ²+R²)⁻¹/²为酉。实际两次记录的效果E=K†FK满足

$$
E=S U^\dagger F U S,\qquad CEC=S U F U^\dagger S,
\qquad U^2=\frac{\lambda-iR}{\lambda+iR} .
\tag{4}
$$

若E=CEC，S的稠密值域给U†FU=UFU†，从而F与U²对易。Cayley函数在R的谱上单射，其Borel逆恢复R，故F应与R的谱投影对易，和式(3)矛盾。因此B=(E−CEC)/2是非零、C奇的有界自伴算符。

实径向光滑Gauss核稠密，可选实、正交、归一的a、b使Im〈a,Eb〉≠0。取

$$
\psi_\pm=(a\pm ib)/\sqrt2,\qquad
\mu_{X,\psi_+}=\mu_{X,\psi_-}\ (X=T,s,Q_T,Q_s),\qquad
p_+-p_-=-2\operatorname{Im}\langle a,Eb\rangle\ne0 .
\tag{5}
$$

两态为正规Gauss有限能源态。原玻色动能／势在此径向表示为实形式；空CAR真空上的Dirac与Majorana期望均为零，尽管相应算符全部保留。因此两态的**原完整平均H**也相同。选择722充分大的λ，使两次实际读取的全部后态仍有限能源。

这证明原模型中四个单参考分布及平均能源不足以决定原顺序记录。它不排除完整量子态、限定任务摘要或近似描述。共轭在这里只用于一节点径向见证，不把它宣布成含任意复Yukawa的全H时间反演对称性。

## 4. 一个可直接积分的原模型见证

在722的ξ、z、S³半密度表示中，取w=1、b=1、ℏ=0.7、λ=64。a为ξ中心−0.5、宽度0.65的归一实Gaussian，b为同宽的一阶实Hermite函数；共同乘归一(1−4z²)²的紧支撑z函数、常数角谐波与空CAR真空。两波形正交；它们的Gaussian尾部使原能源有限。

一节点R=1−ℏ²∂ξ²，实际K不是新的Hamiltonian。其整线核为

$$
(Ku)(\xi)=\frac{-i\lambda}{2\hbar^2\kappa}
 \int_{\mathbb R}e^{-\kappa|\xi-\eta|}u(\eta)\,d\eta,
\qquad \kappa=\sqrt{1-i\lambda}/\hbar,\quad\operatorname{Re}\kappa>0 .
\tag{6}
$$

本参数Re κ≈8.14460，大于原负ξ端能源所需的最大指数速率M/w=2。卷积及其一阶导数有同样的指数尾界，z支撑保持；正ξ端的原动能／势只增加多项式权。因此这些**指定态**的K和I−K后态直接有限能源，不靠有限网格证明，也不把λ=64当成已计算出的任意全图充分常数。

原五维测度积分给：两初态平均H都约39.687017385；Q后总平均H分别约38.804027168、41.243314197。再读字段T时，指定两次结果的联合概率约0.574609776801、0.570362982046，差0.004246794755。一般非热准备的平均能源差可以有符号。

数值使用整线核的Fourier乘子作截断积分，1024／2048节点的记录差≤1.12×10⁻¹⁶。网格校准不能替代尾部证明，也不是完整图或完整热谱模拟。

## 5. 对新菜单使用能源夹权的共同近似

固定原背景γ₀，A=H(γ₀)+c≥1，𝓔=D(A¹/²)。对态或有符号迹类输入定义

$$
\|X\|_{A,1}=\|A^{1/2}XA^{1/2}\|_1,\qquad
\|\Omega\|_{A,1}=\sum_{\boldsymbol r}
 \|A^{1/2}\Omega_{\boldsymbol r}A^{1/2}\|_1 .
\tag{7}
$$

第二式把所有真实结果块保留。它比仅比较概率更强；对正态，范数就是总移位能源。

沿625取包含完整简并簇的P_Λ=1_{A≤Λ}、Z_Λ=I−P_Λ，以及固定最低Gauss本征向量g、Ag=a₀g；可在简并簇中同时选宇称特征向量。使用同一通道

$$
\mathcal C_\Lambda(X)=P_\Lambda X P_\Lambda
+\operatorname{Tr}(Z_\Lambda X)|g\rangle\langle g| .
\tag{8}
$$

这不是新压缩方法；新问题是它对722菜单所需的能源范数是否足够。令Y=A¹/²XA¹/²，夹权后的通道明确为

$$
\widetilde{\mathcal C}_\Lambda(Y)
=P_\Lambda YP_\Lambda
+a_0\operatorname{Tr}(Z_\Lambda A^{-1/2}YA^{-1/2})|g\rangle\langle g|,
\qquad \|\widetilde{\mathcal C}_\Lambda\|_{1\to1}\le1,
\quad \|\mathcal C_\Lambda X-X\|_{A,1}\to0 .
\tag{9}
$$

收缩来自CP和不增正输入的迹；收敛来自迹类谱尾及P_Λ强收敛。可逐被动参考恒等扩展：尾项保其边缘，不需要把未知参考重置。此处没有给无资源约束的统一误差率。

对式(1)任一个实际菜单𝓜(X)=⊕_r K_rXK_r†，由722及字段梯度界，设N_r=‖K_r‖_{𝓑(𝓔)}，则

$$
\|\mathcal M X\|_{A,1}\le\Big(\sum_r N_r^2\Big)\|X\|_{A,1} .
\tag{10}
$$

所以这些记录通道在共同夹权迹空间上有界。不必先证明K保持D(A)或D(A²)，也不把𝓔界误写成这种更强的域结论。

## 6. 有限过程必须保尾项、实际相位与原等待

对P_Λ支撑的输入，有限分支正是原K后接式(8)：

$$
\mathcal J_{r,\Lambda}(X)=B_{r,\Lambda}XB_{r,\Lambda}^\dagger
+\operatorname{Tr}(D_{r,\Lambda}X)|g\rangle\langle g|,\quad
B_{r,\Lambda}=P_\Lambda K_rP_\Lambda,\quad
D_{r,\Lambda}=P_\Lambda K_r^\dagger Z_\Lambda K_rP_\Lambda\succeq0 .
\tag{11}
$$

Σ_rJ_r,Λ保迹，且仍有原结果标签。D一般非零，不能以压缩后K的单一产品替代。这里压缩原实际K；不是先平方有限速度再自定义另一种读口。有限表示可依赖原算符数据，不声称已经给出高效计算算法或自主复位装置。

先固定背景与H。P_Λ与原U(t)=exp(−itH/ℏ)对易，有限等待U_Λ=exp(−itP_ΛHP_Λ/ℏ)在有限子空间恰为原U的限制。对任意有限历史，记

$$
\Omega_{\boldsymbol r}=K_{k,r_k}U_{k-1}\cdots K_{1,r_1}\rho
 K_{1,r_1}^\dagger\cdots U_{k-1}^\dagger K_{k,r_k}^\dagger,
\qquad \Omega_{\boldsymbol r,\Lambda}
=\mathcal J_{k,r_k,\Lambda}\mathcal U_{k-1,\Lambda}\cdots
 \mathcal J_{1,r_1,\Lambda}(\rho_\Lambda) .
\tag{12}
$$

可取ρ_Λ=C_Λρ；也可取704已经控制的自身Gibbs准备。只要求其嵌回原空间后在式(7)中收敛。有限乘积展开、式(9)—(10)以及原等待的能源等距性给

$$
\sum_{\boldsymbol r}\|\Omega_{\boldsymbol r,\Lambda}-
 \Omega_{\boldsymbol r}\|_{A,1}\longrightarrow0 .
\tag{13}
$$

该极限同时保全部记录概率和完整末态，包括与被动参考的关联。它是固定有限历史、固定图的强收敛，不是统一覆盖无限历史或任意来源的diamond范数结论。

若O是原𝓔有界形式，便有

$$
\left|\operatorname{Tr}O(\Omega_{\boldsymbol r,\Lambda}-\Omega_{\boldsymbol r})\right|
\le\|A^{-1/2}OA^{-1/2}\|\,
 \|\Omega_{\boldsymbol r,\Lambda}-\Omega_{\boldsymbol r}\|_{A,1}\to0 .
\tag{14}
$$

这包括原H、652四参考的一阶均值、623末端几何来源G。所有有界末端效果也共同匹配。**末端来源期望不是含变化等待的总几何导数**，两者不能混称。652已控制Q的绝对一阶矩，不借此要求任意高阶矩。

## 7. 立即联合读取的几何一阶过程也可共同输送

现在让原正几何γ变化，P_Λ仍固定在γ₀。式(1)字段读口固定；722给平方Kraus在同一𝓔上C¹。原ρ(γ)及其有限自身Gibbs族按704在夹权迹空间C¹收敛。一次记录的导数是

$$
\mathcal M_\gamma'(X)=\bigoplus_r
 \{K_{r,\gamma}'XK_{r,\gamma}^\dagger+
 K_{r,\gamma}XK_{r,\gamma}'{}^\dagger\},
\qquad \mathcal M_\gamma\in C^1(\mathcal B(\mathfrak T_A)) .
\tag{15}
$$

𝔗_A表示式(7)的夹权迹空间。保持相同有限顺序、无等待插入，令Ω_γ=𝓜_{k,γ}⋯𝓜_{1,γ}ρ_γ。导数保留每个实际instrument以及准备：

$$
\Omega_\gamma'=
\sum_{j=1}^k\mathcal M_{k,\gamma}\cdots\mathcal M_{j,\gamma}'\cdots
 \mathcal M_{1,\gamma}\rho_\gamma
+\mathcal M_{k,\gamma}\cdots\mathcal M_{1,\gamma}\rho_\gamma' .
\tag{16}
$$

有限族仍用式(11)，其中D的导数也由原K及K′共同产生。由于C_Λ与γ无关，式(9)对参数紧集上的连续迹类像一致收敛，逐项展开给

$$
\|\Omega_{\gamma,\Lambda}-\Omega_\gamma\|_{A,1}
+\|\partial_\gamma\Omega_{\gamma,\Lambda}-\partial_\gamma\Omega_\gamma\|_{A,1}
\longrightarrow0 .
\tag{17}
$$

因此不仅记录概率，立即联合后态本身的一阶来源也共同匹配。沿623原H、G有𝓔形式界，总加能及其导数可直接配对：

$$
\partial_\gamma\Delta E
=\sum_{\boldsymbol r}\left\{
 \operatorname{Tr}G\Omega_{\boldsymbol r}
 +\operatorname{Tr}H\Omega_{\boldsymbol r}'\right\}
-\operatorname{Tr}G\rho-\operatorname{Tr}H\rho' .
\tag{18}
$$

式(17)保证有限族的同一式(18)趋向原结果。新仪器在D(A)上的图范数控制以及含变化等待的高阶导数未在此证明；704对原旧菜单的二阶定理不被撤销，也不被无条件扩大。

## 8. 核验与实际限制

第二组采用旧64维径向诊断及722条件化邻居W，依次读T、Q_T、s、Q_s，中间等待0.19，保全部16个结果。每个有限Hamiltonian使用自身Gibbs准备；P固定，简并簇完整保留。字段、平方读口及来源共用同一原表示。

|有限秩|完整记录后态的能源夹权误差|末端来源期望误差|删CP尾项造成的概率损失|
|---|---:|---:|---:|
|16|9.411632143|31.522569894|0.088878846|
|32|6.564433142|18.322923581|0.059013283|
|48|2.908552647|6.904756807|0.017746988|
|64|0|0|舍入误差内为0|

这些误差并不小；不能以概率接近或“有限近似存在”冒称低秩已经准确。64为同一诊断的完整维数，其误差归零只是代数校准，不是原无限维极限的数值证明。原无限维结论由式(9)—(14)承担。

第三组在固定投影下同时改变原H、平方仪器与其自身Gibbs准备，检查立即四次读取。解析后态导数与独立差分的能源夹权误差≤1.64×10⁻⁸；完整64维的总加能导数为−17.637491082，独立差分为−17.637491085。未冻结仪器、尾项或准备。无图像检查。

## 9. 合并结果与停止条件

本轮将C03记录、C09参考、C19准备、C20同图有限表示、C21能源及C22来源进一步放入同一过程。四条边缘的独立拟合被原Gauss反例排除；保完整过程的正面路线成立。没有新增粒子、作用或认知公理，未生成三维、标准模型群、现实耦合或Einstein领先项。

至此停止单参考读口、尾部、低秩和二态例子的扩写。后续接[724](../../724/drafts/STATUS.md)：将实际物质参考接回647—651的关系区域和边界任务，先核哪些选面只需可共同读取的字段、哪些因果类型还需导数／混合资料；优先利用既有同一经典解和同一量子过程，不能以抽象联合记录代替实际空间区域。

382—386、425、522—524直接复用；384已消去的Lipschitz不重开，386／425保持替代桥。707／708的真实分块和任务相关保留条件不被当前全谱近似替代。固定图、指定连续物质、经典动态几何及辅助手征分支仍未自动相同，699限定失败保持，统一目标继续开放。
