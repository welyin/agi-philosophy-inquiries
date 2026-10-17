# 边界大小与时间偏置有限更新：独立审阅

2026-10-17。**PASS_WITH_SCOPE**。没有必修项；通过范围为给定经典EH因果菱形约化字典的有限保辛检查，亦即代码中的PASS_CONDITIONAL_FINITE_BOUNDARY_DICTIONARY。未新增科学轮次、认知角色、公理或实际仪器资格。

## 1. 来源及模型范围

已直接读取[Bub等2408.11094v2](https://arxiv.org/html/2408.11094v2)与[650](../../../archive_629_652/research_note_650.md)、[651](../../../archive_629_652/research_note_651.md)。首篇采用EH协变相空间、球对称与规范限制，并给定$\kappa_0=\kappa_0(L)$后取得式(21)的辛结构；式(22)的括号符号与本稿相容。本项既不是从治理机制导出这份辛结构，也不是重新证明一般边界可积性定理。

文稿附列的另外三篇扩展工作是主线登记的对照线索；本独审不以其未在本次独立复核的量子或热学扩展作为有限计算的证明前提，亦不据其授予本项目任何实际权限。

## 2. 连续公式的独立核对

在$L>0$、$G>0$固定的完整二维相空间上，

$$
\Omega=\frac LG\,d\alpha\wedge dL=d\alpha\wedge dq,
\qquad q=\frac{L^2}{2G},\qquad\{q,\alpha\}=1.
$$

对$X=v(L)\partial_L$，直接取场空间外微分给

$$
\iota_X\Omega=-\frac{Lv}{G}\,d\alpha,\qquad
\mathcal L_X\Omega=\frac{(Lv)'}G\,d\alpha\wedge dL.
$$

故$\dot\alpha=0$时，保辛等价于$Lv=c$；$H=c\alpha/G$给$q'=q+cs/G$。正文正确保留$\delta\alpha$，没有将一维固定$\alpha$曲线上恒为零的拉回形式拿来作非退化检验。

对共同输送映射，直接有

$$
d\alpha'=\frac{L}{L+s}\,d\alpha
+\frac{\alpha s}{(L+s)^2}\,dL,\qquad dL'=dL,
$$

因此$\widetilde F_s^*\Omega=\Omega$。按正文的Hamilton约定，

$$
H=\frac{\alpha L}{G}
\quad\Longrightarrow\quad
\dot L=1,\qquad\dot\alpha=-\frac{\alpha}{L},
$$

故$L'=L+s$且$\alpha'L'=\alpha L$，符号正确。反之只改半径、固定偏置的拉回比确为$(L+s)/L$。

这些有限流仅在初末及复合所需的共同正半径域内使用。正文没有宣称它们在整个正半轴构成完备双向流或全局量子酉。常数面积步长与常数半径步长只在指定起点匹配大小终点，并非全部$L$上的同一映射。

## 3. 实际复算与有限量

默认只读执行作者检查，回执b29072，exit 0，read_only_match，54项精确分数检查通过；所列六份历史来源和原目标哈希由该入口核验，作者文件未被改写。

另外独立以有理数核验有限不变量、复合及积分，回执a54ec4，exit 0，未导入作者实现。得到

$$
b=\frac98,\qquad
\alpha_{\rm radius}=\frac{12}{35},\qquad
\alpha_{\rm area}-\alpha_{\rm radius}=\frac3{35},
\qquad
\frac{L_1}{L_0}=\frac54.
$$

对正文矩形，直接在原变量积分，初始辛体积为$5/2$；错误映射增加$s(3-2)/G=1/2$而成为$3$；共同输送的拉回密度仍为$L/G$，故体积保持$5/2$。这些数值不依赖数值求积或点样本外推。54项检查是已给符号身份的实现诊断，不是物理实验或全域证明的替代。

正文的本地引用已核，无非预期缺失链接；审阅及随后发布核验链接在文件生成前的占位合理。

## 4. 通过的方向结论及停止线

本项具体化旧650—651的共同边界运输要求：只比较最终面积不能认定完整边界后态或过程相同。但正确输送偏置后，固定半径步长仍可相容，因此没有唯一选择面积步长、大小更新律、$\kappa(L)$或温标。

保辛要求适用于声称在这份约化相空间内闭合的Hamilton更新。真实开放区域可以与被省略的物质、控制及边界部门交换辛通量；约化不保辛不构成对这类操作的普遍否决。正文已明确这一例外。

$\alpha$不是已认证的钟读数，参数$s$不是未经核对的实际固有时；大小变化也不自动产生曲率。本项没有物质与记录的完整交换通道、未知伙伴资格、有限实际读取误差或共同热学身份。其净值是旧原则的新文献例化与有限字典诊断，不能改称治理强猜想被证伪或EH被独立选出。停止增加该候选的器件与无截断量子化工作、保留原目标active及CT未采用，范围恰当。

## 5. 最终文件绑定

| 文件 | SHA256 |
|---|---|
| [boundary_charge_screen.md](boundary_charge_screen.md) | bf8a9936ebadf5e1528fa2abe52f63bb3cb6c9b00cdd0d8977be1c253ab63956 |
| [boundary_charge_check.py](boundary_charge_check.py) | 9ee2cef60953aef07e22607d8f261f67ce55646213250138e2b76a4b01517be9 |
| [boundary_charge_results.json](boundary_charge_results.json) | b828d8f8becd8211c37629523297fb9bfa1d159c01142dcb3d953c4ed86b718c |
| [boundary_charge_sources.json](boundary_charge_sources.json) | e09b8670e516a9995d6134f9c768cced2a0a02cb1893be937376714df51b64af |

不自绑定本审阅，不预绑定发布核验。**最终结论：PASS_WITH_SCOPE。**
