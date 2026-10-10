# 1067之后的中心筛查：共同会合事件没有自动消掉关系部门

日期：2026-10-10。性质：下一轮准入前筛查，**计0；不创建编号研究，不改冻结文件或导航**。

## 1. 判断和术语

在同一有限维复量子载体上，两个投影复核 $P,Q$ 可以连续、对称地会合，双方互换时仍得到同一事件，而且所有过程保持同一完整任务代数；这些要求**仍不迫使任务代数中心为标量**。最小非交换反例为三维载体上的 $M_2\oplus\mathbb C$，没有经典部门的反例为四维载体上的 $M_2\oplus M_2$。

需要区分三个意思：

- **中心固定：** 连续且从恒等出发的任务代数自同构逐点固定其有限个中心投影。452及已有连续修订审计已得到这一点。
- **中心为标量：** $Z(\mathcal A)=\mathbb CI$，不存在多个不同简单部门。
- **表示不可约：** $\mathcal A'=\mathbb CI$，连同型重数也没有。中心为标量只给因子 $M_2\otimes I_K$，并不推出整个载体只有二维。

“固定每个部门”与“只存在一个部门”不同。共同事件也不自动删除事件前已有的未知资料和被动参考。

## 2. 回查与问题范围

已回读 [两投影来源筛选](direction_source_after1066.md)、[连续保判断修订](nondemolition_revision_after1066.md)、[1063](../research_note_1063.md)、[1064](../research_note_1064.md) 和 [452](../../archive_429_466/research_note_452.md)。

两投影的二维分块、持续过程固定有限中心、同型组织不能独自选qubit、单成员接触不等于整段任务二维均直接复用。本次只问：**加入对称共同会合和角色互换，是否会令多个中心部门不再可能？** 不把“共同事件”定义成“已经只剩一个块”再得循环论证。

下面给它一份明确且较强的内部操作解释：同一载体的两种投影通过相反的连续酉路径到达同一投影事件；存在可连续实现的角色互换，交换两个原判断及相应路径，固定会合事件；所有操作均属于或归一化同一任务代数。这里没有先生成真实空间的端点身份，这只是内部关系接口的筛查。

## 3. 对称会合的两部门反例

在一个二维块中取

$$
P_0=(I+Z)/2,\qquad Q_\beta=(I+\sin\beta X+\cos\beta Z)/2,
\quad0<\beta<\pi.
$$

令

$$
A_\beta(t)=e^{-it\beta Y/4},\qquad
P_\beta(t)=A_\beta(t)P_0A_\beta(t)^\dagger,
\quad Q_\beta(t)=A_\beta(t)^\dagger Q_\beta A_\beta(t).
$$

两方各移动一半，$t=1$ 时严格会合于

$$
M_\beta=\frac{I+\sin(\beta/2)X+\cos(\beta/2)Z}{2}.
$$

角色互换用

$$
S_\beta=\sin(\beta/2)X+\cos(\beta/2)Z,
\qquad S_\beta^2=I.
$$

它满足 $S_\beta P_0S_\beta=Q_\beta$、$S_\beta Q_\beta S_\beta=P_0$、$S_\beta M_\beta S_\beta=M_\beta$，并将 $A_\beta(t)$ 变成其逆。因此它不仅交换两个末端名字，也交换整个对称会合过程。$e^{-i\pi tS_\beta/2}$ 给从恒等通道到该互换通道的连续路径；整体相位无影响。

现在在 $\mathcal H=\mathbb C^2\oplus\mathbb C^2$ 上取

$$
P=P_0\oplus P_0,\quad Q=Q_{\pi/3}\oplus Q_{2\pi/3},
\quad M=M_{\pi/3}\oplus M_{2\pi/3},
\quad S=S_{\pi/3}\oplus S_{2\pi/3}.
$$

直和后的连续路径仍有全部上述对称、会合、互换和保任务性质。两个部门都非交换。可是

$$
D=(P-Q)^2=\tfrac14 I_2\oplus\tfrac34 I_2,\qquad
Z_1=\frac{\tfrac34I-D}{\tfrac12}=I_2\oplus0,
\quad Z_2=I-Z_1
$$

均属于 $\mathcal A=C^*(I,P,Q)$，且为非平凡中心投影。每块中的 $P_0,Q_\beta$ 生成 $M_2$，故

$$
\mathcal A=M_2\oplus M_2,\qquad Z(\mathcal A)=\mathbb C^2.
$$

没有借外界另放一个标签：部门与关系角已经在同一个量子载体内。真实仪器可取原来的Lüders读取，未知后态和参考全部保留；两投影相位、上述路径及互换都不扩大这个代数。

这些部门有后续操作意义。在两个部门的P成功态上再读Q，概率分别为 $3/4$ 和 $1/4$。即使都已获得同一个会合事件M的成功记录，再读Q的概率仍分别为

$$
\cos^2(\pi/12)=\frac{2+\sqrt3}{4},\qquad
\cos^2(\pi/6)=\frac34.
$$

所以“双方达成了同一个二元事件”不保证所有未知输入的未来行为相同。只保当前事件记录并丢弃部门，会损失实际继续能力。

允许一个经典部门时，$P=P_0\oplus0$、$Q=Q_\beta\oplus0$ 在 $\mathbb C^3$ 上已给 $M_2\oplus\mathbb C$。会合取 $M_\beta\oplus0$，互换取 $S_\beta\oplus(-1)$，路径按块延拓；$(P-Q)^2/\sin^2(\beta/2)=I_2\oplus0$ 显示中心。三维是包含非交换块且仍有非平凡中心的最小载体维数；若要求每个部门非交换，最小为四维。

以上反例不否定“共同事件具有更强含义”的可能性，只说明上述公开操作版本不足。若“同一完整任务”另指全部 $B(\mathcal H)$ 已不可约，那是在输入端加入结论，不能记成对称会合的推论。

## 4. 一项确有选择力的新认知候选

按[研究方法论](../../../猜想/从认知问题到社会机制再到物理检验的研究方法论.md)，社会对应可以是：**同一类复核结果交给下一方时，规则不应暗中随未声明的关系部门改变。** 这个要求比“双方互换公平”强，也并非所有社会判断都必须满足。它适用于声称有统一转接语义的指定接口。

可把它写成精确、可被反例检验的“结果内统一转接率”合同：对所有支持在P部门的输入，Q成功率都是同一 $a$；对所有支持在 $P^\perp$ 部门的输入，Q成功率都是同一 $b$，且 $0<a,b<1$。即

$$
PQP=aP,\qquad P^\perp QP^\perp=bP^\perp.\tag{U}
$$

这里仍要求P、Q是sharp投影，且两种P结果部门非零。它不是只对两份特定准备测到相同频率，也不是只说正反交接概率对称；上面的两块反例本来就逐块互惠。

**自含结论。** 按 $P\mathcal H\oplus P^\perp\mathcal H$ 写

$$
Q=\begin{pmatrix}aI&B\\B^\dagger&bI\end{pmatrix}.
$$

由 $Q^2=Q$ 得

$$
BB^\dagger=a(1-a)I,\quad B^\dagger B=b(1-b)I,
\quad(a+b-1)B=0.
$$

前两项严格正定，所以B同时单射、满射，两个部门维数相同。又 $B\ne0$，故 $b=1-a$，且 $B=\sqrt{a(1-a)}V$，V为部门间酉映射。经块内换基，

$$
P\cong\begin{pmatrix}1&0\\0&0\end{pmatrix}\otimes I_K,
\quad
Q\cong\begin{pmatrix}a&\sqrt{a(1-a)}\\\sqrt{a(1-a)}&1-a\end{pmatrix}\otimes I_K,
\quad\mathcal A\cong M_2\otimes I_K.
$$

因此(U)可以推出中心为标量，**不必先写中心为标量**；$a$ 也无须指定成1/2。但K和未知相关资料仍存在，不能由此宣布完整载体二维。若想进一步消去K，需要另一份对全部相关能力与资料的审计。

这项合同独立于对称会合：第3节反例会合无问题，却因两个部门的转接率不同而违反(U)。它有具体可观察的失败见证，而不是仅给不可约性换一个名字。不过它本身尚未由共识必然性推出；精确全输入合同的实验认证、有限窗口近似版以及实际仪器权限均另需说明。微小不等角可以使精确中心不平凡，不能把有限误差率相等直接当作精确中心已消失。

## 5. 成熟来源及去重边界

已核读[Böttcher、Simon、Spitkovsky，*Similarity between two projections*](https://arxiv.org/pdf/1705.08937)的导言、Proposition 2及§3。该文明确研究同时互换两个投影的酉，且给出一般两投影块公式；互换本身并不要求单一主角。本文件第3节把这一成熟现象写成当前共同事件问题的显式反例，不声称新发现投影互换分类。

另核读[Kribs、Mammarella、Pereira，*Isoclinic Subspaces and Quantum Error Correction*](https://arxiv.org/pdf/1912.10100)的§2、Theorem 4：统一投影夹乘系数与等主角结构有成熟对应。本文件用上面两部门块计算直接证明(U)的结论，不依赖该文记号中的角系数约定；这里 $a$ 是概率，即相应主角的余弦平方。

正式报告检索“PQP／isoclinic／等倾”等，未见(U)作为当前接口认知来源的同等应用；820旧二维主角方法直接复用，不重复计新数学。452的有限中心自同构事实也不另立轮。

## 6. 可复算有限检查

使用已有Python与NumPy执行了下段代码。两个部门都通过对称会合、互换、共同事件固定；中心投影秩2，D的谱为 $(1/4,1/4,3/4,3/4)$。两份P成功后的Q概率为 $0.75,0.25$；共同M成功后的Q概率为 $0.9330127018922193,0.75$。最大矩阵残差 `3.269012620426328e-16`。这是有限反例核算，非普遍分类的数值证明。

```python
import numpy as np
I=np.eye(2); P=np.diag([1.,0.])
X=np.array([[0.,1.],[1.,0.]])
Y=np.array([[0.,-1j],[1j,0.]])
Z=np.diag([1.,-1.]); res=[]; qs=[]; ms=[]
for beta in [np.pi/3,2*np.pi/3]:
    Q=(I+np.sin(beta)*X+np.cos(beta)*Z)/2
    S=np.sin(beta/2)*X+np.cos(beta/2)*Z
    M=(I+S)/2
    A=np.cos(beta/4)*I-1j*np.sin(beta/4)*Y
    res += [np.linalg.norm(S@P@S-Q),np.linalg.norm(S@Q@S-P),
            np.linalg.norm(S@M@S-M),
            np.linalg.norm(A@P@A.conj().T-M),
            np.linalg.norm(A.conj().T@Q@A-M)]
    qs.append(Q); ms.append(M)
def bd(a,b):
    return np.block([[a,np.zeros((2,2))],[np.zeros((2,2)),b]])
p=bd(P,P);q=bd(*qs);D=(p-q)@(p-q)
z1=(.75*np.eye(4)-D)/.5
res += [np.linalg.norm(z1@z1-z1),np.linalg.norm(z1@p-p@z1),
        np.linalg.norm(z1@q-q@z1)]
print(np.linalg.eigvalsh(D),np.trace(z1),max(res))
print([Q[0,0] for Q in qs])
print([np.trace(M@Q) for M,Q in zip(ms,qs)])
```

## 7. 准入判断

本次可以保留两个清楚结论：对称会合／角色互换不独自消掉中心；统一转接率可在明确投影任务内推出单一M₂类型，而不预先假定它。两者都只是成熟结构在新认知候选上的筛查，本次计0，**不建议仅为中心块分类单开编号轮次**。

如果后续真的从共同操作与未知资料保护推出(U)，或以独立可检验条件选出原子事件及完整载体维数，那会改变三维主线的上游输入账，值得重新准入。现阶段继续优化部门角度、增加块数或反复验证同一会合矩阵没有选择收益。
