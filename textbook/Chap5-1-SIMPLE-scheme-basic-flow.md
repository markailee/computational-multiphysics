# 第5章 单场求解算法：SIMPLE 与基本流动计算

第4章讨论了单个通用输运方程的有限体积离散。对于实际流动，连续性方程、动量方程、能量方程以及湍流模型方程相互耦合，不能简单地把每个方程完全独立地求解。特别是不可压缩流动的压力没有独立的随时间演化的控制方程，它通过连续性约束影响速度场。因此，压力—速度耦合是流场求解的核心问题。

本章按由基础到复杂的顺序介绍压力基分离式求解方法。首先推导 SIMPLE 算法及压力修正方程，讨论交错网格、同位网格和动量插值；随后介绍 RANS、$k$–$\varepsilon$、$k$–$\omega$ 和 LES 湍流模拟；最后把能量方程、湍流传热、热膨胀及可压缩性纳入统一的迭代框架。

## 5.1 压力基求解方法：SIMPLE 算法流程及动量插值

### 5.1.1 联立求解与分离求解

离散后的流场方程可以采用两类基本求解方式。

- **联立式求解**：将全部或部分未知量的代数方程组成一个大型系统同时求解。变量之间的耦合较强，但矩阵规模大、结构复杂，求解成本和内存需求较高。
- **分离式求解**：按照一定顺序求解若干较小的问题，通过反复迭代实现变量耦合。每个子问题通常对应一个速度分量、压力或其他标量，便于复用第4章的稀疏求解器。

分离式求解又包括原始变量法和非原始变量法。涡量—流函数方法通过交叉求导消去压力，并引入流函数自动满足二维不可压缩连续性方程，但其三维推广和边界处理不够方便。原始变量法则直接求解速度和压力，其中 SIMPLE、SIMPLEC 和 PISO 等压力修正算法应用最广。

### 5.1.2 控制方程与压力—速度耦合

可压缩或变密度牛顿流体的连续性和动量方程可写为

$$
\frac{\partial\rho}{\partial t}
+\frac{\partial(\rho u_k)}{\partial x_k}=0,
\tag{5.1}
$$

$$
\frac{\partial(\rho u_i)}{\partial t}
+\frac{\partial(\rho u_k u_i)}{\partial x_k}
=-\frac{\partial p}{\partial x_i}
+\frac{\partial}{\partial x_k}
\left[\mu\left(\frac{\partial u_i}{\partial x_k}
+\frac{\partial u_k}{\partial x_i}\right)\right]
+\frac{\partial}{\partial x_i}
\left(\lambda\frac{\partial u_k}{\partial x_k}\right)
+\rho f_i .
\tag{5.2}
$$

对于不可压缩流动，$\rho$ 为常数且 $\nabla\cdot\mathbf U=0$。压力在动量方程中以梯度出现，而连续性方程只含速度。若给定一个不满足连续性约束的速度场，就需要调整压力，使压力梯度驱动出的速度修正消除质量不平衡。

### 5.1.3 压力棋盘振荡与网格布置

以一维、常系数动量方程为例，若在同一个节点集合上存储压力和速度，并用跨越一个节点的中心差分近似压力梯度，则

$$
\left(\frac{\mathrm dp}{\mathrm dx}\right)_i
\approx\frac{p_{i+1}-p_{i-1}}{2\delta x}.
\tag{5.3}
$$

该梯度不包含当前节点压力 $p_i$。在平滑压力场上叠加形如

$$
p_i=p_i^{\mathrm{smooth}}+(-1)^i\widehat p
\tag{5.4}
$$

的交替压力扰动，离散压力梯度可能完全不变。这样的压力振荡不能被常规动量方程识别，称为**棋盘压力振荡**。

解决这一问题的一种经典方法是采用**交错网格**：标量压力位于控制体中心，速度分量位于相应方向的控制体界面上。东侧界面速度 $u_e$ 两侧的压力正好是 $p_P$ 和 $p_E$，因此压力梯度由相邻的 $1$–$\delta$ 压差直接构造：

$$
-\int_{V_e}\frac{\partial p}{\partial x}\,\mathrm dV
\approx A_e(p_P-p_E).
\tag{5.5}
$$

交错网格能够自然地抑制棋盘振荡，并使连续性方程中的界面速度直接可用，但在曲线坐标和非结构化网格上实现较复杂。同位网格则把所有变量存储在单元中心，需要借助动量插值恢复界面压力—速度耦合。

### 5.1.4 动量方程的试探速度解

在交错网格上，东侧速度控制体的离散动量方程可写为

$$
a_e u_e=\sum_{NB}a_{NB}u_{NB}+b_e+A_e(p_P-p_E).
\tag{5.6}
$$

系数 $a_e$、$a_{NB}$ 和源项 $b_e$ 由当前迭代的速度、密度、黏度、时间步长以及对流—扩散格式决定。用当前压力估计 $p^*$ 求解动量方程，得到试探速度：

$$
a_e u_e^*=\sum_{NB}a_{NB}u_{NB}^*+b_e
+A_e(p_P^*-p_E^*).
\tag{5.7}
$$

一般情况下，$u_e^*$ 不满足离散连续性方程。定义修正量

$$
p=p^*+p',\qquad u_e=u_e^*+u_e'.
\tag{5.8}
$$

假设修正后的速度和压力仍满足同一层次的线性化动量方程，则

$$
a_e(u_e^*+u_e')
=\sum_{NB}a_{NB}(u_{NB}^*+u_{NB}')+b_e
+A_e(p_P^*+p_P'-p_E^*-p_E').
\tag{5.9}
$$

减去式（5.7），得到严格的速度修正关系

$$
a_eu_e'=\sum_{NB}a_{NB}u_{NB}'+A_e(p_P'-p_E').
\tag{5.10}
$$

若直接求解式（5.10），所有邻点速度修正将全场联立，失去了分离求解的优势。SIMPLE 的核心近似是**忽略邻点速度修正的影响**，从而得到

$$
u_e'=d_e(p_P'-p_E'),
\qquad d_e=\frac{A_e}{a_e}.
\tag{5.11}
$$

因此，改进速度为

$$
u_e=u_e^*+d_e(p_P'-p_E'),
\qquad
v_n=v_n^*+d_n(p_P'-p_N'),
\qquad d_n=\frac{A_n}{a_n}.
\tag{5.12}
$$

压力修正一旦求得，速度便可用局部代数关系更新。下一步是利用连续性方程构造 $p'$ 的方程。

## 5.2 压力修正方程与 SIMPLE 算法

### 5.2.1 压力修正方程的推导

对固定二维控制体积分连续性方程，并采用一阶隐式时间离散，可得

$$
\frac{\rho_P-\rho_P^n}{\Delta t}V_P
+(\rho u)_eA_e-(\rho u)_wA_w
+(\rho v)_nA_n-(\rho v)_sA_s=0.
\tag{5.13}
$$

将 $u=u^*+u'$、$v=v^*+v'$ 及式（5.12）代入。对于不可压缩流动，整理为标准的压力修正方程：

$$
a_Pp_P'=a_Ep_E'+a_Wp_W'+a_Np_N'+a_Sp_S'+b_P,
\tag{5.14}
$$

其中

$$
a_E=\rho_ed_eA_e,\quad
a_W=\rho_wd_wA_w,\quad
a_N=\rho_nd_nA_n,\quad
a_S=\rho_sd_sA_s,
\tag{5.15}
$$

$$
a_P=a_E+a_W+a_N+a_S,
\tag{5.16}
$$

$$
\begin{aligned}
b_P={}&-\frac{\rho_P-\rho_P^n}{\Delta t}V_P\\
&+(\rho u^*)_wA_w-(\rho u^*)_eA_e
+(\rho v^*)_sA_s-(\rho v^*)_nA_n .
\end{aligned}
\tag{5.17}
$$

右端 $b_P$ 是试探场的**负质量残差**。若试探速度已满足连续性方程，则 $b_P=0$。因此，压力修正方程的源项可直接用于构造质量守恒收敛指标。

式（5.14）在结构上与稳态扩散方程相似，但物理意义不同：它不是独立的物理压力输运方程，而是由动量修正关系和离散质量守恒共同导出的压力修正方程。

### 5.2.2 压力修正方程的边界条件

压力修正方程必须与流场边界条件一致。

**给定边界法向速度。** 边界速度已经确定，不需要修正，即 $u_b'=0$。在组装压力修正方程时，相当于令穿过该边界的压力修正扩散系数为零。这对应齐次 Neumann 型条件

$$
\frac{\partial p'}{\partial n}=0.
\tag{5.18}
$$

**给定压力。** 边界压力不再修正，因此

$$
p_b'=0.
\tag{5.19}
$$

这对应 Dirichlet 条件，也为压力场提供参考水平。

若整个流场边界都只给定速度，不可压缩压力仅由梯度决定，压力修正矩阵会具有常数零空间。此时需要在一个单元固定压力参考值，或采用等价的零均值约束。参考值不改变速度场，但会影响压力的绝对数值。

### 5.2.3 SIMPLE 算法流程

SIMPLE 是 *Semi-Implicit Method for Pressure-Linked Equations* 的缩写，即求解压力耦合方程的半隐方法。一个典型的稳态外迭代过程如下。

1. 给定压力 $p^*$ 和速度初场，计算物性、动量系数及松弛后的对角项。
2. 分别求解各方向动量方程，得到试探速度 $\mathbf U^*$。
3. 由 $\mathbf U^*$ 计算各单元质量残差，组装压力修正方程。
4. 求解 $p'$ 方程。
5. 修正速度，使其在当前修正关系下满足离散连续性方程。
6. 按 $p=p^*+\alpha_pp'$ 更新压力，其中 $\alpha_p$ 为压力亚松弛因子。
7. 求解与流场耦合的温度、湍流量或其他标量方程。
8. 更新系数和物性，检查质量、动量及其他方程残差；未收敛则返回步骤1。

速度更新通常直接采用式（5.12），压力则可亚松弛。动量方程本身也常在组装前施加亚松弛。松弛因子过小会使收敛缓慢，过大则可能造成外迭代振荡或发散。

```text
给定 p*、U* 和其他场的初值
重复执行外迭代
    更新物性和动量离散系数
    求解各方向动量方程，得到试探速度 U*
    计算界面流量和单元质量残差
    组装并求解压力修正方程，得到 p'
    用 p' 修正界面速度和单元速度
    p = p* + alpha_p p'
    求解能量、湍流及其他标量方程
    更新延迟项和物性
    检查残差、全局质量平衡和目标物理量
直到收敛
```

在一次修正后，速度满足的是采用 SIMPLE 近似构造的离散连续性方程，而动量方程因系数更新和压力亚松弛仍可能不完全满足。只有外迭代收敛、$p'\rightarrow0$ 且各方程残差均足够小时，才同时逼近离散动量和连续性方程的解。

### 5.2.4 同位网格上的动量插值

同位网格将速度和压力都存储在单元中心。以 $x$ 方向为例，离散动量方程写为

$$
a_P^u u_P=\sum_{NB}a_{NB}^u u_{NB}+b_P^u
-V_P\left(\frac{\partial p}{\partial x}\right)_P .
\tag{5.20}
$$

定义假拟速度 $\widehat u_P$ 和动量系数 $D_P^u$：

$$
u_P=\widehat u_P-D_P^u\left(\frac{\partial p}{\partial x}\right)_P,
\qquad
\widehat u_P=\frac{\sum_{NB}a_{NB}^uu_{NB}+b_P^u}{a_P^u},
\qquad
D_P^u=\frac{V_P}{a_P^u}.
\tag{5.21}
$$

若只是把 $u_P$、$u_E$ 和单元中心压力梯度分别线性插值到界面，界面速度仍可能无法识别交替压力扰动。Rhie–Chow 动量插值的关键是：**插值假拟速度和动量系数，但在界面上重新使用相邻单元压力构造 $1$–$\delta$ 压差。**

在均匀一维网格上，可写为

$$
u_e=\overline{\widehat u}_e-D_e^u\frac{p_E-p_P}{\delta x_e},
\tag{5.22}
$$

其中上横线表示由两侧单元插值。等价地，

$$
u_e=\overline u_e
+D_e^u\left[
\overline{\left(\frac{\partial p}{\partial x}\right)}_e
-\frac{p_E-p_P}{\delta x_e}
\right].
\tag{5.23}
$$

第二项消除了普通线性插值中的压力解耦，恢复与交错网格相似的界面压力响应。速度修正仍具有局部形式

$$
u_e'=d_e(p_P'-p_E'),
\qquad d_e=\frac{D_e^u}{\delta x_e}.
\tag{5.24}
$$

单元中心速度可根据压力梯度修正，界面压力修正值则可插值到单元中心。用于对流输运的界面质量流量必须由动量插值后的界面速度构造，而不能只由单元中心速度简单平均得到。

### 5.2.5 非结构化网格上的推广

在非结构化网格上，连续性方程写为

$$
\frac{(\rho_P-\rho_P^n)V_P}{\Delta t}
+\sum_{f\in\mathcal F_P}\rho_f\mathbf U_f\cdot\mathbf A_f=0.
\tag{5.25}
$$

令 $\mathbf d_f=\mathbf x_N-\mathbf x_P$、$d_f=\lVert\mathbf d_f\rVert$、$\mathbf e_f=\mathbf d_f/d_f$。一种常用的动量插值形式是

$$
\mathbf U_f=
\overline{\widehat{\mathbf U}}_f
-D_f\left[
\overline{\nabla p}_f+
\left(\frac{p_N-p_P}{d_f}
-\overline{\nabla p}_f\cdot\mathbf e_f\right)\mathbf e_f
\right].
\tag{5.26}
$$

在正交网格上取法向分量即可恢复式（5.22）。非正交项、曲率项和网格畸变修正可通过延迟方式处理。压力修正方程则由 $\mathbf U_f'$ 与压力法向梯度的关系代入式（5.25）获得。实现时应保证压力方程系数、界面速度修正和最终质量流量采用同一套几何与动量系数。

### 5.2.6 SIMPLEC 与收敛诊断

SIMPLE 忽略了 $\sum a_{NB}u_{NB}'$。SIMPLEC 通过近似 $u_{NB}'\approx u_e'$，把该项的一部分移至左端：

$$
(a_e-\sum_{NB}a_{NB})u_e'
=A_e(p_P'-p_E').
\tag{5.27}
$$

因此可采用

$$
d_e^{\mathrm{SIMPLEC}}
=\frac{A_e}{a_e-\sum_{NB}a_{NB}}
\tag{5.28}
$$

代替 SIMPLE 的 $d_e=A_e/a_e$。这一处理减轻了忽略邻点修正造成的影响，通常可减弱对压力亚松弛的依赖，但实际效果仍取决于系数定义、松弛方式和问题类型。

压力修正收敛不代表整个流场已经收敛。至少应同时检查：

- 单元及全局质量不平衡；
- 各方向动量残差；
- 压力和速度在外迭代之间的变化；
- 能量、湍流量及其他耦合方程残差；
- 流量、压降、升阻力、换热量等目标物理量是否稳定。

## 5.3 雷诺平均模拟与标准 $k$–$\varepsilon$ 模型

### 5.3.1 湍流的多尺度性与模拟层级

湍流具有非线性、多尺度、不规则脉动和强耗散等特征。大尺度涡从平均流动获得能量，能量通过级串向小尺度传递，最终在黏性尺度上转化为内能。在惯性子区，经典理论给出能谱的 $-5/3$ 波数标度。

按解析尺度的不同，主要数值方法包括：

| 方法 | 直接解析的内容 | 需要建模的内容 | 主要特点 |
|---|---|---|---|
| DNS | 全部湍流尺度 | 无湍流封闭模型 | 网格需达到 Kolmogorov 耗散尺度，成本极高 |
| LES | 大尺度及惯性区以上的主要涡结构 | 亚格子尺度 | 瞬态三维计算，适合研究涡结构和非定常过程 |
| RANS | 统计平均流动 | 整个湍流谱的影响 | 工程成本低，但依赖封闭模型 |

URANS 在 RANS 封闭基础上求解随时间变化的平均场，可以描述大尺度非定常运动，但并不等同于直接解析全部湍流脉动。

### 5.3.2 雷诺分解与 RANS 方程

将瞬时变量分解为平均量和脉动量：

$$
u_i=U_i+u_i',\qquad p=P+p',
\qquad \overline{u_i'}=0.
\tag{5.29}
$$

平均算子满足

$$
\overline{U_i}=U_i,\qquad
\overline{\varphi_1\varphi_2}
=\Phi_1\Phi_2+\overline{\varphi_1'\varphi_2'},
\tag{5.30}
$$

并且在通常的连续时间平均意义下可与时间、空间微分交换顺序。对不可压缩 Navier–Stokes 方程取平均，得到

$$
\frac{\partial U_i}{\partial x_i}=0,
\tag{5.31}
$$

$$
\rho\left(\frac{\partial U_i}{\partial t}
+U_j\frac{\partial U_i}{\partial x_j}\right)
=-\frac{\partial P}{\partial x_i}
+\frac{\partial}{\partial x_j}
\left[2\mu S_{ij}-\rho\overline{u_i'u_j'}\right]+S_i,
\tag{5.32}
$$

其中平均应变率张量为

$$
S_{ij}=\frac12\left(
\frac{\partial U_i}{\partial x_j}
+\frac{\partial U_j}{\partial x_i}\right).
\tag{5.33}
$$

新增项 $-\rho\overline{u_i'u_j'}$ 称为雷诺应力。该张量对称，共有六个独立分量，但 RANS 方程本身不能确定它们，这就是封闭问题。

### 5.3.3 Boussinesq 涡黏假设

Boussinesq 假设把雷诺应力的偏量部分与平均应变率联系起来：

$$
-\rho\overline{u_i'u_j'}
=2\mu_t\left(S_{ij}-\frac13S_{kk}\delta_{ij}\right)
-\frac23\rho k\delta_{ij},
\tag{5.34}
$$

其中

$$
k=\frac12\overline{u_i'u_i'}
=\frac12\left(\overline{u'^2}+\overline{v'^2}+\overline{w'^2}\right)
\tag{5.35}
$$

为湍动能，$\mu_t$ 为湍流或涡黏系数。不可压缩时 $S_{kk}=0$。将式（5.34）代入 RANS 方程后，偏应变部分表现为有效黏度

$$
\mu_{\mathrm{eff}}=\mu+\mu_t,
\tag{5.36}
$$

而 $-\frac23\rho k\delta_{ij}$ 的散度可并入修正压力。涡黏模型根据确定 $\mu_t$ 所需的附加方程数量，可分为零方程、一方程和双方程模型等。

### 5.3.4 标准 $k$–$\varepsilon$ 模型

标准 $k$–$\varepsilon$ 模型通过湍动能 $k$ 和耗散率 $\varepsilon$ 构造湍流时间与长度尺度：

$$
\mu_t=C_\mu\rho\frac{k^2}{\varepsilon}.
\tag{5.37}
$$

其输运方程为

$$
\frac{\partial\rho k}{\partial t}
+\frac{\partial\rho U_jk}{\partial x_j}
=\frac{\partial}{\partial x_j}
\left[\left(\mu+\frac{\mu_t}{\sigma_k}\right)
\frac{\partial k}{\partial x_j}\right]
+P_k-\rho\varepsilon,
\tag{5.38}
$$

$$
\frac{\partial\rho\varepsilon}{\partial t}
+\frac{\partial\rho U_j\varepsilon}{\partial x_j}
=\frac{\partial}{\partial x_j}
\left[\left(\mu+\frac{\mu_t}{\sigma_\varepsilon}\right)
\frac{\partial\varepsilon}{\partial x_j}\right]
+C_{\varepsilon1}\frac{\varepsilon}{k}P_k
-C_{\varepsilon2}\rho\frac{\varepsilon^2}{k}.
\tag{5.39}
$$

$P_k$ 为由平均剪切产生的湍动能生成项。采用涡黏封闭时，可写为

$$
P_k=-\rho\overline{u_i'u_j'}\frac{\partial U_i}{\partial x_j}.
\tag{5.40}
$$

常用高雷诺数常数为 $C_\mu=0.09$、$C_{\varepsilon1}=1.44$、$C_{\varepsilon2}=1.92$、$\sigma_k=1.0$、$\sigma_\varepsilon=1.3$。不同模型变体可调整常数并增加阻尼函数、可实现性约束或浮力项。

$k$ 方程具有清晰的生成—输运—耗散结构。$\varepsilon$ 方程虽可从脉动方程推导，但包含大量高阶相关项，工程形式主要通过量纲类比、实验和 DNS 数据封闭。

### 5.3.5 近壁处理与适用范围

标准 $k$–$\varepsilon$ 模型适用于充分发展的高雷诺数湍流，不宜直接积分到黏性底层。采用壁面函数时，第一层单元中心布置在对数区：

$$
U^+=\frac{U_t}{u_\tau}
=\frac1\kappa\ln(Ey^+),
\qquad
y^+=\frac{\rho u_\tau y_P}{\mu},
\tag{5.41}
$$

其中 $U_t$ 为平行壁面的速度，$u_\tau=\sqrt{\tau_w/\rho}$ 为摩擦速度。可用 $u_\tau=C_\mu^{1/4}\sqrt{k}$ 估算近壁湍流尺度，近壁耗散率常给为

$$
\varepsilon_P=\frac{C_\mu^{3/4}k_P^{3/2}}{\kappa y_P}.
\tag{5.42}
$$

标准模型稳健、成本低，在自由剪切流和许多工程流动中表现良好，但也存在明显限制：

- 近壁低雷诺数区域需要壁面函数或专门修正；
- 对强逆压梯度、分离、回流和曲率影响的预测有限；
- 各向同性涡黏假设难以描述复杂湍流应力；
- 在部分强梯度区域可能高估涡黏和扩散。

这些限制促使 $k$–$\omega$、SST、雷诺应力模型以及 LES 等方法的发展。

## 5.4 $k$–$\omega$ 与 SST 湍流模型

### 5.4.1 比耗散率与 Wilcox $k$–$\omega$ 模型

$k$–$\omega$ 模型以比耗散率 $\omega$ 替代 $\varepsilon$。两者的常用关系为

$$
\varepsilon=\beta^*k\omega,
\qquad
\mu_t=\rho\frac{k}{\omega},
\qquad
\ell\sim\frac{\sqrt k}{\omega}.
\tag{5.43}
$$

一种便于说明模型结构的守恒形式为

$$
\frac{\partial\rho k}{\partial t}
+\frac{\partial\rho U_jk}{\partial x_j}
=P_k-\beta^*\rho k\omega
+\frac{\partial}{\partial x_j}
\left[(\mu+\sigma_k\mu_t)\frac{\partial k}{\partial x_j}\right],
\tag{5.44}
$$

$$
\frac{\partial\rho\omega}{\partial t}
+\frac{\partial\rho U_j\omega}{\partial x_j}
=\frac{\gamma}{\nu_t}P_k-\beta\rho\omega^2
+\frac{\partial}{\partial x_j}
\left[(\mu+\sigma_\omega\mu_t)\frac{\partial\omega}{\partial x_j}\right]
+D_\omega ,
\tag{5.45}
$$

其中 $\nu_t=\mu_t/\rho$，$D_\omega$ 表示具体版本中的交叉扩散项。Wilcox 模型可直接积分到壁面附近，对边界层和逆压梯度通常比标准 $k$–$\varepsilon$ 更合适；其不足是对自由流 $\omega$ 边界值较敏感。

固壁处 $k=0$，$\omega$ 采用与第一层网格距离相适应的壁面值。进口和自由流中的 $k$、$\omega$ 不应随意置零，否则会产生过大涡黏或奇异尺度。模型验证时，壁面第一层网格和远场条件必须与所采用的模型版本一致。

### 5.4.2 Menter BSL 与 SST 模型

Menter BSL/SST 模型利用混合函数，在近壁区域采用 $k$–$\omega$ 形式，在远离壁面区域逐渐转向转换后的 $k$–$\varepsilon$ 形式。模型常数采用

$$
\phi=F_1\phi_1+(1-F_1)\phi_2 .
\tag{5.46}
$$

将 $k$–$\varepsilon$ 方程通过 $\varepsilon=\beta^*k\omega$ 转换为 $\omega$ 方程时，会出现交叉扩散项

$$
D_\omega=2(1-F_1)\rho\sigma_{\omega2}
\frac1\omega
\frac{\partial k}{\partial x_j}
\frac{\partial\omega}{\partial x_j}.
\tag{5.47}
$$

令 $d$ 为到壁面的距离，定义

$$
CD_{k\omega}=\max\left(
2\rho\sigma_{\omega2}\frac1\omega
\frac{\partial k}{\partial x_j}
\frac{\partial\omega}{\partial x_j},
\varepsilon_{CD}\right),
\tag{5.48}
$$

其中 $\varepsilon_{CD}$ 为防止除零的小正数。混合函数可写为

$$
F_1=\tanh(\arg_1^4),\qquad
F_2=\tanh(\arg_2^2),
\tag{5.49}
$$

$$
\arg_1=\min\left\{
\max\left[\frac{\sqrt k}{\beta^*\omega d},
\frac{500\nu}{d^2\omega}\right],
\frac{4\rho\sigma_{\omega2}k}{CD_{k\omega}d^2}
\right\},
\qquad
\arg_2=\max\left[\frac{2\sqrt k}{\beta^*\omega d},
\frac{500\nu}{d^2\omega}\right].
\tag{5.50}
$$

SST 还限制涡黏，以改善剪切层中的输运行为：

$$
\mu_t=\frac{\rho a_1k}
{\max(a_1\omega,\Omega F_2)},
\tag{5.51}
$$

其中 $\Omega$ 为涡量不变量，常取 $a_1=0.31$。常见参数包括 $\beta^*=0.09$、$\kappa=0.41$、$\sigma_{k1}=0.85$、$\sigma_{\omega1}=0.5$、$\beta_1=0.075$、$\sigma_{k2}=1.0$、$\sigma_{\omega2}=0.856$、$\beta_2=0.0828$。系数应按所采用的原始模型版本完整配套使用，不宜混用不同版本的定义。

本项目实现 Menter 等人 2003 年提出的 SST-2003 版本。该版本在式（5.51）中以应变率不变量 $S$ 代替涡量 $\Omega$，并在 $k$ 和 $\omega$ 两个方程中均采用限制生成项

$$
\widetilde P_k=\min(P_k,10\beta^*\rho k\omega).
$$

相应的生成项系数为 $\gamma_1=5/9$、$\gamma_2=0.44$。采用原始 SST 的涡量形式时，应明确标注为另一模型变体，不能与 SST-2003 的常数及限制器混称。

SST 兼具近壁积分和自由流稳健性，通常比标准 $k$–$\varepsilon$ 更适合逆压梯度和分离流，但它仍是基于涡黏假设的 RANS 模型，不能保证准确预测所有分离、转捩和强各向异性流动。

## 5.5 大涡模拟

### 5.5.1 RANS 的不足与 LES 的基本思想

RANS 对所有湍流尺度的影响进行统计建模，适合许多稳态工程计算，但在以下问题中可能不足：

- 需要瞬时涡结构、湍流噪声和非定常载荷；
- 大尺度分离、涡脱落、旋转流和强二次流；
- 混合层、射流和尾迹的发展对涡黏耗散过于敏感；
- 需要区分大尺度相干运动与随机小尺度脉动。

LES 直接求解网格能够分辨的大尺度运动，只对亚格子尺度建模。由于更小尺度的湍流更接近各向同性，亚格子模型通常比全谱 RANS 模型承担较少的物理描述任务，但 LES 必须采用足够细的三维网格和时间分辨率。

### 5.5.2 空间滤波与滤波方程

瞬时变量分解为滤波量与亚格子量：

$$
f=\overline f+f'.
\tag{5.52}
$$

物理空间滤波定义为

$$
\overline f(\mathbf x,t)=
\int_\Omega G(\mathbf x-\mathbf x';\Delta)
f(\mathbf x',t)\,\mathrm d\mathbf x',
\qquad
\int_\Omega G(\mathbf x;\Delta)\,\mathrm d\mathbf x=1.
\tag{5.53}
$$

盒式滤波器计算简单，常与有限体积网格结合；高斯滤波器具有较光滑的谱响应；截断 Fourier 滤波器在谱方法中常用。有限体积离散的网格本身也起到一定滤波作用。

与雷诺平均不同，空间滤波一般不满足“再次滤波不变”，也不一定与空间微分交换顺序。在交换误差可忽略或可作一致处理的条件下，不可压缩滤波方程为

$$
\frac{\partial\overline u_i}{\partial x_i}=0,
\tag{5.54}
$$

$$
\rho\left(
\frac{\partial\overline u_i}{\partial t}
+\frac{\partial\overline u_i\overline u_j}{\partial x_j}
\right)
=-\frac{\partial\overline p}{\partial x_i}
+\mu\frac{\partial^2\overline u_i}{\partial x_j\partial x_j}
-\frac{\partial\tau_{ij}^{sgs}}{\partial x_j},
\tag{5.55}
$$

其中亚格子应力为

$$
\tau_{ij}^{sgs}
=\rho\left(\overline{u_iu_j}-\overline u_i\overline u_j\right).
\tag{5.56}
$$

该未知项需要亚格子封闭模型。

### 5.5.3 Smagorinsky 与动态模型

亚格子涡黏模型采用

$$
\left(\tau_{ij}^{sgs}\right)^d
=\tau_{ij}^{sgs}-\frac13\tau_{kk}^{sgs}\delta_{ij}
=-2\mu_t\overline S_{ij},
\tag{5.57}
$$

其中 $\overline S_{ij}$ 为可解尺度应变率。Smagorinsky 模型给出

$$
\mu_t=\rho L_s^2|\overline S|,
\qquad
L_s=C_s\Delta,
\qquad
|\overline S|=\sqrt{2\overline S_{ij}\overline S_{ij}},
\qquad
\Delta\approx V^{1/3}.
\tag{5.58}
$$

理论估计 $C_s$ 约为0.18，工程中常取0.1–0.2。固定系数模型在近壁区域可能产生过大耗散，可使用

$$
L_s=\min(\kappa d,C_s\Delta)
\tag{5.59}
$$

或近壁阻尼函数。

动态 Smagorinsky 模型利用两种滤波尺度及 Germano 恒等式，根据局部流场动态确定系数。一种常用投影形式为

$$
C_s^2\approx
\frac{\langle L_{ij}M_{ij}\rangle}
{\langle M_{kl}M_{kl}\rangle},
\tag{5.60}
$$

其中 $L_{ij}$ 为两种尺度之间的可解析应力，$M_{ij}$ 为相应的模型张量，尖括号表示必要的局部或统计平均。动态过程还需处理分母过小和负系数等数值问题。

### 5.5.4 WALE 模型

WALE（Wall-Adapting Local Eddy-Viscosity）模型同时利用应变率和旋转信息，使涡黏在壁面附近自然衰减。定义速度梯度 $g_{ij}=\partial\overline u_i/\partial x_j$ 及

$$
s_{ij}^d=\frac12(g_{ik}g_{kj}+g_{jk}g_{ki})
-\frac13\delta_{ij}g_{kl}g_{lk}.
\tag{5.61}
$$

其涡黏可写为

$$
\mu_t=\rho L_s^2
\frac{(s_{ij}^ds_{ij}^d)^{3/2}}
{(\overline S_{ij}\overline S_{ij})^{5/2}
+(s_{ij}^ds_{ij}^d)^{5/4}},
\qquad
L_s=\min(C_w\Delta,\kappa d).
\tag{5.62}
$$

WALE 不需要额外的传统壁面函数即可获得合理的近壁衰减，适合壁湍流、旋转流和分离流 LES，但计算形式比 Smagorinsky 模型复杂。

LES 的结果验证不能只展示瞬时云图，还应检查网格滤波尺度、时间步、采样时间、平均量和二阶统计量。若网格和时间分辨率不足，再复杂的亚格子模型也无法弥补未解析的大尺度运动。

## 5.6 能量方程与湍流传热

### 5.6.1 焓形式与温度形式

低速、忽略压力功和黏性耗散等附加项时，焓形式的能量方程可写为

$$
\frac{\partial\rho h}{\partial t}
+\nabla\cdot(\rho\mathbf Uh)
=\nabla\cdot(k_f\nabla T)+q_T,
\tag{5.63}
$$

其中 $h$ 为比焓，$k_f$ 为导热系数，$q_T$ 为体积热源。若 $dh=c_p\,\mathrm dT$ 且 $c_p$ 为常数，可写为

$$
\rho c_p\left(\frac{\partial T}{\partial t}
+\mathbf U\cdot\nabla T\right)
=\nabla\cdot(k_f\nabla T)+q_T.
\tag{5.64}
$$

在常 $c_p$ 条件下，还可整理成第4章的通用温度输运形式：

$$
\frac{\partial\rho T}{\partial t}
+\nabla\cdot(\rho\mathbf UT)
=\nabla\cdot\left(\frac{k_f}{c_p}\nabla T\right)
+\frac{q_T}{c_p}.
\tag{5.65}
$$

选择总能量、总焓、内能、焓或温度作为求解变量，取决于可压缩性、状态方程和边界条件。变量形式改变时，源项和扩散系数也必须一致转换，不能只替换未知量符号。

### 5.6.2 有限体积离散与热边界条件

能量方程可直接采用通用对流—扩散离散：

$$
a_PT_P=\sum_{NB}a_{NB}T_{NB}+b_P.
\tag{5.66}
$$

界面质量流量应使用压力—速度耦合后满足离散连续性的流量，导热系数界面插值和热物性更新则按第4章处理。非线性辐射、相变或变物性源项可进行局部线性化。

设壁面热流 $q_w$ 以进入流体单元为正，典型边界包括：

- **给定壁温**：$T_w=T_s$，界面法向温度梯度由 $T_P$、$T_w$ 和中心到壁面的距离近似；
- **给定热流**：直接给定 $q_w=q_0$，加入能量方程右端；
- **对流传热边界**：

$$
q_w=h_{\mathrm{conv}}(T_P-T_w),
\tag{5.67}
$$

其中 $h_{\mathrm{conv}}$ 为壁面与流体之间的对流传热系数。采用湍流壁面函数时，热壁面函数必须与速度壁面处理和第一层网格位置相容。

### 5.6.3 雷诺平均能量方程与 SGDH

对能量方程进行雷诺平均后，会出现湍流热通量 $-\rho\overline{\mathbf U'h'}$。简单梯度扩散假设（SGDH）将其写为

$$
-\rho\overline{u_i'h'}=\lambda_t\frac{\partial T}{\partial x_i},
\qquad
\lambda_t=\frac{\mu_tc_p}{Pr_t},
\tag{5.68}
$$

其中 $Pr_t$ 为湍流普朗特数。于是平均能量方程中的有效导热系数为

$$
k_{\mathrm{eff}}=k_f+\lambda_t.
\tag{5.69}
$$

普通空气和水的许多工程计算常分别取 $Pr_t\approx0.9$ 和 $Pr_t\approx1.0$，但它并不是普适常数。SGDH 假定热通量方向与温度梯度方向相反，并通过标量涡黏建立动量与热量输运的比拟。

### 5.6.4 低、高及变 $Pr$ 数传热

分子普朗特数为

$$
Pr=\frac{\nu}{\alpha}
=\frac{\mu c_p}{k_f},
\qquad Pe=Pr\,Re.
\tag{5.70}
$$

液态金属具有很低的 $Pr$，例如钠和铅铋合金可分别处于0.01附近和0.02附近；超临界水和超临界 CO$_2$ 在拟临界区域的物性及 $Pr$ 会显著变化。此时动量与热量边界层不再具有简单相似性，固定 $Pr_t$ 的 SGDH 可能产生明显误差。

湍流热通量封闭可分为以下层级：

| 方法 | 基本处理 | 特点与限制 |
|---|---|---|
| 变 $Pr_t$ 的 SGDH | 用 $Pr_t=f(Pr,Pe)$ 等经验关系 | 实现简单，但依赖工质、几何和工况 |
| 代数热通量模型 AHFM | 直接建立 $\overline{u_i'T'}$ 的代数关系 | 可描述非梯度和各向异性输运，形式更复杂 |
| 微分热通量模型 DHFM | 求解 $\overline{u_i'T'}$ 的输运方程 | 物理内容丰富，方程和封闭项更多 |
| 四方程模型 | 求解 $k$、$\varepsilon$、$k_\theta$、$\varepsilon_\theta$ | 通过热时间尺度计算湍流热扩散，不依赖固定 $Pr_t$ |
| 机器学习模型 | 用数据拟合热通量或模型系数 | 依赖训练数据覆盖范围和物理约束 |

四方程模型引入温度脉动强度及其耗散率，常用定义为 $k_\theta=\tfrac12\overline{T'^2}$。其典型方程结构为

$$
\frac{\partial k_\theta}{\partial t}
+U_i\frac{\partial k_\theta}{\partial x_i}
=\frac{\partial}{\partial x_i}
\left[\left(\alpha+\frac{\alpha_t}{\sigma_{k\theta}}\right)
\frac{\partial k_\theta}{\partial x_i}\right]
+P_\theta-\varepsilon_\theta,
\tag{5.71}
$$

$$
\frac{\partial\varepsilon_\theta}{\partial t}
+U_i\frac{\partial\varepsilon_\theta}{\partial x_i}
=\frac{\partial}{\partial x_i}
\left[\left(\alpha+\frac{\alpha_t}{\sigma_{\varepsilon\theta}}\right)
\frac{\partial\varepsilon_\theta}{\partial x_i}\right]
+\frac{\varepsilon_\theta}{k_\theta}(C_{P1}P_\theta-C_{D1}\varepsilon_\theta)
+\frac{\varepsilon_\theta}{k}(C_{P2}P_k-C_{D2}\varepsilon).
\tag{5.72}
$$

湍流热扩散率通过局部热时间尺度计算：

$$
\overline{u_i'T'}=-\alpha_t\frac{\partial T}{\partial x_i},
\qquad
\alpha_t=C_\theta k\tau_l^\theta,
\qquad
\tau_l^\theta=\tau_l^\theta(k,\varepsilon,k_\theta,\varepsilon_\theta,Pr).
\tag{5.73}
$$

不同四方程和代数模型的时间尺度及常数定义并不完全相同，实现时应完整采用同一模型的方程组。低 $Pr$ 铅铋后台阶等算例表明，先进热通量模型可显著改善温度场和 Nu 数预测；但在分子扩散影响很小、强对流占主导的情况下，不同模型的差别可能减弱。

## 5.7 热膨胀与可压缩流动

### 5.7.1 热膨胀流体与低速模型

当密度随温度显著变化时，不能再把 $\rho$ 视为常数。热膨胀低速模型采用

$$
\rho=\rho(T)
\tag{5.74}
$$

并联合连续性、动量和能量方程。将连续性方程展开，有

$$
\nabla\cdot\mathbf U
=-\frac1\rho\left(\frac{\partial\rho}{\partial t}
+\mathbf U\cdot\nabla\rho\right)
=-\frac1\rho\frac{\mathrm D\rho}{\mathrm Dt}
=-\frac1\rho\frac{\mathrm d\rho}{\mathrm dT}
\frac{\mathrm DT}{\mathrm Dt}.
\tag{5.75}
$$

因此，速度散度不再为零，而由温度变化引起。低速求解时可将密度物质导数作为源项处理，用上一时间步或当前外迭代的温度场估计，再在能量方程更新后重新计算密度。整体仍可沿用不可压缩压力修正框架，但连续性残差必须包含密度变化。

### 5.7.2 Boussinesq 近似

若密度变化很小，可只在浮力项中保留密度变化，其余位置使用参考密度 $\rho_0$。线性热膨胀关系为

$$
\rho=\rho_0[1-\beta(T-T_0)],
\qquad
\beta=-\frac1{\rho_0}\left(\frac{\partial\rho}{\partial T}\right)_p,
\qquad
\beta|T-T_0|\ll1.
\tag{5.76}
$$

除密度外，其余物性也近似为常数。定义静水修正压力

$$
p_{rgh}=p-\rho_0\mathbf g\cdot\mathbf x,
\tag{5.77}
$$

Boussinesq 控制方程写为

$$
\nabla\cdot\mathbf U=0,
\tag{5.78}
$$

$$
\rho_0\frac{\mathrm D\mathbf U}{\mathrm Dt}
=-\nabla p_{rgh}-\rho_0\beta(T-T_0)\mathbf g
+\mu\nabla^2\mathbf U,
\tag{5.79}
$$

$$
\rho_0c_p\frac{\mathrm DT}{\mathrm Dt}
=\nabla\cdot(k_f\nabla T)+S_q.
\tag{5.80}
$$

与常密度方程相比，主要差别是增加了浮力源项。热流体通常因密度降低而产生与重力相反的浮力，从而驱动自然对流。Boussinesq 近似不适合大温差、强变物性或显著可压缩效应问题。

### 5.7.3 可压缩控制方程

可压缩流动必须联立质量、动量、能量和状态方程：

$$
\frac{\partial\rho}{\partial t}
+\nabla\cdot(\rho\mathbf U)=0,
\tag{5.81}
$$

$$
\frac{\partial\rho\mathbf U}{\partial t}
+\nabla\cdot(\rho\mathbf U\otimes\mathbf U)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf g,
\tag{5.82}
$$

$$
\frac{\partial\rho E}{\partial t}
+\nabla\cdot[\mathbf U(\rho E+p)]
=\nabla\cdot(k_f\nabla T)
+\nabla\cdot(\boldsymbol\tau\cdot\mathbf U)+S_q,
\tag{5.83}
$$

$$
\rho=\rho(p,T).
\tag{5.84}
$$

总比能 $E=e+\tfrac12|\mathbf U|^2$。能量方程使温度、密度和压力形成封闭耦合。对理想气体，状态方程为

$$
\rho=\frac{p}{RT}.
\tag{5.85}
$$

压力基方法可以推广到可压缩流动，但压力修正方程必须同时反映速度修正和密度修正。

### 5.7.4 低马赫数压力分解

低马赫数条件下，常将压力分解为空间均匀的热力学压力与机械压力：

$$
P(\mathbf x,t)=P_{th}(t)+p(\mathbf x,t),
\qquad |p|\ll P_{th}.
\tag{5.86}
$$

动量方程只保留机械压力梯度：

$$
\frac{\partial\rho\mathbf U}{\partial t}
+\nabla\cdot(\rho\mathbf U\otimes\mathbf U)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf g,
\tag{5.87}
$$

而状态方程使用热力学压力：

$$
\rho=\frac{P_{th}}{RT}.
\tag{5.88}
$$

于是，机械压力负责驱动运动，热力学压力和温度负责更新密度。连续性方程仍可整理为式（5.75）的形式，其中密度物质导数由状态方程及上一层次的温度、热力学压力计算；若 $P_{th}$ 随时间变化，还应同时包含 $(\partial\rho/\partial p)_T\,\mathrm DP_{th}/\mathrm Dt$ 的贡献。该模型适合低速热膨胀和浮力驱动流，不用于必须解析声学压力波或强压缩波的问题。

### 5.7.5 真实密度模型的压力耦合方程

真实密度模型直接采用 $\rho=\rho(p,T)$。令

$$
\rho^{n+1}=\rho^*+\rho',
\qquad
\mathbf U_f^{n+1}=\mathbf U_f^*+\mathbf U_f',
\qquad
p^{n+1}=p^*+p'.
\tag{5.89}
$$

在固定温度下线性化状态方程：

$$
\rho'\approx\left(\frac{\partial\rho}{\partial p}\right)_Tp'
=\Psi p',
\qquad
\Psi=\left(\frac{\partial\rho}{\partial p}\right)_T.
\tag{5.90}
$$

理想气体等温线性化给出 $\Psi=1/(RT)$。界面质量流量线性化为

$$
\dot m_f^{n+1}
\approx\dot m_f^*+\rho_f^*Q_f'+Q_f^*\rho_f',
\qquad
Q_f=\mathbf U_f\cdot\mathbf A_f,
\qquad
\dot m_f=\rho_fQ_f.
\tag{5.91}
$$

速度修正产生压力扩散项，密度修正产生压力存储项和压力对流项。令 $Q_f'=-C_f(\nabla p')_f\cdot\mathbf A_f$，则压力耦合方程可概括为

$$
\frac{\Psi V_P}{\Delta t}p_P'
+\sum_{f\in\mathcal F_P}
\left[Q_f^*\Psi_fp_f'
-\rho_f^*C_f(\nabla p')_f\cdot\mathbf A_f\right]
=-\frac{\rho^*-\rho^n}{\Delta t}V_P
-\sum_{f\in\mathcal F_P}\dot m_f^*.
\tag{5.92}
$$

$p_f'$ 可按迎风或其他相容格式插值。与不可压缩压力修正方程相比，式（5.92）额外包含：

- $\Psi V_P/\Delta t$ 对应的**可压缩存储项**；
- $Q_f^*\Psi_fp_f'$ 对应的**压力对流项**。

因此它不再是纯 Poisson 型方程，而是具有压力对流—扩散及存储结构的耦合方程。温度变化时，完整的密度增量还应包含 $(\partial\rho/\partial T)_pT'$；该项可在能量—压力外迭代中延迟处理。

### 5.7.6 可压缩能量方程与求解框架

可压缩流动必须求解能量方程。以理想气体、常 $c_p$ 等简化条件为例，焓形式可写成物质导数形式

$$
\rho c_p\frac{\mathrm DT}{\mathrm Dt}
=\frac{\mathrm Dp}{\mathrm Dt}
+\nabla\cdot(k_f\nabla T)
+\Phi+S_q,
\tag{5.93}
$$

其中 $\Phi$ 为黏性耗散。不同总能、总焓和温度形式的压力功及黏性项应保持等价。温度更新后，通过状态方程重新计算密度，并更新黏度、导热系数等热物性。

一种真实密度、压力基 SIMPLE 的时间步流程为：

1. 给定 $\mathbf U^n$、$\dot m_f^n$、$\rho^n$、$p^n$、$T^n$ 及当前外迭代初值。
2. 组装并求解动量方程，得到 $\mathbf U^*$。
3. 由状态方程计算 $\rho^*$，通过 Rhie–Chow 插值形成 $\dot m_f^*$。
4. 组装并求解含可压缩存储和压力对流项的 $p'$ 方程。
5. 修正速度、界面质量流量、密度和压力。
6. 求解能量方程，更新温度、密度和热物性。
7. 求解湍流或其他标量输运方程。
8. 检查所有残差和守恒平衡；未收敛则继续本时间步外迭代。
9. 收敛后更新旧时间层并推进 $t\leftarrow t+\Delta t$。

压力、温度、密度和界面流量必须采用相容的更新次序。只更新压力而忘记密度，或使用未修正的界面质量流量，都可能破坏能量与质量守恒。

## 本章小结

SIMPLE 通过“动量预测—压力修正—速度修正”的分离迭代实现压力—速度耦合。其核心近似是在速度修正中忽略邻点速度修正，并利用连续性方程把压力修正与质量残差联系起来。交错网格天然采用相邻压力构造界面压差；同位网格则需要 Rhie–Chow 动量插值抑制棋盘压力振荡。SIMPLEC 通过改进速度修正系数减轻邻点项省略的影响。

湍流计算可在同一分离框架中增加封闭方程。RANS 求解平均场并用模型表示全部湍流尺度，常用双方程模型包括 $k$–$\varepsilon$ 和 $k$–$\omega$；SST 结合两者优势并限制涡黏。LES 直接解析大尺度瞬态运动，只对亚格子应力建模，常用 Smagorinsky、动态 Smagorinsky 和 WALE 模型。

能量方程可作为另一个对流—扩散方程接入迭代。湍流热通量可用 SGDH、代数热通量模型或四方程模型封闭；在低、高或强变 $Pr$ 条件下，固定 $Pr_t$ 的可靠性需要特别检验。

热膨胀和可压缩性会改变连续性约束及压力方程。低速模型把密度物质导数作为源项，Boussinesq 近似只在浮力中保留小密度变化，低马赫模型区分热力学压力和机械压力，真实密度模型则通过 $\rho'=\Psi p'$ 在压力方程中引入可压缩存储和压力对流项。

## 习题

1. 说明不可压缩流动中压力没有独立时间演化方程时，为什么仍能通过压力修正获得速度场。
2. 证明式（5.4）的交替压力扰动在某些中心差分压力梯度中不可见，并解释 $1$–$\delta$ 压差为何能抑制棋盘振荡。
3. 从式（5.6）和式（5.7）推导速度修正关系，并指出 SIMPLE 近似省略了哪一项。
4. 将式（5.12）代入二维离散连续性方程，完整推导式（5.14）—（5.17）。
5. 讨论全 Neumann 压力修正系统的奇异性。为什么固定一个压力参考值不影响不可压缩速度？
6. 比较普通线性插值与 Rhie–Chow 插值的界面速度，说明两者对交替压力场的响应有何不同。
7. 推导 SIMPLEC 的速度修正系数，并分析分母变小时可能出现的数值问题。
8. 对不可压缩 Navier–Stokes 方程进行雷诺平均，说明雷诺应力的来源及其六个独立分量。
9. 从 $\mu_t=C_\mu\rho k^2/\varepsilon$ 和 $\varepsilon=\beta^*k\omega$ 推导 $k$–$\omega$ 模型中的 $\mu_t=\rho k/\omega$。
10. 比较标准 $k$–$\varepsilon$、$k$–$\omega$ 和 SST 模型的近壁处理、自由流敏感性及分离流适用性。
11. 推导滤波动量方程中的亚格子应力，并说明它与 RANS 雷诺应力在物理意义上的差别。
12. 分析 Smagorinsky 模型在层流区域或壁面附近可能产生过大涡黏的原因，并比较动态模型和 WALE 的改进思路。
13. 将式（5.63）转换为温度形式，列出所需假设，并说明变 $c_p$ 时应如何修改。
14. 为什么固定 $Pr_t$ 的 SGDH 在低 $Pr$ 液态金属中可能失效？比较 SGDH 与四方程模型需要求解的未知量。
15. 从 $\rho=\rho(T)$ 的连续性方程推导式（5.75），并说明低速模型如何把该关系纳入压力修正。
16. 推导 Boussinesq 动量方程中的浮力项，并解释为什么大温差下该近似失效。
17. 从 $\rho=\rho(p,T)$ 出发推导真实密度压力方程中的存储项和压力对流项。
18. 设计一个可压缩 SIMPLE 算例的收敛检查表，至少包括质量、动量、能量、状态方程和界面流量一致性。

## 延伸阅读

1. Patankar, S. V. *Numerical Heat Transfer and Fluid Flow*. Hemisphere, 1980.
2. Rhie, C. M., and Chow, W. L. Numerical study of the turbulent flow past an airfoil with trailing edge separation. *AIAA Journal*, 1983, 21: 1525–1532.
3. Versteeg, H. K., and Malalasekera, W. *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*. 2nd ed., Pearson, 2007.
4. Ferziger, J. H., and Perić, M. *Computational Methods for Fluid Dynamics*. 3rd ed., Springer, 2002.
5. Wilcox, D. C. *Turbulence Modeling for CFD*. 3rd ed., DCW Industries, 2006.
6. Menter, F. R. Two-equation eddy-viscosity turbulence models for engineering applications. *AIAA Journal*, 1994, 32: 1598–1605.
7. Pope, S. B. *Turbulent Flows*. Cambridge University Press, 2000.
8. Smagorinsky, J. General circulation experiments with the primitive equations. *Monthly Weather Review*, 1963, 91: 99–164.
9. Germano, M., Piomelli, U., Moin, P., and Cabot, W. H. A dynamic subgrid-scale eddy viscosity model. *Physics of Fluids A*, 1991, 3: 1760–1765.
