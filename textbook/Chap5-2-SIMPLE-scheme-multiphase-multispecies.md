# 第5章（续）多相流、多组分流动与化学反应

多相流和多组分流动都涉及多种物质，但二者采用的描述方式不同。多相流关注不同相之间的界面、相对运动和相间质量—动量—能量交换；多组分流动通常假设各组分在分子尺度上充分混合，重点描述浓度扩散、混合物物性和化学反应。两类问题都可以在压力基有限体积框架中求解，但需要增加体积分数、颗粒轨迹、组分质量分数和化学反应源项等方程。

本章依次介绍 VOF、代数滑移、Euler–Euler、Euler–Lagrange、多组分输运及化学反应模型。公式编号承接第5章前半部分。

## 5.8 多相流与 VOF 模型

### 5.8.1 多相流型与模拟方法

常见气液流型包括：

- **泡状流**：小气泡分散在连续液体中；
- **液滴流**：小液滴分散在连续气体中；
- **弹状流**：气泡聚集成较大的气团，气液呈块状交替分布；
- **环状流**：液体在壁面形成液膜，气体主要在中部流动；
- **分层流和自由表面流**：互不相溶的流体由较清晰的界面分隔。

流型取决于流量、物性、管径、重力和历史发展过程，工程中常用流型图作初步判断。相应的数值方法包括经验关联、拉格朗日颗粒追踪、代数滑移、多流体、离散元及全耦合界面解析方法。选择模型时，首先要判断相界面是否需要且能够被网格解析。

### 5.8.2 VOF 的适用范围与体积分数

流体体积函数（Volume of Fluid，VOF）方法适用于具有清晰界面的互不相溶流体，例如自由表面、溃坝、射流破裂和大气泡运动。当气泡或液滴远小于控制体、数量众多且界面无法逐个解析时，VOF 不再合适。

定义第 $k$ 相的指示函数 $\chi_k$：该相占据的位置为1，其余位置为0。单元 $V_P$ 中的体积分数为

$$
C_k(V_P)=\frac1{V_P}\int_{V_P}\chi_k(\mathbf x)\,\mathrm dV,
\qquad
\sum_{k=1}^{N_{\mathrm{phase}}}C_k=1.
\tag{5.94}
$$

因此，

$$
C_k=\begin{cases}
0, & \text{单元内不含第 }k\text{ 相},\\
1, & \text{单元被第 }k\text{ 相充满},\\
0<C_k<1, & \text{单元包含相界面}.
\end{cases}
\tag{5.95}
$$

每个相的体积分数满足

$$
\frac{\partial C_k}{\partial t}
+\nabla\cdot(C_k\mathbf U)=S_{C_k},
\tag{5.96}
$$

其中 $S_{C_k}$ 可表示相间传质等体积分数源项。无源且速度散度为零时，式（5.96）也可写为 $\partial C_k/\partial t+\mathbf U\cdot\nabla C_k=0$。

### 5.8.3 混合物控制方程与表面张力

VOF 假设同一单元内各相共享速度和压力，只求解一组混合物动量方程：

$$
\frac{\partial\rho_m\mathbf U}{\partial t}
+\nabla\cdot(\rho_m\mathbf U\otimes\mathbf U)
=-\nabla p+\nabla\cdot\boldsymbol\tau
+\mathbf F_T+\rho_m\mathbf g,
\tag{5.97}
$$

其中常用体积分数加权物性为

$$
\rho_m=\sum_kC_k\rho_k,
\qquad
\mu_m=\sum_kC_k\mu_k.
\tag{5.98}
$$

不可压缩、无相间净质量生成时，混合物连续性方程为 $\nabla\cdot\mathbf U=0$。能量方程可写为

$$
\frac{\partial\rho_mh_m}{\partial t}
+\nabla\cdot(\rho_m\mathbf Uh_m)
=\nabla\cdot(\lambda_m\nabla T)+q,
\tag{5.99}
$$

并满足体积平均关系 $\rho_mh_m=\sum_kC_k\rho_kh_k$。因此，比焓和比热容若按单位质量定义，应采用质量加权，而不能在所有情况下都直接使用体积分数算术平均。导热系数等输运物性的混合规则还取决于界面形态和热流方向。

连续表面力（CSF）模型把表面张力写成界面附近的体积力：

$$
\mathbf F_T=\sigma\kappa\mathbf n\delta_s,
\qquad
\mathbf n=\frac{\nabla C}{|\nabla C|},
\qquad
\kappa=-\nabla\cdot\mathbf n,
\qquad
\delta_s\approx|\nabla C|.
\tag{5.100}
$$

于是 $\mathbf F_T\approx\sigma\kappa\nabla C$。曲率估计、法向计算和界面涂抹会影响寄生流等数值误差，通常需要适当的光滑和界面锐化。

### 5.8.4 体积分数的有限体积方程

对固定单元离散式（5.96），得到

$$
\frac{C_{k,P}^{n+1}-C_{k,P}^n}{\Delta t}V_P
+\sum_{f\in\mathcal F_P}C_{k,f}Q_f
=V_PS_{C_k,P},
\qquad
Q_f=\mathbf U_f\cdot\mathbf A_f.
\tag{5.101}
$$

关键问题是如何构造界面体积分数 $C_{k,f}$。普通一阶迎风虽然有界，但会把界面迅速抹平；无约束的高阶格式又可能产生 $C_k<0$ 或 $C_k>1$。VOF 通量格式必须同时考虑守恒、有界性、低数值扩散和界面锐度。

时间步还受到体积分数库朗数限制：

$$
Co_f=\frac{|Q_f|\Delta t}{V_P}.
\tag{5.102}
$$

部分显式界面格式以 $Co<1$ 为基本要求，实际计算常采用更小的时间步以控制界面误差。

### 5.8.5 代数界面捕捉与 CICSAM

代数方法直接使用单元体积分数构造界面通量。除 Downwind 和 STACS 等格式外，CICSAM 是任意网格上常用的压缩界面捕捉方法。

设 $U$ 为远上游单元、$D$ 为直接上游供体单元、$A$ 为下游受体单元。归一化变量定义为

$$
\widetilde C=\frac{C-C_U}{C_A-C_U},
\qquad
\widetilde C_U=0,\qquad \widetilde C_A=1.
\tag{5.103}
$$

归一化界面值 $\widetilde C_f$ 必须满足有界性约束；当 $0<\widetilde C_D<1$ 时，应使界面值处于供体值与1之间。HYPER-C 具有高压缩性：

$$
\widetilde C_{f,\mathrm{HYPER-C}}=
\begin{cases}
\widetilde C_D, & \widetilde C_D<0\ \text{或}\ \widetilde C_D>1,\\[2pt]
\min\left(1,\dfrac{\widetilde C_D}{Co_f}\right), & 0<\widetilde C_D<1.
\end{cases}
\tag{5.104}
$$

但界面与流动方向接近相切时，过强压缩可能引起界面失真。CICSAM 在 HYPER-C 与较光滑的 ULTIMATE–QUICKEST（UQ）格式之间切换：

$$
\widetilde C_{f,\mathrm{CICSAM}}
=\gamma_f\widetilde C_{f,\mathrm{HYPER-C}}
+(1-\gamma_f)\widetilde C_{f,\mathrm{UQ}},
\tag{5.105}
$$

$$
\gamma_f=\min\left(k_\gamma\frac{1+\cos2\theta_f}{2},1\right),
\qquad
\cos\theta_f=\frac{|\mathbf d\cdot\mathbf n|}{|\mathbf d||\mathbf n|}.
\tag{5.106}
$$

其中 $\mathbf d$ 为界面两侧单元中心连线，$\mathbf n$ 为供体单元体积分数梯度方向，$k_\gamma$ 常取1；$k_\gamma=0$ 时退化为 UQ。得到 $\widetilde C_f$ 后，再反变换为 $C_f$ 并计算体积分数通量。

### 5.8.6 几何界面重构：SLIC 与 PLIC

几何方法先在混合单元内重构界面，再计算一个时间步内穿过界面的相体积。

- **SLIC** 用平行于某一坐标方向的直线或平面近似界面，构造简单，但方向性较强。
- **PLIC** 用任意方向的直线或平面近似界面，精度和界面形状表达能力更好。

PLIC 在单元内构造

$$
\mathbf n\cdot\mathbf x=d,
\qquad
\mathbf n\approx-\nabla C,
\tag{5.107}
$$

并选择距离 $d$，使平面一侧的相体积与已知 $C_kV_P$ 相等。可以枚举平面经过各顶点时的 $d_{\min}$、$d_{\max}$，再用割线法或求根算法确定 $d$。

随后，根据界面速度构造时间步内穿过某一单元面的输运多面体，计算其中属于第 $k$ 相的几何体积。该体积就是界面相通量。PLIC 通常比代数方法更锐利，在高曲率和较大时间步下更有利于保持相体积，但三维任意多面体上的几何求交、体积计算和并行实现明显更复杂。

### 5.8.7 VOF 求解流程与检查

一个典型的 VOF 时间步为：

1. 由当前 $C_k$ 重构界面并计算体积分数界面通量；
2. 求解体积分数方程，检查 $\sum_kC_k=1$ 和 $0\leq C_k\leq1$；
3. 更新混合物密度、黏度及其他物性；
4. 计算界面法向、曲率和表面张力；
5. 用 SIMPLE、PIMPLE 或其他压力基方法求解混合物速度和压力；
6. 求解能量及其他耦合方程；
7. 检查总体积、总质量、界面位置和目标物理量，随后推进时间步。

对越界值进行简单裁剪只能作为很小数值误差的处理手段，不能代替有界通量格式，也不能用来掩盖显著的非守恒。溃坝等验证算例应同时比较自由面形状、测点压力、总液体体积和时间步敏感性。

## 5.9 代数滑移与 Euler–Euler 多流体模型

### 5.9.1 VOF 的局限与其他模型

VOF 需要网格解析界面，三维瞬态计算成本较高，不适合大量亚网格气泡、液滴或颗粒。对于高度分散体系，可以采用：

- **代数滑移模型（ASM）**：求解混合物方程和各相体积分数，用代数关系估计相间滑移；
- **Euler–Euler 多流体模型**：分别求解各相连续性、动量和能量方程，通过相间交换项耦合；
- **Euler–Lagrange 模型**：连续相在欧拉网格上求解，离散相逐个或按颗粒包追踪。

### 5.9.2 代数滑移模型

定义混合密度、质量平均速度和漂移速度：

$$
\rho_m=\sum_kC_k\rho_k,
\qquad
\mathbf U_m=\frac{\sum_kC_k\rho_k\mathbf U_k}{\rho_m},
\qquad
\mathbf U_{dr,k}=\mathbf U_k-\mathbf U_m.
\tag{5.108}
$$

混合物连续性和动量方程为

$$
\frac{\partial\rho_m}{\partial t}
+\nabla\cdot(\rho_m\mathbf U_m)=0,
\tag{5.109}
$$

$$
\begin{aligned}
\frac{\partial\rho_m\mathbf U_m}{\partial t}
+\nabla\cdot(\rho_m\mathbf U_m\otimes\mathbf U_m)
={}&-\nabla p+\nabla\cdot\boldsymbol\tau_m+\rho_m\mathbf g\\
&-\nabla\cdot\left(\sum_kC_k\rho_k
\mathbf U_{dr,k}\otimes\mathbf U_{dr,k}\right).
\end{aligned}
\tag{5.110}
$$

次相 $p$ 的体积分数方程可写为

$$
\frac{\partial C_p\rho_p}{\partial t}
+\nabla\cdot(C_p\rho_p\mathbf U_m)
=-\nabla\cdot(C_p\rho_p\mathbf U_{dr,p})
+\sum_q(\dot m_{qp}-\dot m_{pq}),
\tag{5.111}
$$

其中 $\dot m_{qp}$ 表示由 $q$ 相转入 $p$ 相的质量率。ASM 适合各相运动方向基本一致、滑移较小的流动，并通过经验关系避免为每一相单独求解完整动量方程。

### 5.9.3 滑移速度与漂移通量

定义次相 $p$ 相对主相 $q$ 的滑移速度

$$
\mathbf U_{pq}=\mathbf U_p-\mathbf U_q.
\tag{5.112}
$$

漂移速度与滑移速度的关系为

$$
\mathbf U_{dr,p}=\mathbf U_{pq}
-\sum_k\frac{C_k\rho_k}{\rho_m}\mathbf U_{kq}.
\tag{5.113}
$$

一种局部平衡漂移通量模型给出

$$
\mathbf U_{pq}
=\frac{\tau_p}{f_{drag}}\frac{\rho_p-\rho_m}{\rho_p}\mathbf a,
\qquad
\mathbf a=\mathbf g-\frac{\mathrm D\mathbf U_m}{\mathrm Dt},
\qquad
\tau_p=\frac{\rho_pd_p^2}{18\mu_q}.
\tag{5.114}
$$

Schiller–Naumann 阻力修正为

$$
f_{drag}=\frac{C_DRe_p}{24}
=\begin{cases}
1+0.15Re_p^{0.687}, & Re_p\lesssim1000,\\
0.0183Re_p, & Re_p\gtrsim1000.
\end{cases}
\tag{5.115}
$$

其中相对雷诺数按颗粒直径、滑移速度和主相物性计算。湍流条件下还需考虑离散相的湍流扩散。漂移通量模型形式简单，但其准确性高度依赖局部平衡和滑移关系是否成立。

### 5.9.4 Euler–Euler 两流体方程

Euler–Euler 模型为每一相分别建立守恒方程，并假设各相共享压力。以 $q$ 相为例：

$$
\frac{\partial C_q\rho_q}{\partial t}
+\nabla\cdot(C_q\rho_q\mathbf U_q)
=\dot m_{pq}-\dot m_{qp}+S_q,
\tag{5.116}
$$

$$
\begin{aligned}
&\frac{\partial C_q\rho_q\mathbf U_q}{\partial t}
+\nabla\cdot(C_q\rho_q\mathbf U_q\otimes\mathbf U_q)\\
&\quad=-C_q\nabla p+\nabla\cdot\boldsymbol\tau_q
+C_q\rho_q\mathbf g+\mathbf R_{pq}
+\dot m_{pq}\mathbf U_{pq}-\dot m_{qp}\mathbf U_{qp}
+\mathbf F_{lift,q}+\mathbf F_{vm,q}+\mathbf F_q .
\end{aligned}
\tag{5.117}
$$

各相体积分数满足

$$
C_p+C_q=1.
\tag{5.118}
$$

传质所携带的动量采用供体相速度。例如 $p$ 相向 $q$ 相传质时，进入 $q$ 相的动量按 $\mathbf U_p$ 计算。相间作用力满足作用—反作用关系

$$
\mathbf R_{pq}=-\mathbf R_{qp},
\qquad
\mathbf F_{*,p}=-\mathbf F_{*,q}.
\tag{5.119}
$$

### 5.9.5 相间阻力、升力与虚拟质量力

相间阻力常写成

$$
\mathbf R_{pq}=K_{pq}(\mathbf U_p-\mathbf U_q),
\qquad
K_{pq}=K_{qp}
=\frac{C_pC_q\rho_pf_{drag}}{\tau_p},
\tag{5.120}
$$

$$
Re_p=\frac{\rho_q|\mathbf U_p-\mathbf U_q|d_p}{\mu_q}.
\tag{5.121}
$$

颗粒或气泡上的升力可示意写为

$$
\mathbf F_{lift,p}
=-C_L C_p\rho_q(\mathbf U_q-\mathbf U_p)
\times(\nabla\times\mathbf U_q),
\tag{5.122}
$$

其中 $C_L$ 依模型和流动条件确定。升力与主相速度梯度及滑移有关，对很小或紧密堆积颗粒不宜直接套用稀相关系。

虚拟质量力表示次相相对主相加速时带动周围流体所产生的惯性效应：

$$
\mathbf F_{vm,p}
=C_{VM}C_p\rho_q
\left(\frac{\mathrm D_q\mathbf U_q}{\mathrm Dt}
-\frac{\mathrm D_p\mathbf U_p}{\mathrm Dt}\right),
\tag{5.123}
$$

其中 $C_{VM}$ 常取约0.5，$\mathrm D_q/\mathrm Dt=\partial/\partial t+\mathbf U_q\cdot\nabla$。当次相密度远小于主相密度时，虚拟质量效应尤其重要。上述封闭项通常非线性，并依赖粒径和流型，是多流体模型收敛较慢和不确定性的主要来源。

### 5.9.6 多相模型的选择

| 模型 | 典型前提 | 主要输出 | 典型问题 |
|---|---|---|---|
| VOF | 界面可由网格解析，各相共享速度 | 界面位置和形状 | 自由表面、分层流、微通道液滴 |
| ASM | 分散较均匀，滑移较小且可代数估计 | 混合物场、体积分数和漂移速度 | 液固悬浮、低到中等浓度鼓泡 |
| Euler–Euler | 各相可具有不同速度，需描述相间作用 | 各相体积分数、速度和温度 | 鼓泡塔、流化床、沉降 |
| Euler–Lagrange | 离散相较稀疏，颗粒间作用可忽略或另建模 | 颗粒轨迹及相间源项 | 喷雾、粉尘、蒸发液滴 |

课件中的体积分数阈值可作为初步经验：Euler–Lagrange 常用于低于约10%的离散相，ASM 常用于一定范围内的中等浓度分散体系，Euler–Euler 可处理更高浓度和更复杂的相间运动。这些数值不是普适判据，粒径、Stokes 数、碰撞和聚并破碎同样重要。

## 5.10 Euler–Lagrange 颗粒追踪

### 5.10.1 颗粒运动方程

Euler–Lagrange 方法适合固体颗粒、液滴或气泡分散在连续流体中的问题，可计算轨迹以及相间质量、动量和热量交换。单个球形颗粒满足

$$
\frac{\mathrm d\mathbf X_p}{\mathrm dt}=\mathbf V_p,
\tag{5.124}
$$

$$
m_p\frac{\mathrm d\mathbf V_p}{\mathrm dt}
=\frac12C_D\rho_g|\mathbf V_g-\mathbf V_p|
(\mathbf V_g-\mathbf V_p)A_p
+V_p(\rho_p-\rho_g)\mathbf g+m_p\mathbf F_{add},
\tag{5.125}
$$

其中

$$
m_p=\rho_pV_p=\rho_p\frac{\pi d_p^3}{6},
\qquad
A_p=\frac{\pi d_p^2}{4}.
\tag{5.126}
$$

将阻力项写成粒子弛豫形式，可得

$$
\frac{\mathrm d\mathbf V_p}{\mathrm dt}
=F_D(\mathbf V_g-\mathbf V_p)
+\frac{\rho_p-\rho_g}{\rho_p}\mathbf g+\mathbf F_{add},
\qquad
F_D=\frac{3\mu_gC_DRe_p}{4\rho_pd_p^2},
\tag{5.127}
$$

$$
Re_p=\frac{\rho_gd_p|\mathbf V_p-\mathbf V_g|}{\mu_g},
\qquad
C_D=\begin{cases}
\dfrac{24}{Re_p}\left(1+\dfrac16Re_p^{2/3}\right), & Re_p\lesssim1000,\\[6pt]
0.424, & Re_p\gtrsim1000.
\end{cases}
\tag{5.128}
$$

附加力可包括升力、虚拟质量力、压力梯度力、Brown 力和热泳力等，应根据粒径和流动条件选取。

### 5.10.2 颗粒传热与传质

令 $\dot m_{evap}>0$ 为蒸发质量率，则颗粒质量满足

$$
\frac{\mathrm dm_p}{\mathrm dt}=-\dot m_{evap},
\qquad
\dot m_{evap}=h_mA_p\rho_\infty\ln(1+B_M),
\tag{5.129}
$$

其中 $B_M$ 为传质驱动势。若把蒸发潜热 $L$ 单独计入能量平衡，颗粒温度满足

$$
m_pc_p\frac{\mathrm dT_p}{\mathrm dt}
=h_qA_p(T_\infty-T_p)-\dot m_{evap}L.
\tag{5.130}
$$

传热传质系数常用 Frössling 或 Ranz–Marshall 关联：

$$
Sh=\frac{h_md_p}{\mathcal D}
=2.0+0.6Re_p^{1/2}Sc^{1/3},
\qquad
Nu=\frac{h_qd_p}{\lambda}
=2.0+0.6Re_p^{1/2}Pr^{1/3}.
\tag{5.131}
$$

蒸发时 $d_p$ 随质量损失减小，组分、温度和直径必须同步更新。若发生沸腾、表面燃烧或凝结，应替换为相应的质量和能量模型。

### 5.10.3 单向与双向耦合

**单向耦合**只考虑连续相对颗粒的阻力、热传导和湍流作用，颗粒不改变连续相，适合极稀薄且负载很小的体系。

**双向耦合**还将颗粒产生的质量、动量和能量源项反馈到连续相：

$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf V_f)=S_\rho,
\tag{5.132}
$$

$$
\frac{\partial\rho y_k}{\partial t}
+\nabla\cdot(\rho\mathbf V_fy_k)
=-\nabla\cdot\mathbf j_k+S_{y_k},
\tag{5.133}
$$

$$
\frac{\partial\rho\mathbf V_f}{\partial t}
+\nabla\cdot(\rho\mathbf V_f\otimes\mathbf V_f)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf g+\mathbf S_V,
\tag{5.134}
$$

$$
\frac{\partial\rho h}{\partial t}+\nabla\cdot(\rho\mathbf V_fh)
=\nabla\cdot(\lambda\nabla T)+S_q.
\tag{5.135}
$$

在单元 $V_c$ 中，颗粒源项由该单元内所有颗粒或颗粒包累加并除以 $V_c$。例如，蒸发使气相质量增加，颗粒对气体的阻力与气体对颗粒的阻力大小相等、方向相反，颗粒与气体之间的换热量也应成对出现。若使用统计颗粒包，还必须乘以每个包代表的实际颗粒数，避免把颗粒包当成单个颗粒而低估耦合强度。

### 5.10.4 湍流分散与湍流调制

湍流对颗粒的影响包括两方面：

- **湍流分散效应**：瞬时速度脉动改变颗粒轨迹、沉积、聚并和破碎；
- **湍流调制效应**：颗粒反过来削弱或增强连续相湍流，需要在湍流方程中增加相应源项。

在 RANS 流场中，离散随机行走模型用随机速度近似湍流脉动：

$$
u'=v'=w'=\lambda\sqrt{\frac{2k}{3}},
\qquad
\tau_e\approx0.3\frac{k}{\varepsilon},
\tag{5.136}
$$

其中 $\lambda$ 为标准正态随机数，各速度分量应独立抽样。颗粒在一个涡寿命内跟随当前脉动速度，寿命结束后重新抽样。LES 中则优先直接使用可解尺度瞬时速度，避免重复加入不适当的随机脉动。

### 5.10.5 颗粒追踪计算流程

稳态双向耦合通常采用以下外迭代：

1. 求解连续相流场、温度和组分；
2. 在连续相场上积分颗粒轨迹、直径和温度；
3. 计算各单元的质量、动量和能量交换源项；
4. 更新湍流分散或调制模型；
5. 将源项反馈连续相并重新求解；
6. 重复至连续相、颗粒轨迹和源项均收敛。

非稳态问题则需要在每个时间步内处理颗粒穿格、边界碰撞、蒸发和源项累加。颗粒时间步与流体时间步可以不同，但必须保证轨迹积分和交换量的时间一致性。

## 5.11 多组分气体流动

### 5.11.1 多组分问题的类型

多组分流动假设组分之间相互扩散并充分混合，浓度梯度驱动的传质通常用 Fick 定律描述。根据组分对流场的影响，可分为：

1. **稀组分被动输运**：主体流动几乎不受组分影响，如低浓度污染物扩散；
2. **不可压缩或低马赫变密度混合物**：各组分影响密度、黏度等物性，进而影响主体流动；
3. **可压缩多组分气体**：状态方程、能量和扩散通量均与组分耦合，还需要考虑组分焓输运和化学反应。

多相流中的“相”可以具有清晰界面和不同速度，而多组分模型中的各组分通常共享同一个混合物速度，并以扩散速度相对混合物运动。

### 5.11.2 组分守恒与扩散通量

设 $y_k$ 为第 $k$ 种组分质量分数，则

$$
\sum_{k=1}^Ny_k=1,
\qquad
\frac{\partial\rho y_k}{\partial t}
+\nabla\cdot(\rho\mathbf Vy_k)
=-\nabla\cdot\mathbf j_k+\dot\omega_k,
\tag{5.137}
$$

其中 $\dot\omega_k$ 为单位时间、单位体积的组分质量净生成率。简单 Fick 扩散写为

$$
\mathbf j_k=-\rho D_{k,mix}\nabla y_k,
\qquad
\sum_k\mathbf j_k=0.
\tag{5.138}
$$

所有组分方程相加应回到总质量方程，因此 $\sum_k\dot\omega_k=0$。实际多组分扩散可能需要 Maxwell–Stefan 方程，简单 Fick 形式还需修正最后一个组分的扩散通量，以保证 $\sum_k\mathbf j_k=0$。数值上可只独立求解 $N-1$ 个组分并计算最后一个组分，或在归一化时保持总质量一致，不能依靠反复裁剪负浓度来获得形式上合理的结果。

### 5.11.3 不可压缩变密度多组分流动

混合物连续性和动量方程为

$$
\frac{\partial\rho}{\partial t}
+\nabla\cdot(\rho\mathbf V)=0,
\tag{5.139}
$$

$$
\rho\frac{\mathrm DV_i}{\mathrm Dt}
=-\frac{\partial p}{\partial x_i}
+\frac{\partial\tau_{ij}}{\partial x_j}
+(\rho-\rho_0)g_i.
\tag{5.140}
$$

密度和黏度由组分决定，例如

$$
\rho=\rho(y_1,\ldots,y_N,T),
\qquad
\mu=\mu(y_1,\ldots,y_N,T),
\qquad
D_{k,mix}=D_{k,mix}(y_1,\ldots,y_N,T,p).
\tag{5.141}
$$

低马赫条件下，可沿用热膨胀流动的散度源项：

$$
\nabla\cdot\mathbf V=-\frac1\rho\frac{\mathrm D\rho}{\mathrm Dt}.
\tag{5.142}
$$

求解时，先用预估压力求速度，再把组分—密度变化和试探速度质量残差共同引入压力修正方程。SIMPLE 的主要流程不变，但密度和扩散系数必须随组分更新。

### 5.11.4 可压缩多组分控制方程

可压缩多组分气体满足总质量、动量、组分和总焓方程。定义

$$
H=\sum_ky_kh_k+\frac12|\mathbf V|^2,
\qquad
h_k=h_{f,k}^0+\int_{T_0}^Tc_{p,k}(T')\,\mathrm dT',
\tag{5.143}
$$

其中 $h_{f,k}^0$ 为参考生成焓。总焓方程可写为

$$
\frac{\partial\rho H}{\partial t}
+\nabla\cdot(\rho\mathbf VH)
=\frac{\partial p}{\partial t}
+\nabla\cdot\left(\lambda\nabla T-\sum_kh_k\mathbf j_k\right)
+\nabla\cdot(\boldsymbol\tau\cdot\mathbf V)
-\sum_kh_{f,k}^0\dot\omega_k.
\tag{5.144}
$$

右端第二项包含热传导和组分扩散携带的焓，最后一项是化学反应引起的能量变化。若使用相对焓 $\widetilde h_k=h_k-h_{f,k}^0$ 构造 $\widetilde H=\sum_ky_k\widetilde h_k+|\mathbf V|^2/2$，生成焓基准移入化学反应源项，能量方程的形式和源项解释会更清晰。

理想气体混合物满足

$$
p=\rho R_{mix}T,
\qquad
R_{mix}=R_u\sum_k\frac{y_k}{M_k}
=\frac{R_u}{M_{mix}},
\tag{5.145}
$$

其中 $M_k$ 为摩尔质量，$R_u$ 为通用气体常数。混合物摩尔质量满足 $M_{mix}^{-1}=\sum_ky_k/M_k$。

### 5.11.5 混合物物性

Wilke 规则可估算混合黏度：

$$
\mu=\sum_k\frac{x_k\mu_k}{\Phi_k},
\qquad
\Phi_k=\sum_lx_l\phi_{kl},
\tag{5.146}
$$

$$
\phi_{kl}=\frac{\left[1+(\mu_k/\mu_l)^{1/2}(M_l/M_k)^{1/4}\right]^2}
{\sqrt8(1+M_k/M_l)^{1/2}},
\tag{5.147}
$$

其中 $x_k$ 为摩尔分数。一种简单的导热系数混合规则为

$$
\lambda=\frac12\left[\sum_kx_k\lambda_k
+\left(\sum_k\frac{x_k}{\lambda_k}\right)^{-1}\right].
\tag{5.148}
$$

由 Chapman–Enskog 理论得到二元扩散系数后，多组分有效扩散系数可近似为

$$
D_{k,mix}=\frac{1-x_k}
{\displaystyle\sum_{l\ne k}\frac{x_l}{D_{kl}}}.
\tag{5.149}
$$

这些混合规则具有各自的适用条件；高压、强非理想或极端物性体系应采用更完整的热物性模型。

### 5.11.6 多组分求解流程

一种 $N$ 组分可压缩分离式求解过程为：

1. 求解可压缩 Navier–Stokes 方程，得到速度和压力；
2. 在密度暂时固定的条件下求解组分输运方程，并保证质量分数归一和扩散通量相容；
3. 更新 $\rho$、$R_{mix}$、$\mu_{eff}$、$\lambda_{eff}$ 和 $D_{k,mix}$；
4. 求解包含组分焓扩散和化学源项的能量方程；
5. 求解湍流方程并更新有效黏度和扩散系数；
6. 重复至速度、压力、组分、温度和物性收敛。

不可压缩变密度驱动方腔、丙烷—空气同轴射流等验证问题表明，除速度和浓度云图外，还应比较截面组分剖面、轴向速度、混合长度和整体组分质量守恒。

### 5.11.7 扩散通量修正的离散一致性

为了保证所有组分方程之和严格回到总质量方程，可以给未修正的 Fick 通量增加一个共同的一致性修正：

$$
\mathbf j_k^0=-\rho D_{k,mix}\nabla y_k,
\qquad
\mathbf j_k=\mathbf j_k^0+y_k\sum_{l=1}^{N}\rho D_{l,mix}\nabla y_l .
\tag{5.149a}
$$

由于 $\sum_k y_k=1$，立即可得

$$
\sum_k\mathbf j_k
=-\sum_k\rho D_{k,mix}\nabla y_k
+\left(\sum_ky_k\right)\sum_l\rho D_{l,mix}\nabla y_l=0.
\tag{5.149b}
$$

该修正是混合物平均近似的一种数值封闭，并不等价于完整 Maxwell–Stefan 方程，也不能产生后者可能出现的逆向扩散、渗透扩散或扩散障壁。若只求解 $N-1$ 个方程，则最后一个组分的通量必须由 $-\sum_{k<N}\mathbf j_k$ 计算，不能同时又随意指定一个独立的末组分扩散通量。

本项目当前实现求解前 $N-1$ 个修正组分方程，并由 $y_N=1-\sum_{k<N}y_k$ 直接计算末组分。末组分不再独立求解，其扩散通量由 $-\sum_{k<N}\mathbf j_k$ 封闭；其扩散系数仍参与其他组分的修正式（5.149a）。修正式的第二项在 SIMPLE 外迭代中显式延迟。对每个内部面只构造一次带方向的修正通量 $J^c_{k,f}$，在拥有者单元右端加入 $-J^c_{k,f}$，在邻居单元右端加入 $+J^c_{k,f}$，因此内部面通量成对抵消。修正通量使用与未修正扩散项相同的法向梯度和非正交分解，湍流扩散系数也包含在同一封闭中。边界处采用与标量方程一致的 Dirichlet 或零法向扩散通量处理，末组分的边界值同样由质量分数约束计算。

这种约束形式避免了同时求解全部 $N$ 个显式修正方程时可能出现的非物理组成定点。只要速度质量通量满足离散连续性，前 $N-1$ 个方程和末组分质量约束相加就与总质量方程一致。外迭代收敛判据仍必须同时包括所有组分的更新量、$\max_P|\sum_k y_{k,P}-1|$ 和质量分数上下界。不能通过每一步大幅裁剪负值来伪造有界性；高阶对流产生的显著越界应通过限制器、网格或格式选择处理。

### 5.11.8 材料模块与多相/多组分软件框架

材料模型与输运方程应分离。材料层只根据局部状态计算物性，不拥有速度、压力、体积分数或组分方程；输运层则不应包含某一套具体气体混合规则。本项目采用以下统一 SI 接口：

| 接口数据 | 含义与单位 |
|---|---|
| `MaterialState::temperature` | 温度 $T$，K |
| `MaterialState::pressure` | 热力学绝对压力，Pa；不同于不可压缩 SIMPLE 的表压或压力修正 |
| `MaterialState::massFractions` | 组分质量分数 $y_k$ |
| `MaterialState::volumeFractions` | 相体积分数 $C_k$ |
| `SpeciesProperties::molarMass` | 摩尔质量，kg/mol |
| `MaterialProperties::density` | 密度，kg/m³ |
| `MaterialProperties::viscosity` | 动力黏度，Pa·s |
| `MaterialProperties::thermalConductivity` | 导热系数，W/(m·K) |
| `MaterialProperties::specificHeatCapacity` | 质量定压比热容，J/(kg·K) |
| `MaterialProperties::speciesDiffusivities` | 分子扩散系数，m²/s |

`src/material/` 中实现三类相互独立的材料模型：

1. `ConstantMaterial`：常密度、常黏度纯流体，保持原有算例兼容；
2. `IdealGasMixture`：根据式（5.145）计算密度，按 Wilke 规则计算黏度，按式（5.148）计算导热系数，按质量加权计算比热容；扩散系数可选择常数或式（5.149）的混合物平均近似；
3. `VolumeFractionMixture`：以 $C_k$ 加权计算 $\rho_m=\sum C_k\rho_k$、$\mu_m=\sum C_k\mu_k$ 和本阶段的算术平均导热系数，以 $\sum C_k\rho_kc_{p,k}/\rho_m$ 计算质量比热容。该类只提供物性，不表示 VOF 或 Euler–Euler 方程已经实现。

式（5.149）在纯组分极限下分子和分母同时趋于零。当前代码在该极限采用对应组分的已配置参考扩散系数，以避免奇异值；这只是数值极限处理，高压或非理想体系仍需更完整的扩散模型。

`src/flow/MultiSpeciesModel.*` 负责质量分数场、边界、时间历史、对流—扩散、通量封闭、湍流 Schmidt 数和收敛诊断。`src/flow/multiphase/` 定义统一的 `MultiphaseModel` 生命周期和 `MultiphaseContext`，并分别给出 `VofModel`、`AlgebraicSlipModel`、`EulerEulerModel`、`DpmModel` 接口。上下文只传递网格、场、质量通量、材料和多组分模块的只读或受控访问，不让某一多相模型继承或修改 SIMPLE 内部矩阵。

当前 `SimpleSolver` 的模块调用顺序为：

1. 初始化材料、可选组分模块和可选多相模块；
2. 外迭代开始，保存各模块上一次的场；
3. 求解动量预测和压力修正，更新混合物质量通量；
4. 调用可选多相模块，随后求解湍流方程；
5. 调用可选多组分模块，使用当前质量通量和湍流扩散系数；
6. 检查速度、压力、连续性、湍流、组分更新量和组分封闭，必要时继续外迭代；
7. 时间步收敛后统一推进各模块的时间历史。

本阶段只有**惰性组分的常密度、常黏度被动耦合**接入 SIMPLE。`rho_mix`、`mu_mix` 等材料输出用于物性验证和未来耦合，不能被误解为当前动量方程已经使用变密度。组分—密度低马赫反馈需要在压力修正中实现式（5.142）的散度源项；可压缩能量、化学反应和完整 Maxwell–Stefan 求解也不在本阶段范围内。VOF、代数滑移、Euler–Euler 和 DPM 的具体求解函数只保留接口，选择这些模型会明确抛出“未实现”异常，而不是静默退化为单相流。

### 5.11.9 本阶段采用的标准验证问题

**二元 Stefan 管。** 对具有相同摩尔质量的两种理想气体、常密度和常二元扩散系数，令 $A$ 从 $x=0$ 向 $x=L$ 扩散，$B$ 的总通量为零。由总质量通量与 Fick 扩散通量之和可得

$$
U=\frac{D}{L}\ln\frac{1-y_{A,L}}{1-y_{A,0}},
\qquad
y_A(x)=1-(1-y_{A,0})\exp\left(\frac{Ux}{D}\right),
\qquad y_B=1-y_A .
\tag{5.149c}
$$

本项目取 $L=0.1$ m、$D=2\times10^{-5}$ m²/s、$y_{A,0}=0.8$、$y_{A,L}=0.1$，因而 $U=3.00815479355\times10^{-4}$ m/s。该解析解同时检验对流、扩散、边界和两组分通量封闭。更一般的丙酮—甲醇—空气 Stefan 管包含不相等的分子量和三元 Maxwell–Stefan 扩散，可参考 [COMSOL Stefan Tube 模型文档](https://doc.comsol.com/6.1/doc/com.comsol.help.models.chem.stefan_tube/stefan_tube.html)；当前简化模型不与该三元实验结果混为一谈。

**三组分瞬态扩散。** 初始均匀介质在 $x=0$ 处突然施加另一固定组成，半无限域的解为

$$
y_k(x,t)=y_{k,b}+(y_{k,0}-y_{k,b})
\operatorname{erf}\left(\frac{x}{2\sqrt{Dt}}\right).
\tag{5.149d}
$$

验证算例取三种等摩尔质量惰性组分和相同 $D=10^{-3}$ m²/s，使解析解与通量修正相容；$t=0.005$ s 时右端边界影响可忽略。三个组分同时求解可检查 BDF2 时间离散、空间精度、质量分数封闭和物性随组成的变化。相同二元扩散系数是该解析算例的必要假设，不是多组分模型的一般限制。

**非零扩散修正的解析二元问题。** 为检验不同扩散系数下的非零修正项，令二元混合物无宏观流动，$D_A\ne D_B$，两端固定质量分数。由式（5.149a）得到

$$
J_A=-\rho\left[D_A(1-y_A)+D_By_A\right]\frac{\mathrm dy_A}{\mathrm dx}=\text{常数}.
$$

记 $y_0=y_A(0)$、$y_L=y_A(L)$，则解析隐式剖面为

$$
\frac{x}{L}=
\frac{D_A(y_A-y_0)+\frac12(D_B-D_A)(y_A^2-y_0^2)}
{D_A(y_L-y_0)+\frac12(D_B-D_A)(y_L^2-y_0^2)}.
\tag{5.149e}
$$

验证算例取 $D_A=2\times10^{-5}$ m²/s、$D_B=5\times10^{-6}$ m²/s、$y_0=0.8$、$y_L=0.1$，并用单调二分反演式（5.149e）与数值结果比较。该问题直接覆盖非零修正通量，不能由相同扩散系数的 Stefan 管替代。

**理想气体与 Wilke 物性。** 采用 300 K、101325 Pa 下摩尔分数为 0.79/0.21 的 N₂/O₂ 均匀混合物，分别核对质量—摩尔分数转换、式（5.145）、式（5.146）—（5.149）和质量加权比热容。均匀场避免了被动载体的近似影响物性验证。

经典三元非 Fick 扩散实验可进一步采用 Duncan 和 Toor 的 H₂–N₂–CO₂ 双球实验（[NASA ADS 文献记录](https://ui.adsabs.harvard.edu/abs/1962AIChE...8...38D/abstract)），但必须先实现并验证 Maxwell–Stefan 矩阵通量，不能用当前独立 Fick 通量替代该实验的物理机制。

## 5.12 化学反应流

### 5.12.1 燃烧形式与控制方程

化学反应通过物质转化和能量释放或吸收，改变温度、组分、压力、密度及流动特性。典型燃烧形式包括：

- **层流燃烧**：火焰锋面较平滑，传播由分子扩散和化学动力学共同控制；
- **湍流燃烧**：湍流增强混合和传热传质，火焰褶皱，并可能局部熄灭和再点燃；
- **爆轰**：燃烧波以超音速传播，激波压缩与快速放热相互维持，压力峰值和网格要求远高于普通燃烧。

化学反应流的基本未知量仍是 $\mathbf V$、$p$、$T$ 和 $y_k$，控制方程为式（5.139）—（5.145）。关键在于计算化学源项 $\dot\omega_k$ 及相应能量源项，并处理流动、扩散和化学时间尺度之间的巨大差异。

### 5.12.2 总包反应、基元反应与 Arrhenius 定律

一个总包反应可表示为

$$
\mathrm F+a\,\mathrm{Ox}\rightarrow b\,\mathrm{Pr}.
$$

若反应速率 $q=k_G(T)[X_F]^n[X_{Ox}]^m$，则各摩尔浓度变化率为

$$
\frac{\mathrm d[X_F]}{\mathrm dt}=-q,
\qquad
\frac{\mathrm d[X_{Ox}]}{\mathrm dt}=-aq,
\qquad
\frac{\mathrm d[X_{Pr}]}{\mathrm dt}=bq.
\tag{5.150}
$$

总包反应是整体计量关系，不一定代表真实基元过程。基元反应则描述实际分子事件，燃烧中常见双分子反应。温度相关速率常数通常采用 Arrhenius 形式：

$$
k_f(T)=AT^b\exp\left(-\frac{E_A}{R_uT}\right).
\tag{5.151}
$$

指前因子 $A$、温度指数 $b$ 和活化能 $E_A$ 必须来自同一套反应机理，并注意浓度单位、压力单位和摩尔质量单位的一致性。

### 5.12.3 多步反应机理的源项

第 $i$ 个可逆基元反应写为

$$
\sum_{j=1}^N\nu_{ji}'X_j
\rightleftharpoons
\sum_{j=1}^N\nu_{ji}''X_j,
\qquad i=1,\ldots,L.
\tag{5.152}
$$

其反应进度率为

$$
q_i=k_{fi}\prod_{j=1}^N[X_j]^{\nu_{ji}'}
-k_{ri}\prod_{j=1}^N[X_j]^{\nu_{ji}''}.
\tag{5.153}
$$

组分 $j$ 的摩尔净生成率和质量净生成率分别为

$$
\dot\omega_{j,mol}=\sum_{i=1}^L\nu_{ji}q_i,
\qquad
\nu_{ji}=\nu_{ji}''-\nu_{ji}',
\qquad
\dot\omega_j=M_j\dot\omega_{j,mol}.
\tag{5.154}
$$

能量方程中的化学热源为

$$
S_{chem}=-\sum_{j=1}^Nh_{f,j}^0\dot\omega_j.
\tag{5.155}
$$

正、逆反应速率常数可以通过热力学平衡常数关联。详细机理中的快速可逆反应会造成显著刚性，通常需要隐式 ODE 积分、雅可比矩阵或适当的分裂方法，而不能简单使用过大的显式化学时间步。

### 5.12.4 层流火焰与 Damköhler 数

层流火焰厚度通常很小，传播速度受压力、温度、混合物组成和反应动力学影响。流动与化学时间尺度之比用 Damköhler 数表示：

$$
Da=\frac{\tau_{flow}}{\tau_{chem}},
\qquad
\tau_{flow}=\frac{\ell_0}{v_m'},
\qquad
\tau_{chem}=\frac{\delta_L}{S_L},
\tag{5.156}
$$

其中 $S_L$ 和 $\delta_L$ 分别为层流火焰速度和厚度。

- $Da\gg1$：化学反应相对流动很快，反应区可在局部快速完成；
- $Da\ll1$：化学反应较慢，组分和温度可在反应前经历较长输运过程；
- $Da\sim1$：流动混合与化学反应均不可忽略，需要同时考虑两者。

预混和非预混是燃料与氧化剂初始混合状态的分类，不能仅凭 $Da$ 的大小定义。非预混燃烧通常更受混合控制，而预混火焰传播对化学动力学和火焰结构十分敏感。

### 5.12.5 湍流燃烧模式与有限速率模型

根据火焰厚度、积分尺度 $\ell_0$ 和 Kolmogorov 尺度 $\ell_K$，可初步区分：

| 模式 | 典型判据 | 主要特征 |
|---|---|---|
| 褶皱层流火焰 | $Da>1$，$\delta_L\ll\ell_K$ | 反应区薄，湍流主要使火焰面褶皱 |
| 旋涡小火焰 | $Da\sim1$，$\ell_0>\delta_L>\ell_K$ | 小火焰结构与涡相互作用，反应区增厚 |
| 分布反应 | $Da<1$，$\delta_L>\ell_0$ | 反应时间超过主要涡时间尺度，反应分布更分散 |

**有限速率模型**直接使用式（5.153）和式（5.154）计算化学源项。它适合层流、复杂有限速率化学以及需要保留不同化学时间尺度的问题。在湍流中直接对平均组分和平均温度代入 Arrhenius 公式，会忽略脉动与强非线性化学之间的耦合，通常需要额外的湍流—化学相互作用模型。

### 5.12.6 涡耗散模型

涡耗散模型（EDM）假设化学反应足够快，速率由湍流混合控制，从而避免直接计算 Arrhenius 动力学。对反应 $r$，其速率幅值可概括为

$$
\mathcal R_r=A\rho\frac{\varepsilon}{k}
\min\left\{
\min_{\mathcal R}\left[
\frac{Y_{\mathcal R}}{\nu_{\mathcal R,r}'M_{\mathcal R}}
\right],
B\frac{\sum_{\mathcal P}Y_{\mathcal P}}
{\sum_{\mathcal P}\nu_{\mathcal P,r}''M_{\mathcal P}}
\right\},
\tag{5.157}
$$

组分质量源项由 $(\nu_{i,r}''-\nu_{i,r}')M_i\mathcal R_r$ 计算。反应物项限制燃烧不能越过缺乏反应物的区域，产物项可用于避免未混合区域过早反应。模型常数 $A$、$B$ 按具体实现选取。

EDM 成本低，适合快速湍流非预混燃烧和单步或两步总包机理。由于各反应共享湍流混合速率，它不适合需要区分多种化学时间尺度的详细机理。LES 中可将混合时间尺度改为亚格子尺度估计，例如

$$
\tau_{sgs}^{-1}=\sqrt{2S_{ij}S_{ij}},
\qquad
S_{ij}=\frac12\left(\frac{\partial u_i}{\partial x_j}
+\frac{\partial u_j}{\partial x_i}\right).
\tag{5.158}
$$

### 5.12.7 涡耗散概念模型

涡耗散概念模型（EDC）同时考虑湍流混合和有限速率化学，特别适合 $Da\sim1$ 的情形。其基本假设为：

1. 化学反应发生在占单元体积一小部分的精细结构中；
2. 精细结构可视为恒压反应器，以单元当前组分和温度为初值；
3. 在精细尺度时间 $\tau^*$ 内用 ODE 求解器积分化学动力学，得到 $Y_i^*$。

组分源项写为

$$
R_i=\rho\kappa\frac{Y_i^*-Y_i}{\tau^*},
\tag{5.159}
$$

其中 $\kappa$ 为精细结构体积分数。经典尺度估计为

$$
\kappa=\frac{\xi^{*2}}{1-\xi^{*2}},
\qquad
\xi^*=C_\xi\left(\frac{\nu\varepsilon}{k^2}\right)^{1/4},
\qquad
\tau^*=C_\tau\left(\frac{\nu}{\varepsilon}\right)^{1/2},
\tag{5.160}
$$

课件采用 $C_\xi=2.1377$、$C_\tau=0.4082$。实际实现还需限制 $\kappa$ 和 $\tau^*$ 的数值范围，并保证与化学积分时间一致。

部分平均搅拌反应器（PaSR）形式的 EDC 可使用

$$
\kappa=\frac{t_c}{t_c+t_{mix}}=\frac1{1+Da},
\qquad
Da=\frac{t_{mix}}{t_c},
\qquad
\tau^*=\min(t_c,t_{mix}),
\tag{5.161}
$$

其中 $t_c$ 由局部有限速率化学估计，$t_{mix}$ 由湍流时间尺度估计。与标准 EDC 相比，PaSR 通过化学与混合时间直接调节反应体积分数和积分时间，但模型定义不能混用。

### 5.12.8 爆燃转爆轰与数值要求

爆燃转爆轰（DDT）涉及强间断、非线性放热和激波—化学耦合。主要困难包括：

- 激波、接触间断和反应区需要高分辨率捕捉；
- 爆轰胞格和诱导区尺度很小，网格不足会改变传播速度甚至误判是否起爆；
- Arrhenius 速率对温度高度敏感，数值误差会被非线性放大；
- 三维大规模计算需要并行和自适应加密。

一种简化的诱导变量 $\delta$ 模型可写为

$$
\frac{\partial\rho\delta}{\partial t}
+\frac{\partial\rho\delta u_j}{\partial x_j}
=-k_\delta\rho\exp\left(-\frac{E_\delta}{RT}\right),
\tag{5.162}
$$

用于描述诱导过程的消耗。实际 DDT 计算通常需要守恒型可压缩 Euler 或 Navier–Stokes 求解器、高分辨率界面通量、合适的激波限制器、严格时间步控制和 AMR。普通不可压缩 SIMPLE 压力修正框架不能直接处理强激波和超音速爆轰波。

### 5.12.9 化学反应流的综合求解与验证

一个含化学反应的时间步可组织为：

1. 求解连续性和动量方程，获得满足质量约束的速度与压力；
2. 求解组分对流—扩散，并与化学源项耦合；
3. 用所选有限速率、EDM、EDC 或其他模型计算化学源项；
4. 求解包含化学热和组分焓扩散的能量方程；
5. 更新密度、温度相关物性、湍流量和扩散系数；
6. 重复外迭代，直至组分、温度、压力和反应速率稳定。

验证时应检查：

- 元素守恒、总质量守恒和能量守恒；
- $0\leq y_k\leq1$ 及 $\sum_ky_k=1$；
- 化学积分误差、时间步和网格敏感性；
- 层流火焰速度、点火延迟、火焰厚度或实验温度剖面；
- 湍流模型与燃烧模型的适用范围是否一致。

喷管两步煤油燃烧等算例可联合使用 PIMPLE、$k$–$\varepsilon$、总焓、组分输运和 EDM，但模型组合能够运行并不等于化学预测已经得到验证。最终仍需与实验或更高保真度结果比较。

## 本章小结

多相流模型的选择首先取决于界面尺度和相间运动。VOF 通过体积分数输运和界面重构解析自由表面，CICSAM 属于代数压缩格式，PLIC 属于几何重构。ASM 用漂移速度补充混合物模型，Euler–Euler 分别求解各相方程并用阻力、升力、虚拟质量力和传质耦合，Euler–Lagrange 则追踪颗粒轨迹并把交换量反馈给连续相。

多组分流动以共享混合物速度和扩散通量为基础。组分方程必须满足质量分数归一、扩散通量总和为零和元素守恒。可压缩多组分能量方程还需包含组分扩散携带的焓及生成焓对应的化学能源。

化学反应源项可由总包或基元机理、Arrhenius 定律和多步反应网络计算。有限速率模型保留化学动力学，EDM 以湍流混合控制快速反应，EDC 通过精细结构反应器结合混合与有限速率化学。DDT 还涉及强激波和反应区耦合，需要高分辨率可压缩算法和 AMR，不能简单套用不可压缩 SIMPLE 框架。

## 习题

1. 比较多相流与多组分流动中“不同物质”的数学描述，说明何时各相可以具有不同速度。
2. 从指示函数推导体积分数输运方程，并证明无源不可压缩条件下纯相单元应保持 $C_k=0$ 或1。
3. 推导式（5.101），说明为什么普通一阶迎风不适合长期保持清晰 VOF 界面。
4. 解释 CICSAM 中 HYPER-C 与 UQ 混合的作用，以及角度因子如何影响界面压缩。
5. 描述 PLIC 的重构和通量计算两个步骤。为什么几何通量通常比单元中心插值更有利于相体积守恒？
6. 从 $\mathbf U_{dr,k}=\mathbf U_k-\mathbf U_m$ 推导 ASM 动量方程中的漂移应力项。
7. 推导 Schiller–Naumann 阻力修正在 $Re>1000$ 时近似为 $0.0183Re$ 的条件。
8. 对 Euler–Euler 两相动量方程求和，说明共享压力和相间作用—反作用关系如何保证总动量一致。
9. 比较阻力、升力和虚拟质量力的物理来源及各自适用条件。
10. 从式（5.125）推导式（5.127），并验证 $F_D$ 的量纲。
11. 写出一个蒸发液滴的质量和能量更新步骤，说明蒸发潜热的符号及气相源项如何与之配对。
12. 为什么使用统计颗粒包时必须乘以颗粒包权重？若遗漏会造成什么结果？
13. 对所有组分方程求和，证明 $\sum_k\mathbf j_k=0$ 和 $\sum_k\dot\omega_k=0$ 的必要性。
14. 推导理想气体混合物的 $R_{mix}=R_u\sum_ky_k/M_k$，并说明质量分数和摩尔分数如何转换。
15. 解释多组分总焓方程中 $-\nabla\cdot\sum_kh_k\mathbf j_k$ 的物理意义。
16. 对反应 $\mathrm F+a\mathrm{Ox}\rightarrow b\mathrm{Pr}$，验证各组分的质量生成率满足总质量守恒，并写出相应元素守恒条件。
17. 为什么 Arrhenius 源项不能直接对平均温度求值后就代表湍流平均反应速率？
18. 比较有限速率、EDM 和 EDC 对化学时间尺度与湍流时间尺度的处理。
19. 分析 $Da\gg1$、$Da\sim1$ 和 $Da\ll1$ 时数值求解中应重点关注的误差来源。
20. 说明强爆轰计算为什么需要高分辨率守恒格式和 AMR，以及仅加密网格为何仍不足以保证正确起爆。

## 延伸阅读

1. Hirt, C. W., and Nichols, B. D. Volume of fluid (VOF) method for the dynamics of free boundaries. *Journal of Computational Physics*, 1981, 39: 201–225.
2. Brackbill, J. U., Kothe, D. B., and Zemach, C. A continuum method for modeling surface tension. *Journal of Computational Physics*, 1992, 100: 335–354.
3. Ubbink, O., and Issa, R. I. A method for capturing sharp fluid interfaces on arbitrary meshes. *Journal of Computational Physics*, 1999, 153: 26–50.
4. Crowe, C. T., Schwarzkopf, J. D., Sommerfeld, M., and Tsuji, Y. *Multiphase Flows with Droplets and Particles*. 2nd ed., CRC Press, 2012.
5. Bird, R. B., Stewart, W. E., and Lightfoot, E. N. *Transport Phenomena*. 2nd ed., Wiley, 2002.
6. Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*. Wiley, 2003.
7. Poinsot, T., and Veynante, D. *Theoretical and Numerical Combustion*. 3rd ed., 2011.
8. Peters, N. *Turbulent Combustion*. Cambridge University Press, 2000.
9. Wilke, C. R. A viscosity equation for gas mixtures. *Journal of Chemical Physics*, 1950, 18: 517–519. [DOI: 10.1063/1.1747673](https://doi.org/10.1063/1.1747673).
10. Taylor, R., and Krishna, R. *Multicomponent Mass Transfer*. Wiley, 1993.
11. Duncan, J. B., and Toor, H. L. An experimental study of three component gas diffusion. *AIChE Journal*, 1962, 8: 38–41. [NASA ADS record](https://ui.adsabs.harvard.edu/abs/1962AIChE...8...38D/abstract).
12. Crank, J. *The Mathematics of Diffusion*. 2nd ed., Oxford University Press, 1975.
13. Krishna, R., and Wesselingh, J. A. The Maxwell–Stefan approach to mass transfer. *Chemical Engineering Science*, 1997, 52: 861–911.
14. [COMSOL Application Library: Stefan Tube](https://doc.comsol.com/6.1/doc/com.comsol.help.models.chem.stefan_tube/stefan_tube.html). 包含三元 Stefan 管的模型方程、二元扩散系数和实验对比设置。
