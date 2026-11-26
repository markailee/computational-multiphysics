# 第5章（续）中子输运与扩散方程的数值求解

中子输运计算的目标是求得中子角通量及其积分量，并进一步确定有效增殖因数、反应率、功率分布、泄漏和动力学参数。与第4章的流体标量输运相比，中子输运多了**飞行方向**这一独立变量；散射和裂变又把不同方向、不同能量群耦合起来。因此，一个完整算法必须同时处理时间、空间、方向、能量和本征值五个层面。

本章从一般 Boltzmann 中子输运方程出发，依次介绍多群离散纵标（$S_N$）形式、$P_1$ 扩散近似、有限体积离散、边界条件、固定源与 $k$ 本征值迭代，以及通量收敛后的物理量计算。重点推导单元中心有限体积法，同时说明其与有限差分法（FDM）和有限元法（FEM）的对应关系。最后讨论 Reed、IAEA、Takeda、TWIGL、LRA、Kobayashi 和 C5G7 等经典数值算例。公式编号承接前一篇声学章节。

本章用 $\psi$ 表示角通量，$\phi$ 表示标量通量，$\mathbf J$ 表示净中子流密度；用 $v$ 表示中子速率，用 $\Sigma$ 表示宏观截面，用 $D$ 表示扩散系数。符号 $\mathsf A$ 和 $\mathsf F$ 分别表示离散后的非裂变算子和裂变算子，粗体 $\boldsymbol\Phi$ 或 $\boldsymbol\Psi$ 表示全部离散未知量。除特别说明外，界面外法向均由当前单元指向外部。

## 5.25 中子输运控制方程与常用近似

### 5.25.1 角通量、标量通量和中子流

设 $\psi(\mathbf r,E,\boldsymbol\Omega,t)$ 为位置 $\mathbf r$、能量 $E$、飞行方向 $\boldsymbol\Omega$ 和时间 $t$ 的中子角通量。若采用“每单位能量、每单位立体角”的定义，其量纲可写为

$$
[\psi]=\mathrm{cm^{-2}\,s^{-1}\,energy^{-1}\,sr^{-1}}.
$$

方向积分得到标量通量和净中子流密度：

$$
\boxed{
\phi(\mathbf r,E,t)=\int_{4\pi}\psi\,\mathrm d\Omega,
\qquad
\mathbf J(\mathbf r,E,t)=\int_{4\pi}\boldsymbol\Omega\psi\,\mathrm d\Omega
}.
\tag{5.268}
$$

标量通量可以理解为单位时间内穿过以该点为中心的小球的所有中子轨迹长度密度；它不是某个特定表面上的单向计数。穿过分界面的净流动由 $\mathbf J\cdot\mathbf n$ 给出，而反应率由相应截面与标量通量的乘积给出。角通量与角密度 $n_\psi$ 的关系是 $\psi=v n_\psi$；使用其他通量定义时，必须同步调整源项、角度权重和量纲。

宏观截面 $\Sigma_x$ 的单位通常为 $\mathrm{cm^{-1}}$，表示中子在单位飞行距离内发生反应 $x$ 的概率。$\Sigma_t$、$\Sigma_a$、$\Sigma_s$ 和 $\Sigma_f$ 分别表示总、吸收、散射和裂变截面。通常 $\Sigma_a$ 包含裂变消失道；$\nu\Sigma_f$ 才表示每次裂变平均产生的中子数所对应的产生截面。

### 5.25.2 一般时间相关 Boltzmann 输运方程

在线性、静止介质且核数据不随中子通量本身改变的条件下，一般中子输运方程写为

$$
\boxed{
\begin{aligned}
&\frac1{v(E)}\frac{\partial\psi}{\partial t}
+\boldsymbol\Omega\cdot\nabla\psi
+\Sigma_t(\mathbf r,E)\psi\\
&\quad=
\int_0^\infty\!\!\int_{4\pi}
\Sigma_s(\mathbf r,E'\!\to E,\boldsymbol\Omega'\!\to\boldsymbol\Omega)
\psi(\mathbf r,E',\boldsymbol\Omega',t)
\,\mathrm d\Omega'\,\mathrm dE'\\
&\qquad+
\frac{\chi_p(E)}{4\pi}
\int_0^\infty \nu_p(E')\Sigma_f(\mathbf r,E')
\phi(\mathbf r,E',t)\,\mathrm dE'\\
&\qquad+
\frac1{4\pi}\sum_{i=1}^{I}\chi_{d,i}(E)\lambda_i C_i(\mathbf r,t)
+Q(\mathbf r,E,\boldsymbol\Omega,t).
\end{aligned}}
\tag{5.269}
$$

左端依次是中子数随时间的变化、定向自由飞行的净流出和总碰撞消失；右端依次为散射进入当前能量—方向状态、瞬发裂变源、缓发中子源和独立外源。这里 $\chi_p$ 与 $\chi_{d,i}$ 分别为瞬发谱和第 $i$ 组缓发谱，$\lambda_i$ 为缓发先驱核衰变常数，$C_i$ 为先驱核浓度。缓发先驱核满足

$$
\boxed{
\frac{\partial C_i}{\partial t}
=\int_0^\infty \beta_i(E')\nu(E')\Sigma_f(\mathbf r,E')
\phi(\mathbf r,E',t)\,\mathrm dE'
-\lambda_i C_i
}.
\tag{5.270}
$$

$\beta_i$ 是第 $i$ 组缓发中子份额，$\nu$ 是每次裂变的总平均产额。式（5.269）采用 $\nu_p$ 表示瞬发产额，以避免与总产额重复计数。若使用单一缓发份额模型，常写成 $\nu_p=(1-\beta)\nu$、$\beta=\sum_i\beta_i$；若核数据库分别给出瞬发和缓发数据，则应直接遵循数据库定义。能谱通常按 $\int\chi(E)\,\mathrm dE=1$ 归一化。

式（5.269）中的 $\boldsymbol\Omega\cdot\nabla\psi$ 是中子沿直线自由飞行产生的**方向性对流项**，并不是流体携带中子的宏观速度项。对一般流体标量方程作形式对应时，可令存储系数为 $1/v$、广义对流矢量为 $\boldsymbol\Omega$；若坚持写成 $\rho\mathbf u$ 的形式，则相应 $\rho=1/v$、$\mathbf u=v\boldsymbol\Omega$。运动介质中的相对论速度变换、碰撞核运动以及先驱核随流体迁移需要额外模型，不能简单在中子输运方程中任意加入一个 $\mathbf u\cdot\nabla\psi$。

### 5.25.3 多群 $S_N$ 方程和 $k$ 本征值形式

确定性计算通常先把能量区间划分为 $G$ 个能群，再对方向变量进行离散。离散纵标法选取 $M$ 个方向 $\boldsymbol\Omega_m$ 和权重 $w_m$，使

$$
\int_{4\pi}f(\boldsymbol\Omega)\,\mathrm d\Omega
\approx\sum_{m=1}^{M}w_m f(\boldsymbol\Omega_m),
\qquad \sum_{m=1}^{M}w_m=4\pi.
$$

群 $g$、方向 $m$ 的稳态多群方程为

$$
\boxed{
\begin{aligned}
&\boldsymbol\Omega_m\cdot\nabla\psi_{g,m}
+\Sigma_{t,g}\psi_{g,m}\\
&\quad=
\sum_{g'=1}^{G}\sum_{n=1}^{M}w_n
\Sigma_{s,g'\to g}(\boldsymbol\Omega_n\to\boldsymbol\Omega_m)
\psi_{g',n}\\
&\qquad+\frac{\chi_g}{4\pi k_{\mathrm{eff}}}
\sum_{g'=1}^{G}\nu\Sigma_{f,g'}\phi_{g'}
+Q_{g,m},
\qquad
\phi_g=\sum_{m=1}^{M}w_m\psi_{g,m}.
\end{aligned}}
\tag{5.271}
$$

这里已假定裂变中子方向各向同性。若保留时间项，只需在左端加入 $v_g^{-1}\partial\psi_{g,m}/\partial t$，并恢复式（5.270）的缓发源。$S_N$ 角度求积应准确积分足够阶数的多项式，并与散射矩展开和反射方向相匹配。柱坐标和球坐标中的方向会随空间位置改变，还会出现角度再分布项，不能直接照搬上述直角坐标形式。

对于无独立外源的临界问题，$Q_{g,m}=0$，$k_{\mathrm{eff}}$ 是为使基波解保持稳态而引入的裂变源缩放因子。$k_{\mathrm{eff}}>1$、$=1$ 和 $<1$ 分别对应超临界、临界和次临界。这里的 $1/k_{\mathrm{eff}}$ 是**本征值问题的数学构造**，不是瞬态计算中每一步都要对真实裂变截面施加的物理修正。

积分整个区域后，还可以把有效增殖因数解释为

$$
\boxed{
k_{\mathrm{eff}}=\frac{\mathcal P_f}{\mathcal L+\mathcal A},
\qquad
\mathcal P_f=\int_V\!\int_0^\infty\nu\Sigma_f\phi\,\mathrm dE\,\mathrm dV
}.
\tag{5.272}
$$

其中 $\mathcal L$ 是穿过外边界的净流出率，$\mathcal A$ 是总吸收率。该等式应与离散本征方程的守恒平衡一致；存在独立外源时，简单用总产生量除以总损失量并不能替代式（5.271）的本征值求解。

### 5.25.4 $P_1$ 近似与 Fick 扩散律

当角通量只保留前两个角矩时，

$$
\boxed{
\psi(\mathbf r,\boldsymbol\Omega)
\approx\frac1{4\pi}
\left[\phi(\mathbf r)+3\boldsymbol\Omega\cdot\mathbf J(\mathbf r)\right]
}.
\tag{5.273}
$$

对单群或已经整理好的群方程分别取零阶和一阶角矩，可得

$$
\boxed{
\frac1v\frac{\partial\phi}{\partial t}
+\nabla\cdot\mathbf J+\Sigma_t\phi=S_0,
\qquad
\frac1v\frac{\partial\mathbf J}{\partial t}
+\frac13\nabla\phi+\Sigma_{\mathrm{tr}}\mathbf J=\mathbf S_1
}.
\tag{5.274}
$$

$S_0$ 包含散射、裂变和外源的零阶矩，$\mathbf S_1$ 是相应的一阶源矩。若一阶中子流比通量更快达到准稳态，可忽略 $v^{-1}\partial\mathbf J/\partial t$，并在线性各向异性散射、各向同性源条件下得到 Fick 定律：

$$
\boxed{
\mathbf J=-D\nabla\phi,
\qquad
D=\frac1{3\Sigma_{\mathrm{tr}}},
\qquad
\Sigma_{\mathrm{tr}}=\Sigma_t-\Sigma_{s,1}
}.
\tag{5.275}
$$

$\Sigma_{s,1}$ 是一阶散射矩；在线性各向异性单群散射下可写为 $\bar\mu_0\Sigma_s$，其中 $\bar\mu_0$ 是平均散射余弦。多群计算应使用与群常数库相一致的输运修正定义，而不能在未核对约定时把全部群间散射机械地代入 $\Sigma_{\mathrm{tr}}$。

由此得到常用的稳态多群扩散本征值方程：

$$
\boxed{
-\nabla\cdot(D_g\nabla\phi_g)
+\Sigma_{r,g}\phi_g
=\sum_{g'\ne g}\Sigma_{s,g'\to g}\phi_{g'}
+\frac{\chi_g}{k_{\mathrm{eff}}}
\sum_{g'=1}^{G}\nu\Sigma_{f,g'}\phi_{g'}
+Q_g
},
\tag{5.276}
$$

其中移除截面为

$$
\boxed{
\Sigma_{r,g}=\Sigma_{t,g}-\Sigma_{s,g\to g}
=\Sigma_{a,g}+\sum_{g'\ne g}\Sigma_{s,g\to g'}
}.
\tag{5.277}
$$

第二个等号适用于通常的二体散射群平衡；$(n,2n)$ 等多中子产生反应若单独建模，应按实际粒子数平衡增加相应产生项。扩散近似把未知量从每个“单元—能群—方向”减少为每个“单元—能群”，所以适合全堆扩散和粗网格计算。扩散项不是对式（5.269）再额外添加的物理机制，而是对定向输运的闭合适似；在同一未知量上同时保留完整对流算子和独立扩散算子会造成重复建模。

在真空、强吸收体、窄缝、点源附近或只有几个平均自由程的区域内，角通量通常高度各向异性，$P_1$ 假设可能失效。此时应使用 $S_N$、高阶 $P_N$/SP$_N$、特征线法或其他输运模型，而不是仅靠加密网格弥补模型误差。

### 5.25.5 FVM、FDM 和 FEM 共用的代数形式

三种确定性空间离散最终都把固定源或本征值问题写成

$$
\boxed{
\mathsf A\mathbf x=\frac1{k_{\mathrm{eff}}}\mathsf F\mathbf x+\mathbf q
}.
\tag{5.278}
$$

$\mathbf x$ 可以是包含所有方向的多群角通量 $\boldsymbol\Psi$，也可以只是多群标量通量 $\boldsymbol\Phi$。$\mathsf A$ 包含存储、对流、扩散、移除以及移入左端的散射项，$\mathsf F$ 包含裂变谱和 $\nu\Sigma_f$。临界本征值问题中 $\mathbf q=0$；固定源问题通常令 $k_{\mathrm{eff}}=1$ 并保留实际外源，或根据明确的次临界增殖模型处理裂变项。

- **FVM** 对每个控制体积分守恒方程，通过共享界面通量保证局部和全局守恒，适合复杂材料分区及非结构网格。
- **FDM** 在结构网格节点或单元中心直接近似导数，模板简洁；规则网格上常可得到与 FVM 相同的系数，但材料界面和边界处理必须与守恒通量一致。
- **FEM** 从弱形式出发，在单元上用形函数表示通量，适合非结构几何和高阶近似。扩散方程常用连续 Galerkin FEM；一阶 $S_N$ 对流方程常用间断 Galerkin、迎风稳定化、最小二乘或自共轭角通量形式，不能把无稳定化的中心 Galerkin 格式当作普遍可靠选择。

因此，“采用 FVM/FEM/FDM”并没有改变中子输运的物理近似层次。必须分别说明：采用完整输运还是扩散模型、如何离散能量和方向，以及最后如何离散空间。

## 5.26 输运与扩散方程的有限体积离散

### 5.26.1 $S_N$ 控制体平衡

在离散方向上，标量通量和当前矩采用求积公式

$$
\boxed{
\phi_{g,P}=\sum_{m=1}^{M}w_m\psi_{g,m,P},
\qquad
\mathbf J_{g,P}=\sum_{m=1}^{M}w_m\boldsymbol\Omega_m\psi_{g,m,P}
}.
\tag{5.279}
$$

对单元 $P$ 积分式（5.271），并使用后向 Euler 时间离散，得到

$$
\boxed{
\frac{V_P}{v_g\Delta t}
\left(\psi_{g,m,P}^{n+1}-\psi_{g,m,P}^{n}\right)
+\sum_{f\in\mathcal F_P}A_f\mu_{m,f}\psi_{g,m,f}^{n+1}
+V_P\Sigma_{t,g,P}\psi_{g,m,P}^{n+1}
=V_Pq_{g,m,P}
},
\tag{5.280}
$$

其中

$$
\mu_{m,f}=\boldsymbol\Omega_m\cdot\mathbf n_f,
\qquad
q_{g,m,P}=S^{\mathrm{scat}}_{g,m,P}
+S^{\mathrm{fis}}_{g,m,P}+Q_{g,m,P}.
$$

$q_{g,m,P}$ 在当前源迭代或时间层内冻结后，左端只含一个能群、一个方向上的空间未知量。本征值计算的 $S^{\mathrm{fis}}$ 含 $1/k_{\mathrm{eff}}$；真实瞬态则按式（5.269）使用瞬发和缓发源，不能继续沿用本征值缩放，详见第5.28.5节。二维问题中 $V_P$ 换成面积、$A_f$ 换成边长，并按模型乘以统一厚度。

### 5.26.2 一阶迎风格式及单元系数

定向对流的速度就是 $\boldsymbol\Omega_m$。当 $\mu_{m,f}>0$ 时中子离开当前单元，界面值由当前单元上游给出；当 $\mu_{m,f}<0$ 时则由相邻单元给出：

$$
\boxed{
\psi_{g,m,f}=
\begin{cases}
\psi_{g,m,P},&\mu_{m,f}\ge0,\\
\psi_{g,m,N},&\mu_{m,f}<0,
\end{cases}}
\tag{5.281}
$$

$N$ 是内部界面另一侧的单元。整理式（5.280）可得到单元代数方程

$$
\boxed{
\begin{aligned}
a_{P,g,m}\psi_{g,m,P}
&=\sum_{N\in\mathcal N_P}a_{N,g,m}\psi_{g,m,N}+b_{g,m,P},\\
a_{N,g,m}&=A_f\max(-\mu_{m,f},0),\\
a_{P,g,m}&=\frac{V_P}{v_g\Delta t}+V_P\Sigma_{t,g,P}
+\sum_{f\in\mathcal F_P}A_f\max(\mu_{m,f},0),\\
b_{g,m,P}&=\frac{V_P}{v_g\Delta t}\psi_{g,m,P}^{n}
+V_Pq_{g,m,P}+b^{\mathrm{in}}_{g,m,P}.
\end{aligned}}
\tag{5.282}
$$

$b^{\mathrm{in}}$ 是已知入流边界的贡献，稳态时令 $V_P/(v_g\Delta t)=0$。内部界面在 $P$ 与 $N$ 两侧各组装一次，通量自然大小相等、符号相反。对于非负源和非负入流，该格式保持非负性，并在冻结源项时给出对角占优的单元系统；代价是一阶数值耗散。后向 Euler 不需要显式 CFL 条件来保证这种隐式平衡的可解性，但大时间步仍会产生明显时间截断误差。

在规则网格上，可以按 $\boldsymbol\Omega_m$ 各分量的符号排列单元，使所有上游值先于当前单元计算，这就是**输运扫描**。每个八分区方向对应一种扫描顺序。非结构网格需要构造上游依赖图；反射边界、循环网格依赖或全耦合散射会使单次扫描不再充分，需要迭代、直接法或 Krylov 求解。

### 5.26.3 菱形差分、指数格式与高精度通量

一维单群、单方向单元中，令入流和出流角通量为 $\psi_{\mathrm{in}}$ 和 $\psi_{\mathrm{out}}$，单元平均值为 $\bar\psi$，单元光学厚度为 $\tau=\Sigma_t h/|\mu|$。菱形差分使用

$$
\boxed{
\bar\psi=\frac12(\psi_{\mathrm{in}}+\psi_{\mathrm{out}})
}.
\tag{5.283}
$$

代入一维平衡方程 $|\mu|(\psi_{\mathrm{out}}-\psi_{\mathrm{in}})+\Sigma_t h\bar\psi=hq$，对 $\mu>0$ 可得

$$
\boxed{
\bar\psi=\frac{|\mu|\psi_{\mathrm{in}}+\dfrac h2q}
{|\mu|+\dfrac{\Sigma_t h}{2}},
\qquad
\psi_{\mathrm{out}}=2\bar\psi-\psi_{\mathrm{in}}
}.
\tag{5.284}
$$

该格式在光滑解上具有较高精度，但粗光学厚度网格上可能出现负通量和振荡。直接把负值截为零会改变粒子数平衡；若使用负通量修正，必须同时检查其守恒性和对相邻上游/下游单元的影响。

对同一单元精确积分常系数特征线方程，可得到指数格式

$$
\boxed{
\psi_{\mathrm{out}}=\psi_{\mathrm{in}}e^{-\tau}
+\frac{q}{\Sigma_t}(1-e^{-\tau}),
\qquad
\bar\psi=\psi_{\mathrm{in}}\frac{1-e^{-\tau}}{\tau}
+\frac{q}{\Sigma_t}\left[1-\frac{1-e^{-\tau}}{\tau}\right]
}.
\tag{5.285}
$$

$\Sigma_t\to0$ 时应使用极限形式或稳定的指数函数求值。指数关系也是特征线法、指数加权和高精度有限体积格式的基础。多维高阶重构可以使用线性源、二阶迎风、间断有限元或限制器，但必须确保界面两侧使用同一个守恒数值通量，并分别进行空间、方向和时间收敛试验。不能只凭局部模板的“二阶”名称宣称整个离散收敛阶。

### 5.26.4 多群扩散方程的有限体积系数

扩散方程式（5.276）在单元 $P$ 上的后向 Euler 平衡为

$$
\boxed{
\frac{V_P}{v_g\Delta t}(\phi_{g,P}^{n+1}-\phi_{g,P}^{n})
+V_P\Sigma_{r,g,P}\phi_{g,P}^{n+1}
+\sum_{f\in\mathcal F_P}A_fJ_{g,n,f}
=V_PS_{g,P}
},
\tag{5.286}
$$

其中 $J_{g,n,f}=-D_g\nabla\phi_g\cdot\mathbf n_f$ 是向外净扩散流。对正交内部界面，令界面两侧单元中心到界面的距离为 $d_{Pf}$ 和 $d_{Nf}$，$d_{PN}=d_{Pf}+d_{Nf}$。由界面通量连续可得

$$
\boxed{
J_{g,n,f}=\frac{D_{g,\mathrm{eff}}}{d_{PN}}(\phi_{g,P}-\phi_{g,N}),
\qquad
D_{g,\mathrm{eff}}
=\frac{d_{PN}}{d_{Pf}/D_{g,P}+d_{Nf}/D_{g,N}}
}.
\tag{5.287}
$$

这相当于两侧扩散导纳串联；等距离时退化为调和平均 $2D_PD_N/(D_P+D_N)$。简单算术平均 $D$ 在扩散系数差异很大时可能产生不合理界面流。非正交网格按第4章采用“主扩散项加非正交修正”，修正项可延迟显式处理，但共享界面的总通量必须守恒。

定义 $G_{g,f}=A_fD_{g,\mathrm{eff}}/d_{PN}$，则

$$
\boxed{
\begin{aligned}
a_{P,g}\phi_{g,P}
&=\sum_N G_{g,f}\phi_{g,N}+b_{g,P},\\
a_{P,g}&=\frac{V_P}{v_g\Delta t}+V_P\Sigma_{r,g,P}
+\sum_fG_{g,f}+G^{\mathrm{ext}}_{g,P},\\
b_{g,P}&=\frac{V_P}{v_g\Delta t}\phi_{g,P}^{n}
+V_PS_{g,P}+b^{\mathrm{ext}}_{g,P}.
\end{aligned}}
\tag{5.288}
$$

$G^{\mathrm{ext}}$ 和 $b^{\mathrm{ext}}$ 分别来自外边界的对角项和源项。群源可以写成

$$
\boxed{
S_{g,P}=\sum_{g'\ne g}\Sigma_{s,g'\to g,P}\phi_{g',P}
+\frac{\chi_{g,P}}{k_{\mathrm{eff}}}
\sum_{g'=1}^{G}\nu\Sigma_{f,g',P}\phi_{g',P}
+Q_{g,P}
}.
\tag{5.289}
$$

若全部群同时组装，群间散射和裂变产生矩阵的非对角块；若按群 Gauss–Seidel 迭代，则这些项暂时放入 $b_{g,P}$。热群中向上散射使群矩阵不再能按能量单向扫过，必须迭代或联立求解。

### 5.26.5 FEM 弱形式及 FDM 对应关系

对扩散方程式（5.276）乘以试探函数 $w_g$，分部积分并采用 Robin 边界 $J_n=\gamma\phi$，得到每组弱形式

$$
\boxed{
\begin{aligned}
&\int_VD_g\nabla w_g\cdot\nabla\phi_g\,\mathrm dV
+\int_V\Sigma_{r,g}w_g\phi_g\,\mathrm dV
+\int_{\partial V}\gamma w_g\phi_g\,\mathrm dA\\
&\quad-\sum_{g'\ne g}\int_V\Sigma_{s,g'\to g}w_g\phi_{g'}\,\mathrm dV\\
&\qquad=\frac1{k_{\mathrm{eff}}}\sum_{g'=1}^{G}
\int_V\chi_g\nu\Sigma_{f,g'}w_g\phi_{g'}\,\mathrm dV
+\int_Vw_gQ_g\,\mathrm dV .
\end{aligned}}
\tag{5.290}
$$

代入 $\phi_g=\sum_aN_a\phi_{g,a}$，即可得到扩散矩阵、移除质量矩阵、散射矩阵和裂变矩阵。Robin 边界在弱形式中自然产生面积分；反射边界 $\gamma=0$ 时该积分消失。瞬态项额外产生 $\int_Vw_gv_g^{-1}\dot\phi_g\,\mathrm dV$。

在规则长方体网格上，单元中心 FVM 与单元中心 FDM 通常具有相同五点和七点模板。顶点中心 FDM、有限差分节点法和 FEM 的自由度位置不同，不能仅按矩阵外观判断其材料界面处理是否等价。对 $S_N$ 方程，间断 FEM 在每个单元界面采用与式（5.281）相同的迎风数值通量，因此可以看作分片多项式有限体积思想的推广。

## 5.27 边界条件与材料界面处理

### 5.27.1 输运方程的入流、真空与反射边界

一阶双曲型输运方程只需在**入射方向**给定边界值。设 $\mathbf n$ 为区域外法向，则一般入流条件为

$$
\boxed{
\psi_{g}(\mathbf r_b,\boldsymbol\Omega_m)
=\psi^{\mathrm{in}}_{g}(\mathbf r_b,\boldsymbol\Omega_m),
\qquad \boldsymbol\Omega_m\cdot\mathbf n<0
}.
\tag{5.291}
$$

真空边界对应 $\psi^{\mathrm{in}}=0$，但出射角通量通常不为零。已知外中子束或相邻区域提供的角通量时，直接给定 $\psi^{\mathrm{in}}$。对式（5.282），边界入流贡献为

$$
b^{\mathrm{in}}_{g,m,P}
=\sum_{f\in\partial V_P\cap\Gamma}
A_f\max(-\mu_{m,f},0)\psi^{\mathrm{in}}_{g,m,f}.
$$

镜面反射要求入射方向 $\boldsymbol\Omega$ 与出射方向 $\boldsymbol\Omega^*$ 满足

$$
\boxed{
\boldsymbol\Omega^*=\boldsymbol\Omega-2(\boldsymbol\Omega\cdot\mathbf n)\mathbf n,
\qquad
\psi_b(\mathbf r_b,\boldsymbol\Omega)
=\beta\psi_b(\mathbf r_b,\boldsymbol\Omega^*)
},
\tag{5.292}
$$

其中 $\boldsymbol\Omega\cdot\mathbf n<0$，$\beta=1$ 为理想反射，$0\le\beta<1$ 可表示返照损失。离散角度集合应尽量使反射方向仍在集合中；否则需要角度插值，并检查插值是否保持入射、出射电流平衡。反射条件把不同方向耦合起来，可以在角度外迭代中使用上一轮出射通量，也可以组装方向间矩阵。

白反射把返回中子在入射半球上各向同性重新分布。更一般的“镜面—白反射”模型可写成

$$
\boxed{
\psi_b(\boldsymbol\Omega)
=\beta\left[(1-\alpha_w)\psi_b(\boldsymbol\Omega^*)
+\frac{\alpha_w}{\pi}
\int_{\boldsymbol\Omega'\cdot\mathbf n>0}
(\boldsymbol\Omega'\cdot\mathbf n)\psi_b(\boldsymbol\Omega')
\,\mathrm d\Omega'\right]
},
\tag{5.293}
$$

其中 $\boldsymbol\Omega\cdot\mathbf n<0$，$\alpha_w$ 是白反射比例。因为平面半球的 $\int_{\Omega\cdot n<0}|\Omega\cdot n|\,\mathrm d\Omega=\pi$，该归一化保证 $\beta=1$ 时入射与出射积分电流相等。离散实现应以同样的方向权重计算分子和分母，而不是直接对所有出射方向做不带权平均。

周期性边界把一侧入射通量映射到配对边界的出射通量，适用于无限晶格的重复单元；它不是一般的零梯度条件。无论是反射还是周期映射，都应检查配对位置、法向、方向权重以及离散守恒。

### 5.27.2 扩散真空、返照和外推边界

扩散模型不再保留入射角通量，因而要用部分电流构造近似边界。由 $P_1$ 展开式（5.273）可得向外和向内的部分电流

$$
\boxed{
J^+=\int_{\boldsymbol\Omega\cdot\mathbf n>0}
(\boldsymbol\Omega\cdot\mathbf n)\psi\,\mathrm d\Omega
=\frac{\phi_b}{4}+\frac{J_n}{2},
\qquad
J^-=\int_{\boldsymbol\Omega\cdot\mathbf n<0}
|\boldsymbol\Omega\cdot\mathbf n|\psi\,\mathrm d\Omega
=\frac{\phi_b}{4}-\frac{J_n}{2}
},
\tag{5.294}
$$

且 $J_n=J^+-J^-=-D\,\partial\phi/\partial n$。Marshak 真空条件令 $J^-=0$，所以 $J_n=\phi_b/2$。若边界返照系数定义为 $\beta=J^-/J^+$，则

$$
\boxed{
J_n=\gamma\phi_b,
\qquad
\gamma=\frac12\frac{1-\beta}{1+\beta}
}.
\tag{5.295}
$$

$\beta=0$ 给出 Marshak 真空，$\beta=1$ 给出零净流反射。Marshak 条件来自半球积分；Mark 条件则从特定离散方向求积得到，常见 $S_2$ 形式对应 $\gamma=1/\sqrt3$。这些系数取决于角度近似和边界模型，不应把不同文献中的 Robin 系数脱离定义直接互换。

工程扩散计算还常令通量在物理边界外的外推距离 $d_{\mathrm{ex}}$ 处为零。各向同性散射平面的输运理论给出近似 $d_{\mathrm{ex}}\simeq0.7104\lambda_{\mathrm{tr}}$，而 Marshak 近似对应 $2D=2\lambda_{\mathrm{tr}}/3$。这两种处理不是同一个边界条件。直接在真实物理表面施加 $\phi_b=0$ 会过强地抑制出射通量，尤其不适合作为一般真空输运边界的默认近似。[Noh 的边界条件比较](https://doi.org/10.7733/jnfcwt.2018.16.2.203)给出了 Marshak、Mark、零通量和返照条件的推导及数值差异。

### 5.27.3 Robin 边界的有限体积消元

设边界单元中心到边界的距离为 $d_{Pb}$，两点扩散通量为

$$
J_n=\frac{D_P}{d_{Pb}}(\phi_P-\phi_b).
$$

与齐次 Robin 条件 $J_n=\gamma\phi_b$ 联立，可消去未知的边界通量：

$$
\boxed{
\phi_b=\frac{D_P\phi_P}{D_P+\gamma d_{Pb}},
\qquad
J_n=G_b\phi_P,
\qquad
G_b=\frac{D_P\gamma}{D_P+\gamma d_{Pb}}
}.
\tag{5.296}
$$

因此只需在式（5.288）的 $G^{\mathrm{ext}}_{g,P}$ 中加入 $A_bG_b$。真空 Marshak 条件取 $\gamma=1/2$，理想反射取 $\gamma=0$。若直接给定 Dirichlet 值 $\phi_b=\bar\phi_b$，则边界对 $a_P$ 的贡献为 $A_bD_P/d_{Pb}$，右端增加 $A_bD_P\bar\phi_b/d_{Pb}$。

若 Marshak 边界给定独立入射部分电流 $\bar J^-$，由式（5.294）有 $J_n=\phi_b/2-\bar J^-$。与两点扩散通量联立得

$$
\boxed{
J_n=G_b\phi_P-H_b\bar J^-,
\qquad
G_b=\frac{D_P}{2D_P+d_{Pb}},
\qquad
H_b=\frac{2D_P}{2D_P+d_{Pb}}
},
\tag{5.297}
$$

这里 $D_P$ 的量纲为长度，$G_b$ 和 $H_b$ 均无量纲，而 $\bar J^-$ 与 $J_n$ 的量纲相同。代入式（5.288）时，$A_bG_b$ 进入 $G^{\mathrm{ext}}_{g,P}$；把 $J_n$ 留在控制体平衡左端时，独立入射源使右端增加 $A_bH_b\bar J^-$。一般 Robin 入射模型应重新推导相应系数，不能不加区分地沿用这组 Marshak 系数。

材料内部界面则要求标量通量和法向净流连续，式（5.287）已经自动实现这两点。若采用非匹配网格、间断有限元或不同物理模型耦合，应通过共同界面通量或 mortar 条件保证守恒。扩散—输运耦合界面还需要保持角矩或部分电流的一致映射，不能把标量通量相等当作唯一条件。

### 5.27.4 边界检查的守恒标准

边界条件实现后，至少进行以下检查：

1. 无裂变和外源、无吸收、全反射且截面恒定时，总中子数应保持不变；
2. 真空边界只有入射方向被置零，出射方向仍可携带泄漏；
3. 稳态固定源问题中，外源与裂变产生率及入射输入之和应等于吸收与出射泄漏之和；
4. 交换单元属主与邻接单元后，内部界面通量只改变符号；
5. 扩散系数跨越几个数量级时，界面通量仍连续且解不依赖单元中心距离的人为选取；
6. 全反射本征值问题若缺少吸收和泄漏，可能出现零空间或不唯一性，应检查物理模型和矩阵可解性。

## 5.28 线性系统、源迭代与 $k_{\mathrm{eff}}$ 求解

### 5.28.1 固定源方程和输运扫描

冻结散射、裂变和边界耦合源后，每个群—方向或每个扩散群都可写成稀疏线性系统

$$
\boxed{
\mathsf T_{g,m}\boldsymbol\psi_{g,m}=\mathbf q_{g,m}
\quad\text{或}\quad
\mathsf T_g\boldsymbol\phi_g=\mathbf q_g
}.
\tag{5.298}
$$

$\mathsf T_{g,m}$ 包含对流、总碰撞和时间存储，$\mathsf T_g$ 包含扩散、移除和时间存储。求解方法可按问题结构选择：

- **显式输运扫描**：适合无循环依赖的直角网格 $S_N$ 对流系统，按上游顺序逐单元求解，内存开销小；
- **稀疏直接求解**：适合中小规模多群扩散或需要反复使用同一矩阵的情况，可用 LU/Cholesky 分解并复用因子；
- **Krylov 迭代**：非对称迎风或群耦合系统可用 GMRES、BiCGStab；对称正定扩散系统可用共轭梯度；
- **多重网格预条件**：椭圆扩散系统可用几何或代数多重网格，Krylov 方法负责处理非对称和群间耦合。

内迭代容差应明显小于外层源迭代或本征值迭代的当前误差，否则代数误差会限制物理残差。反过来，在外迭代刚开始时把每个线性系统解到机器精度通常浪费时间，可使用随外残差下降而收紧的容差策略，并在最后检查真实方程残差。

### 5.28.2 散射源迭代与多群迭代

把式（5.271）的散射和裂变源放在上一轮，可写成

$$
\boxed{
\mathsf T\boldsymbol\Psi^{(r+1/2)}
=\mathsf S\boldsymbol\Psi^{(r)}
+\frac1{k_{\mathrm{eff}}}\mathsf F\boldsymbol\Psi^{(r)}
+\mathbf q
}.
\tag{5.299}
$$

每次用扫描或线性求解得到新的角通量，再更新散射矩和裂变源，称为源迭代。纯散射单群问题中，其收敛速度通常随散射比 $c=\Sigma_s/\Sigma_t$ 接近 1 而恶化；在光学厚、强散射介质中，未加速源迭代可能极慢。

多群扩散也可按“快群到热群”的顺序求解。若只有向下散射，群方程具有单向依赖；存在向上散射时，需要在群间反复迭代或求解完整块矩阵。群内迭代、群间散射迭代、角通量源迭代和裂变外迭代应在日志中区分，否则无法判断收敛瓶颈。

### 5.28.3 $k$ 本征值的幂迭代

幂迭代是确定性临界计算最常用的外层算法。给定正的非零初始通量和 $k^{(0)}$，第 $l$ 次迭代求解

$$
\boxed{
\mathsf A\mathbf x^{(l+1)}
=\frac1{k^{(l)}}\mathsf F\mathbf x^{(l)}
}.
\tag{5.300}
$$

令 $\mathcal P_l$ 表示与 $\mathbf x^{(l)}$ 对应的体积积分裂变中子产生率，则更新公式为

$$
\boxed{
k^{(l+1)}=k^{(l)}\frac{\mathcal P_{l+1}}{\mathcal P_l},
\qquad
\mathcal P_l=\sum_{P,g}V_P\nu\Sigma_{f,g,P}\phi_{g,P}^{(l)}
}.
\tag{5.301}
$$

也可以在角度—能量—空间积分意义下写 $\mathcal P_l=\mathbf w^T\mathsf F\mathbf x^{(l)}$。此时 $\mathbf w$ 必须包含体积、方向权重并正确处理 $\chi_g$ 的归一化；对完整 $S_N$ 未知量直接重复累加裂变源而不做方向积分，会多出 $4\pi$ 因子。每次更新后可按 $\mathcal P_{l+1}$ 或总功率重新归一化通量，避免数值过大或过小，归一化不改变本征值。

典型收敛指标包括

$$
\boxed{
\eta_k=\frac{|k^{(l+1)}-k^{(l)}|}{k^{(l+1)}},
\qquad
\eta_R=\frac{\left\|\mathsf A\mathbf x-k^{-1}\mathsf F\mathbf x\right\|_2}
{\left\|k^{-1}\mathsf F\mathbf x\right\|_2+\varepsilon},
\qquad
\eta_x=\left\|\widetilde{\mathbf x}^{(l+1)}-\widetilde{\mathbf x}^{(l)}\right\|_\infty
},
\tag{5.302}
$$

其中 $\widetilde{\mathbf x}$ 是同一规则归一化后的通量。只看 $k$ 的稳定并不充分：在优势比接近 1 时，$k$ 可能先趋于稳定，而裂变源空间分布仍在缓慢漂移。应同时检查本征残差和通量/功率分布。

若 $\mathsf A^{-1}\mathsf F$ 的最大两个特征值模满足 $|\lambda_1|>|\lambda_2|$，幂迭代误差渐近按 $|\lambda_2/\lambda_1|^l$ 衰减。该比值称为优势比，越接近 1 收敛越慢。可采用 Chebyshev 加速、Wielandt 移位、Anderson 加速、Krylov 子空间本征求解器或 Jacobian-free Newton–Krylov，但必须保持基波正通量和全局裂变平衡。

### 5.28.4 扩散合成加速与一致离散

扩散合成加速（DSA）用低阶扩散问题修正输运源迭代的长波误差。其抽象步骤为

$$
\boxed{
\mathsf A_D\,\delta\boldsymbol\phi
=\mathsf S_D\left(\boldsymbol\phi^{(r+1/2)}-\boldsymbol\phi^{(r)}\right),
\qquad
\boldsymbol\Psi^{(r+1)}
=\boldsymbol\Psi^{(r+1/2)}+\mathsf L^\dagger\delta\boldsymbol\phi
},
\tag{5.303}
$$

$\mathsf A_D$ 是低阶扩散算子，$\mathsf L^\dagger$ 把标量修正提升为角通量或角矩修正。DSA 只加速迭代，不应改变最终的高阶离散方程。

低阶算子与高阶对流离散之间的一致性十分重要。直接套用连续 $P_1$ 扩散系数，而高阶部分采用菱形差分、间断有限元或特殊边界处理，可能导致加速失效甚至发散。经典理论见 [Alcouffe 的 DSA 方法](https://doi.org/10.13182/NSE77-1)和 [Adams–Larsen 的迭代方法综述](https://doi.org/10.1016/S0149-1970(01)00023-3)。实现时应分别验证固定源散射迭代、DSA 修正、反射边界和粗网格强散射极限。

### 5.28.5 瞬态通量与缓发先驱核耦合

把多群扩散或输运方程写成紧凑形式，后向 Euler 给出

$$
\boxed{
\left(\frac{\mathsf M_v}{\Delta t}+\mathsf A\right)\boldsymbol\Phi^{n+1}
=\frac{\mathsf M_v}{\Delta t}\boldsymbol\Phi^n
+\mathsf F_p\boldsymbol\Phi^{n+1}
+\sum_i\mathsf D_i\mathbf C_i^{n+1}
+\mathbf Q^{n+1}
},
\tag{5.304}
$$

其中 $\mathsf M_v$ 存储 $1/v_g$，$\mathsf F_p$ 为瞬发裂变算子，$\mathsf D_i$ 包含 $\chi_{d,i}\lambda_i$。先驱核方程可写为

$$
\boxed{
\left(\frac1{\Delta t}+\lambda_i\right)\mathbf C_i^{n+1}
=\frac1{\Delta t}\mathbf C_i^n
+\mathsf F_{d,i}\boldsymbol\Phi^{n+1}
}.
\tag{5.305}
$$

可把 $\mathbf C_i^{n+1}$ 消元代入式（5.304），得到带历史源和附加对角项的通量系统；也可在每个时间步内进行通量—先驱核耦合迭代。初始稳态应满足 $\mathbf C_i=\mathsf F_{d,i}\boldsymbol\Phi/\lambda_i$，否则会在无扰动时产生人为瞬态。

隐式格式稳定并不意味着大步长准确。快速反应性插入、控制棒移动和反馈变化需要总功率、峰值功率和先驱核浓度的时间步收敛试验。Crank–Nicolson、二阶后差或指数时间积分可提高精度，但应分别检查振荡、非负性和刚性源处理。完整空间—时间求解与点堆/改进准静态近似也不能混用同一组动力学参数而不核对定义。

## 5.29 通量、$k_{\mathrm{eff}}$ 与其他重要物理量

### 5.29.1 反应性、反应率和功率归一化

求得基波本征值后，反应性定义为

$$
\boxed{
\rho=\frac{k_{\mathrm{eff}}-1}{k_{\mathrm{eff}}}=1-\frac1{k_{\mathrm{eff}}}
}.
\tag{5.306}
$$

$\rho>0$、$=0$ 和 $<0$ 分别对应超临界、临界和次临界。$1\ \mathrm{pcm}=10^{-5}$ 的反应性单位；以 $\beta_{\mathrm{eff}}$ 为单位时，$\rho[\$]=\rho/\beta_{\mathrm{eff}}$。文献也常把 $10^5(k-k_{\mathrm{ref}})$ 简称为“pcm 本征值差”，它与严格的反应性差 $10^5(1/k_{\mathrm{ref}}-1/k)$ 略有区别，报告时应注明采用哪一种。

区域 $V_j$ 内反应 $x$ 的反应率和由裂变能量释放决定的功率为

$$
\boxed{
R_{x,j}=\int_{V_j}\!\int_0^\infty\Sigma_x(\mathbf r,E)
\phi(\mathbf r,E)\,\mathrm dE\,\mathrm dV,
\qquad
P_{\mathrm{th}}=\int_V\!\int_0^\infty
\kappa(E)\Sigma_f(\mathbf r,E)\phi(\mathbf r,E)\,\mathrm dE\,\mathrm dV
},
\tag{5.307}
$$

$\kappa$ 是每次裂变可沉积的能量，单位可为 J；若核数据库直接提供 $\kappa\Sigma_f$ 或发热截面，应优先保持同一数据来源和群结构。裂变功率不是所有中子吸收发热，也不应把 $\nu\Sigma_f$ 错当成 $\Sigma_f$。

本征函数幅值是任意的。给定目标总功率 $P_{\mathrm{target}}$ 后，通量归一化因子为

$$
\boxed{
c_P=\frac{P_{\mathrm{target}}}{P_{\mathrm{th}}},
\qquad
\boldsymbol\Phi\leftarrow c_P\boldsymbol\Phi,
\qquad
P_j=\frac{P_{\mathrm{target}}P_{\mathrm{th},j}}{P_{\mathrm{th}}}
}.
\tag{5.308}
$$

$k_{\mathrm{eff}}$ 和相对通量形状不因该缩放改变。若只需平均值为 1 的棒功率，可使用 $p_j=N_{\mathrm{pin}}P_{\mathrm{th},j}/\sum_{\ell=1}^{N_{\mathrm{pin}}}P_{\mathrm{th},\ell}$。燃耗、热工水力或探测器响应若依赖绝对通量，必须在后处理前完成一致的功率归一化。

### 5.29.2 泄漏平衡、功率误差和收敛判据

对收敛后的无外源本征解，可独立计算裂变产生率 $\mathcal P_f$、吸收率 $\mathcal A$ 和外表面净流出率 $\mathcal L$，并检查

$$
\boxed{
\varepsilon_B=\frac{\mathcal P_f}{k_{\mathrm{eff}}}-\mathcal A-\mathcal L,
\qquad
k_{\mathrm{balance}}=\frac{\mathcal P_f}{\mathcal A+\mathcal L}
}.
\tag{5.309}
$$

$\varepsilon_B$ 相对于 $\mathcal P_f/k_{\mathrm{eff}}$ 的大小以及 $k_{\mathrm{balance}}$ 与幂迭代本征值的差异，能发现边界通量符号、角度权重、散射矩阵方向和漏装源项等错误。固定外源问题还应把外源率和独立入射率写入平衡，不能把此时得到的比值直接当作 $k_{\mathrm{eff}}$。

C5G7 等基准常同时报告最大/最小棒功率、组件功率以及全部分布的误差。令第 $j$ 根燃料棒的相对百分误差为 $e_j=100(p_j-p_j^{\mathrm{ref}})/p_j^{\mathrm{ref}}$，可定义

$$
\boxed{
\mathrm{AVG}=\frac1N\sum_{j=1}^{N}|e_j|,
\qquad
\mathrm{RMS}=\sqrt{\frac1N\sum_{j=1}^{N}e_j^2},
\qquad
\mathrm{MRE}=100\frac{\sum_j|p_j-p_j^{\mathrm{ref}}|}
{\sum_jp_j^{\mathrm{ref}}}
}.
\tag{5.310}
$$

AVG 和 RMS 给每根棒相同权重，MRE 按参考功率加权总偏差。Monte Carlo 参考值还含统计不确定度，确定性结果与参考值的差不能简单全部解释为确定性离散误差；比较时应同时给出参考值版本、置信区间和本方网格/角度/群结构。

### 5.29.3 伴随通量、有效缓发份额和中子代时间

伴随通量描述不同位置和能量的中子对链式反应或目标响应的重要性。与式（5.278）对应的离散伴随本征问题为

$$
\boxed{
\mathsf A^T\boldsymbol\Phi^\dagger
=\frac1{k_{\mathrm{eff}}}\mathsf F^T\boldsymbol\Phi^\dagger
}.
\tag{5.311}
$$

对非对称 FVM/$S_N$ 离散，直接转置实际离散矩阵通常比“仅把连续方程中的对流项改号再重新离散”更可靠，因为后者还需保证数值边界和权重构成真正的代数伴随。

令 $\mathsf F_{d,i}$ 为第 $i$ 组缓发中子产生算子，$\mathsf M_v$ 为 $1/v$ 加权矩阵，常用的伴随加权动力学参数为

$$
\boxed{
\beta_{\mathrm{eff},i}
=\frac{(\boldsymbol\Phi^\dagger)^T\mathsf F_{d,i}\boldsymbol\Phi}
{(\boldsymbol\Phi^\dagger)^T\mathsf F\boldsymbol\Phi},
\qquad
\beta_{\mathrm{eff}}=\sum_i\beta_{\mathrm{eff},i},
\qquad
\Lambda_{\mathrm{eff}}
=\frac{(\boldsymbol\Phi^\dagger)^T\mathsf M_v\boldsymbol\Phi}
{(\boldsymbol\Phi^\dagger)^T\mathsf F\boldsymbol\Phi}
}.
\tag{5.312}
$$

这里的内积必须包含与正问题一致的体积、能量群和方向权重。$\beta_{\mathrm{eff}}$ 不同于未经伴随加权的总缓发份额；反射层、不同能谱和局部吸收会显著改变中子重要性。$\Lambda_{\mathrm{eff}}$ 是点堆方程约定下的中子代时间，常见关系为瞬发寿命 $\ell=k_{\mathrm{eff}}\Lambda_{\mathrm{eff}}$。不同程序可能把 $k$ 因子放在不同的 tally 定义中，比较输出时应核对实际公式，而不只核对参数名称。[OpenMC 动力学参数文档](https://docs.openmc.org/en/stable/usersguide/kinetics.html)说明了迭代裂变概率法及其实现约定。

小扰动下，若 $\mathsf A$ 和 $\mathsf F$ 分别变为 $\mathsf A+\delta\mathsf A$、$\mathsf F+\delta\mathsf F$，一阶微扰给出

$$
\boxed{
\frac{\delta k}{k}
\approx
\frac{(\boldsymbol\Phi^\dagger)^T
(\delta\mathsf F-k\,\delta\mathsf A)\boldsymbol\Phi}
{(\boldsymbol\Phi^\dagger)^T\mathsf F\boldsymbol\Phi}
}.
\tag{5.313}
$$

该式可用于材料截面微扰和控制棒价值估计。大扰动或局部强扰动会改变通量形状，一阶近似不再可靠，应重新求解扰动后的本征问题或使用精确广义微扰方法。

### 5.29.4 离散收敛与不确定度

设同一物理量在三级等比加密网格上的值为 $Q_h,Q_{h/2},Q_{h/4}$，可用

$$
\boxed{
p_{\mathrm{obs}}\approx
\frac{\ln\left|(Q_h-Q_{h/2})/(Q_{h/2}-Q_{h/4})\right|}{\ln2}
}
\tag{5.314}
$$

估计渐近收敛阶。该公式要求三级解都进入渐近区，且内迭代误差远小于相邻网格解之差。$S_N$ 计算还要独立细化方向数；瞬态还要独立减小 $\Delta t$。空间、角度、能量群、时间和核数据不确定度不能通过一次“网格无关性”试验全部排除。

推荐把误差分成四层报告：线性求解误差、外迭代误差、离散误差、模型与核数据误差。若采用 Monte Carlo 参考解，还应增加统计误差，并检查其裂变源是否已收敛。[OpenMC 本征值方法文档](https://docs.openmc.org/en/stable/methods/eigenvalue.html)特别强调 $k$ 与裂变源分布的收敛并不等价，这一点同样适用于确定性幂迭代。

## 5.30 经典数值算例

### 5.30.1 解析解、制造解与 Reed 非均匀平板

在编写复杂堆芯求解器之前，应先做以下基础验证：

1. **纯吸收无限介质**：给定入射通量，检验 $\psi(s)=\psi_0e^{-\Sigma_ts}$，检查真空入流、光学厚度和指数衰减；
2. **均匀扩散平板或立方体**：与解析特征值和正弦通量比较，检查扩散系数、移除项和外推边界；
3. **制造解方法（MMS）**：先选定光滑通量，再代回方程生成散射、裂变和外源项，以同时验证空间算子、群耦合和边界实现；
4. **Reed 问题**：使用强吸收、强散射、真空和多材料分区的一维平板，检验空间格式在非光滑解和不同光学厚度区域中的鲁棒性。

Reed 问题源自 W. H. Reed 关于输运差分格式的工作（*Nuclear Science and Engineering*, 1971, [DOI: 10.13182/NSE46-309](https://doi.org/10.13182/NSE46-309)）。它常被用来比较菱形差分、指数格式、迎风格式和高阶 $S_N$ 方法。应使用同一文献版本的材料分区、源强和参考点通量；Reed 问题的参考解本身也可能是高精度数值解，不能笼统称为解析解。

### 5.30.2 IAEA 2D/3D PWR 与 Takeda 三维堆芯

**IAEA 2D PWR** 是两群、多材料区的经典稳态扩散算例，含控制棒和水反射层造成的强局部扰动。它适合检验：

- 多群移除、向下散射和裂变矩阵的组装；
- 材料界面的调和扩散导纳；
- 幂迭代和组件功率归一化；
- 细网 FDM/FVM 与粗网高阶或节点法的误差。

常用细网参考值为 $k_{\mathrm{eff}}=1.029585$，可参见 [Kromar 的 IAEA 基准计算报告](https://inis.iaea.org/records/xpmrb-n9d92/files/25044526.pdf)。该数值属于特定二维两群扩散模型，不是连续能量输运的“真实值”。三维 IAEA、Biblis 等相关算例的几何和边界不能由二维名称推断，必须使用对应说明书。

**Takeda 基准**包含小型三维非均匀堆芯及控制棒工况，适合比较三维扩散、$P_N$、$S_N$、节点输运和 Monte Carlo 的 $k_{\mathrm{eff}}$、控制棒价值与区域通量。其核心价值是揭示三维扩散近似与输运近似之间的差别，而不是只验证二维模板向三维复制是否正确。问题书可追溯到 Takeda 等的 NEACRP-L-330；[三维输运结果汇总](https://inis.iaea.org/records/32nwz-06137/files/22085401.pdf?download=1)给出了不同程序和方法的比较。

[Argonne Benchmark Problem Book](https://rsicc.ornl.gov/codes/mis/mis0/mis-004.html)及 ANL-7416 各补编是 IAEA 以及后文 TWIGL、LRA 等经典问题的重要原始来源。复现时应引用具体题号和补编，而不是只写“IAEA 题”或“LRA 题”；Takeda 问题则应引用其独立的 NEACRP 报告。

### 5.30.3 TWIGL 与 LRA 时空动力学问题

**TWIGL 二维种子—包层堆问题**采用两能群和缓发先驱核模型，通过角部种子区热吸收截面的阶跃或斜坡变化激发瞬态。经典斜坡在 $0.2\ \mathrm{s}$ 内完成，常用比较量包括总功率随时间的变化、区域功率和最终渐近形状。它主要检验：

- 初始稳态通量与先驱核平衡；
- 阶跃/斜坡材料数据的准确时间位置；
- 后向 Euler、ADI 或其他隐式格式的时间收敛；
- 每个时间步内空间迭代的容差。

相关节点法研究还表明，零通量边界与真空边界、不同横向泄漏近似以及网格处理会给出不同参考结果。因此，文献中的 TWIGL 初始 $k$ 和功率曲线只有在边界、对称性和扰动定义完全相同时才可比较。

**LRA 二维 BWR 动力学问题**进一步包含控制棒移动、快速功率变化和温度反馈，适合验证通量—先驱核—热反馈耦合。比较内容应包括功率峰值及其到达时刻、温度历程、空间功率/温度分布和时间步内迭代次数。只比较最终功率会遗漏瞬态峰值误差。TWIGL 与 LRA 的规范及大量 FDM、节点法结果均可在上述 ANL 问题书中查找。

### 5.30.4 Kobayashi 真空通道与 dog-leg 问题

Kobayashi 三维固定源问题包含三类几何：带方形空腔的屏蔽体、直空腔通道和折转的 dog-leg 通道。每类都有纯吸收和 50% 散射工况；纯吸收可用数值积分获得高精度参考，散射工况使用 Monte Carlo 参考。[OECD/NEA 原始报告](https://www.oecd-nea.org/upload/docs/application/pdf/2020-01/nsc-doc2000-4.pdf)给出了完整几何、截面、参考点和边界定义。

这组问题特别适合检验：

- 大真空区中的自由飞行和远距离定向泄漏；
- $S_N$ 射线效应、方向求积阶数和不同求积族；
- 强材料—真空界面处的负通量、振荡和角度收敛；
- 固定源扫描、DSA 和三维并行实现；
- 扩散模型的适用范围。

原始规范为兼容部分二阶方法，把“真空”区的总截面设为 $10^{-4}\ \mathrm{cm^{-1}}$，而不是严格零。复现时必须保留这一数值。Kobayashi 报告和后续 Ardra 计算均显示，低阶方向离散可能产生射线效应，而扩散解在空腔通道内会出现根本性错误。比较量应是规定位置和截面的标量通量，而不是把该固定源问题错误地整理成 $k_{\mathrm{eff}}$ 算例。

### 5.30.5 C5G7 与 C5G7-TD 高保真输运问题

**C5G7** 是 OECD/NEA 的七群 UO$_2$/MOX 燃料组件输运基准，目标是检验**燃料棒以上不做空间均匀化**的确定性方法。原始模型采用 $17\times17$ 棒格、四种组件构成的四分之一堆芯、水反射层和七群输运修正截面，要求输出本征值、棒功率和组件功率。主要难点是强非均匀性、尖锐通量梯度和高空间—角度分辨率。

[NEA/NSC/DOC(2003)16](https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2003-16.pdf)给出了原始规格和多群 Monte Carlo 参考解。该报告的表 3 和表 10 分别列出二维 $k_{\mathrm{eff}}=1.186550$、原始三维 $k_{\mathrm{eff}}=1.183810$；表中同时报告的是统计百分误差/置信区间，而不是把这些数字视为无限精度常数。后续文献和三维扩展控制棒问题会使用不同几何、参考计算和数值，不应与上述两个值混用。

建议按以下顺序复现 C5G7：

1. 核对七群截面单位、散射矩阵方向、$\chi_g$ 和材料编号；
2. 在单组件或较粗模型上验证固定源与 $k$ 幂迭代；
3. 分别进行网格和角度收敛试验，再比较 $k$、最大/最小棒功率和组件功率；
4. 使用式（5.310）报告 AVG、RMS、MRE，并注明参考解统计误差；
5. 最后才扩大到全三维、并行加速和多物理耦合。

**C5G7-TD** 在此基础上检验无空间均匀化的时间相关输运，包括 Cartesian/六角几何、中子噪声以及后续带反馈的动力学阶段。[OECD/NEA C5G7-TD 官方页面](https://oecd-nea.org/jcms/pl_32145/deterministic-time-dependent-neutron-transport-benchmark-without-spatial-homogenisation-c5g7-td)说明了各阶段和数据包。复现瞬态时必须指定所选 Phase/Part、数据版本和参考解形成方式；该基准仍在分阶段发展，不能用静态 C5G7 的结果代替瞬态验证。

### 5.30.6 算例选择对照

| 算例 | 主要模型 | 首要验证目标 | 容易误判之处 |
|---|---|---|---|
| 均匀平板/制造解 | 扩散或单群输运 | 系数、源项、边界和收敛阶 | 把代数收敛当成离散收敛 |
| Reed | 一维多材料固定源输运 | 光学厚度、非光滑解、空间格式 | 混用不同材料分区或参考解版本 |
| IAEA 2D/3D | 多群扩散本征值 | $k$、组件功率、材料界面 | 把扩散参考值当作输运真值 |
| Takeda | 三维多群堆芯 | 三维输运、控制棒价值、伴随 | 只比较 $k$，不比较区域通量 |
| TWIGL | 两群时空扩散动力学 | 阶跃/斜坡、先驱核、时间推进 | 忽略边界定义和时间步内残差 |
| LRA | 两群动力学与反馈 | 功率峰值、热反馈、耦合迭代 | 只比较终态，不比较峰值时刻 |
| Kobayashi | 三维固定源输运 | 真空通道、射线效应、方向收敛 | 使用扩散模型或把空腔严格置零 |
| C5G7 | 七群棒分辨输运 | $k$、棒功率、空间—角度误差 | 混用原始三维与三维扩展参考值 |
| C5G7-TD | 非均匀瞬态输运 | 时间—空间—角度耦合 | 未说明 Phase/Part 和数据版本 |

## 5.31 实施流程与常见错误

一个可复现的多群 FVM 求解流程可整理为：

```text
读取几何、材料、群常数、缓发数据和边界说明
检查长度/时间/能量单位、群索引、chi 归一化和散射矩阵方向
构造单元—界面连接、体积、面积向量、材料界面和上游依赖图
选择扩散或 SN 模型、角度求积、空间格式和时间格式
组装对流/扩散/移除矩阵，施加真空、反射、入流或 Robin 边界
固定源：进行群内求解和散射源迭代，必要时使用 DSA
本征值：进行幂迭代，更新 k 并检查通量形状和本征残差
瞬态：每个时间步耦合瞬发/缓发源及先驱核，检查时间层残差
重新计算反应率、泄漏和功率，执行独立全局平衡
按目标功率归一化，输出通量、k、功率、反应性和动力学参数
分别完成线性、外迭代、空间、角度、时间和基准题验证
```

常见实现错误包括：

- 把 $\phi=\int\psi\,\mathrm d\Omega$ 实现成不带 $w_m$ 的方向求和，或在 $S_N$ 裂变源中遗漏/重复 $4\pi$；
- 把 $\Sigma_{s,g'\to g}$ 的下标方向写反，使能量从高能群向低能群的转移错误；
- 用 $\Sigma_a$ 代替多群移除截面 $\Sigma_{r,g}$，从而漏掉群内外的散射平衡；
- 把瞬发产额、总产额和缓发份额重复计入裂变源；
- 在临界本征值迭代之外仍机械保留 $1/k$，导致瞬态物理被额外缩放；
- 用算术平均扩散系数处理强材料界面，或两侧界面没有使用同一个数值通量；
- 把真空边界理解为所有方向通量为零，或在真实表面上不加说明地令 $\phi_b=0$；
- 只按 $k$ 的尾数判断收敛，不检查裂变源、功率形状和真实本征残差；
- 在强散射粗网格上使用不一致 DSA，出现“越加速越发散”；
- 对负通量简单截零而不检查粒子数平衡，或把射线效应当成空间网格误差；
- 用绝对通量计算功率前忘记本征向量归一化；
- 把不同版本基准题的几何、边界、群常数和参考值拼接比较。

## 本章小结

一般中子输运方程在时间、空间、方向和能量上描述中子的自由飞行、碰撞、散射、裂变和外源。多群 $S_N$ 方法保留方向变量，有限体积迎风通量把每个方向的上游单元与下游单元守恒连接；$P_1$ 近似和 Fick 定律则把问题化为只含标量通量的多群扩散方程。规则网格上的 FDM 可得到类似模板，FEM 则通过弱形式、迎风界面或稳定化构造同类代数系统。

真空、反射、白反射、入流和返照边界必须与未知量的角度性质一致。扩散 Robin 边界可通过单元中心—边界两点关系消元，得到明确的矩阵对角项和入射源项。固定源计算依赖输运扫描、稀疏直接或 Krylov 求解；散射源迭代负责角矩和群间耦合，幂迭代负责 $k_{\mathrm{eff}}$ 与基波通量，DSA 等方法用于加速强散射问题。

本征向量的绝对幅值需通过总功率归一化。可信结果不仅要给出通量和 $k_{\mathrm{eff}}$，还应给出反应率、泄漏、功率、反应性以及必要的伴随加权动力学参数，并通过独立粒子平衡、真实残差、空间—角度—时间收敛和经典基准验证。扩散与输运、不同边界约定以及不同版本参考解不能混为一谈。

## 习题

1. 从式（5.269）出发，说明每一项的物理意义，并检查角通量、散射源和裂变源的量纲。
2. 在单一缓发份额假设下，将式（5.269）的瞬发裂变源改写为 $(1-\beta)\nu\Sigma_f$ 形式，并说明何时这种改写不适用。
3. 推导无源、无缓发修正的稳态多群 $k$ 本征值方程，并解释为什么裂变项含有 $1/k_{\mathrm{eff}}$。
4. 从 $P_1$ 展开推导零阶和一阶角矩方程，说明忽略哪一项后才得到通常的 Fick 扩散律。
5. 证明式（5.277）的移除截面与多群散射平衡一致，并讨论 $(n,2n)$ 反应需要怎样额外处理。
6. 从一维控制体平衡推导迎风系数式（5.282），验证其对角占优和非负源下的非负性。
7. 推导菱形差分式（5.284）和指数格式式（5.285），讨论粗光学厚度网格上两者的差异。
8. 由界面通量连续推导式（5.287），证明等距离时得到调和平均扩散系数。
9. 推导 Marshak 条件 $J_n=\phi_b/2$ 和一般返照系数式（5.295），说明为什么物理表面上的 $\phi_b=0$ 不是一般真空条件。
10. 从 Marshak 入射电流关系与两点扩散通量推导式（5.297），并给出矩阵对角项和右端项。
11. 证明幂迭代更新式（5.301），并说明通量归一化为什么不会改变本征值。
12. 对单群、单区模型解释为什么 $k$ 可能比空间裂变源分布更早表现收敛，并提出两个补充判据。
13. 推导后向 Euler 先驱核式（5.305）的消元形式，并给出无扰动稳态所需的初始先驱核浓度。
14. 从一阶广义本征方程微分推导式（5.313），讨论局部强扰动下一阶公式的误差来源。
15. 设计一套由纯吸收平板、IAEA 2D、Kobayashi 和 C5G7 组成的验证路线，分别为每一级列出输入、输出和失败判据。

## 延伸阅读

1. Duderstadt, J. J., and Martin, W. R. *Transport Theory*. Wiley, 1979.
2. Lewis, E. E., and Miller, W. F. *Computational Methods of Neutron Transport*. Wiley, 1984；后由 American Nuclear Society 重印。
3. Bell, G. I., and Glasstone, S. *Nuclear Reactor Theory*. Van Nostrand Reinhold, 1970.
4. Stacey, W. M. *Nuclear Reactor Physics*. Wiley, 2007.
5. Adams, M. L., and Larsen, E. W. Fast iterative methods for discrete-ordinates particle transport calculations. *Progress in Nuclear Energy*, 2002, 40(1): 3–159. [DOI](https://doi.org/10.1016/S0149-1970(01)00023-3).
6. Alcouffe, R. E. Diffusion synthetic acceleration methods for the diamond-differenced discrete-ordinates equations. *Nuclear Science and Engineering*, 1977, 64: 344–355. [DOI](https://doi.org/10.13182/NSE77-1).
7. Reed, W. H. New difference schemes for the neutron transport equation. *Nuclear Science and Engineering*, 1971, 46: 309–314. [DOI](https://doi.org/10.13182/NSE46-309).
8. Noh, T. A study on diffusion approximations to neutron transport boundary conditions. *Journal of Nuclear Fuel Cycle and Waste Technology*, 2018, 16(2): 203–209. [DOI](https://doi.org/10.7733/jnfcwt.2018.16.2.203).
9. Argonne Code Center. *Benchmark Problem Book*, ANL-7416 及 Supplements 1–3. [RSICC 数据包说明](https://rsicc.ornl.gov/codes/mis/mis0/mis-004.html)。
10. Kobayashi, K., Sugimura, N., and Nagaya, Y. 3D radiation transport benchmark problems and results for simple geometries with void region. *Progress in Nuclear Energy*, 2001, 39(2): 119–144. [DOI](https://doi.org/10.1016/S0149-1970(01)00007-5)；[OECD/NEA 报告](https://www.oecd-nea.org/upload/docs/application/pdf/2020-01/nsc-doc2000-4.pdf)。
11. OECD/NEA. *Benchmark on Deterministic Transport Calculations Without Spatial Homogenisation: A 2-D/3-D MOX Fuel Assembly Benchmark*, NEA/NSC/DOC(2003)16. [报告全文](https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nsc-doc2003-16.pdf).
12. OECD/NEA. *Deterministic Time-Dependent Neutron Transport Benchmark without Spatial Homogenisation (C5G7-TD)*. [官方说明与数据入口](https://oecd-nea.org/jcms/pl_32145/deterministic-time-dependent-neutron-transport-benchmark-without-spatial-homogenisation-c5g7-td).
13. Wang, Y., Schunert, S., and Baker, B. A. *Rattlesnake: User Manual*, INL/EXT-15-37337, 2015/2016. [INL 报告](https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_8838.pdf)。包含 $S_N$、$P_N$、扩散、边界、瞬态和多个经典算例的实现说明。
14. OpenMC. [Eigenvalue Methods](https://docs.openmc.org/en/stable/methods/eigenvalue.html)、[Tallies](https://docs.openmc.org/en/stable/methods/tallies.html) 和 [Kinetics Parameters](https://docs.openmc.org/en/stable/usersguide/kinetics.html)。用于理解参考解、反应率统计、源收敛和伴随加权量。
15. 李卫、刘晓晶、柴翔、赵鹏程. 基于有限体积法的中子扩散和输运算法开发. *核技术*, 2025, 48(1): 010603. [DOI](https://doi.org/10.11889/j.0253-3219.2025.hjs.48.240037).
16. Takeda, T., et al. *3-D Neutron Transport Benchmarks*, NEACRP-L-330, OECD/NEA, 1991；相关[三维输运程序结果汇总](https://inis.iaea.org/records/32nwz-06137/files/22085401.pdf?download=1).
