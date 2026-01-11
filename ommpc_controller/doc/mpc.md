# 1. 微分平坦推导

我来详细解释 `computeFlatInputwithHopfFibration` 函数，这是四旋翼控制中的**微分平坦映射**核心函数。

## 函数概述

该函数实现了从**期望加速度**到**姿态和角速度**的映射，这是四旋翼**微分平坦特性**的关键应用。使用**Hopf纤维化**方法避免了传统欧拉角表示的奇异性问题。

## 数学背景

### 四旋翼微分平坦性

四旋翼系统的**平坦输出**是位置和偏航角：$\sigma = [x, y, z, \psi]^T$

所有状态和控制输入都可以由平坦输出及其导数代数表示：
$$
\begin{aligned}
\text{位置} &: \mathbf{p} = [x, y, z]^T \\
\text{速度} &: \mathbf{v} = \dot{\mathbf{p}} \\
\text{加速度} &: \mathbf{a} = \ddot{\mathbf{p}} \\
\text{推力} &: f = m \|\mathbf{a} + g\mathbf{e}_3\| \\
\text{姿态} &: R = f(\mathbf{a}, \psi) \\
\text{角速度} &: \boldsymbol{\omega} = f(\mathbf{a}, \dot{\mathbf{a}}, \psi, \dot{\psi})
\end{aligned}
$$

## 参数说明

```cpp
输入:
  thr_acc    // 期望加速度向量 (包含重力补偿)
  jer        // 期望加加速度 (jerk)
  yaw        // 期望偏航角
  yawd       // 期望偏航角速度
  att_est    // 估计姿态 (奇异性时的备用值)

输出:
  att        // 计算得到的期望姿态四元数
  omg        // 计算得到的期望角速度
```

## 详细实现分析

### 1. **推力方向向量归一化** (第853-855行)

```cpp
Eigen::Vector3d abc = thr_acc.normalized();
double a = abc(0), b = abc(1), c = abc(2);
```

将加速度向量归一化得到机体Z轴方向：
$$
\mathbf{z}_B = \frac{\mathbf{a}_{\text{des}} + g\mathbf{e}_3}{\|\mathbf{a}_{\text{des}} + g\mathbf{e}_3\|} = \begin{bmatrix} a \\ b \\ c \end{bmatrix}
$$

**物理意义**：四旋翼的推力方向定义了机体坐标系的Z轴。

### 2. **推力方向时间导数** (第856-857行)

```cpp
Eigen::Vector3d abc_dot = (thr_acc.dot(thr_acc) * Eigen::MatrixXd::Identity(3, 3) 
                           - thr_acc * thr_acc.transpose()) 
                           / thr_acc.norm() / thr_acc.squaredNorm() * jer;
```

这是单位向量的时间导数公式：
$$
\dot{\mathbf{z}}_B = \frac{d}{dt}\left(\frac{\mathbf{v}}{\|\mathbf{v}\|}\right) = \frac{I - \mathbf{z}_B\mathbf{z}_B^T}{\|\mathbf{v}\|} \dot{\mathbf{v}}
$$

其中：
- $(I - \mathbf{z}_B\mathbf{z}_B^T)$ 是投影到垂直于 $\mathbf{z}_B$ 的平面
- $\dot{\mathbf{v}} = \text{jer}$ 是加加速度

**物理意义**：计算推力方向的变化率，用于后续计算角速度。

### 3. **奇异性检查** (第859行)

```cpp
if(1.0 + c > 1e-3 && thr_acc.norm() > 0.1)
```

**两个条件**：

#### 条件1: `1.0 + c > 1e-3`

检查是否接近**下方奇点** (gimbal lock)：
- $c = z_B \cdot \mathbf{e}_3$ 是Z轴在世界坐标系Z方向的投影
- $c = -1$ 表示机体完全倒置（Z轴向下）
- $1 + c \approx 0$ 时，Hopf纤维化公式会出现除零

#### 条件2: `thr_acc.norm() > 0.1`

检查加速度幅值是否足够大：
- 加速度太小时，方向不可靠
- 0.1 m/s² 约为重力的1%

**失败处理** (第868-871行)：
```cpp
else {
    std::cout << "Near singularity!!!!!" << std::endl;
    omg = omg_old;      // 使用上次的角速度
    att = att_est;      // 使用估计的姿态
}
```

### 4. **Hopf纤维化姿态计算** (第860-863行)

这是函数的核心！

#### Step 1: 计算倾斜四元数

```cpp
double norm = sqrt(2 * (1 + c));
Eigen::Quaterniond q((1 + c) / norm, -b / norm, a / norm, 0);
```

**Hopf纤维化公式**：
$$
q_{\text{tilt}} = \frac{1}{\sqrt{2(1+c)}} \begin{bmatrix} 1+c \\ -b \\ a \\ 0 \end{bmatrix}
$$

**数学原理**：这个四元数将世界Z轴 $\mathbf{e}_3 = [0,0,1]^T$ 旋转到推力方向 $\mathbf{z}_B = [a,b,c]^T$。

**验证**：
$$
q_{\text{tilt}} \cdot [0,0,1]^T \cdot q_{\text{tilt}}^* = [a,b,c]^T
$$

**为什么叫Hopf纤维化？**

Hopf纤维化描述了从3D空间到4D四元数空间的映射关系。给定期望的推力方向 $\mathbf{z}_B$，存在**无穷多个**姿态满足要求（绕 $\mathbf{z}_B$ 轴旋转的自由度）。这个公式选择了其中**标准的一个**。

#### Step 2: 计算偏航四元数

```cpp
Eigen::Quaterniond q_yaw(cos(yaw / 2), 0, 0, sin(yaw / 2));
```

绕世界Z轴旋转偏航角：
$$
q_{\text{yaw}} = \begin{bmatrix} \cos(\psi/2) \\ 0 \\ 0 \\ \sin(\psi/2) \end{bmatrix}
$$

#### Step 3: 组合得到完整姿态

```cpp
att = q * q_yaw;
```

$$
q_{\text{total}} = q_{\text{tilt}} \otimes q_{\text{yaw}}
$$

**顺序很重要**！
- 先应用倾斜（将Z轴对准推力）
- 再绕新Z轴旋转偏航角

### 5. **角速度计算** (第864-867行)

```cpp
double syaw = sin(yaw), cyaw = cos(yaw);
omg(0) = syaw * a_dot - cyaw * b_dot - (a * syaw - b * cyaw) * c_dot / (c + 1);
omg(1) = cyaw * a_dot + syaw * b_dot - (a * cyaw + b * syaw) * c_dot / (c + 1);
omg(2) = (b * a_dot - a * b_dot) / (1 + c) + yawd;
```

这些公式来自于**四元数微分方程**：
$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \boldsymbol{\omega} \end{bmatrix}
$$

结合 $\dot{\mathbf{z}}_B$ 和 $\dot{\psi}$，推导得到：

**Roll轴角速度** ($\omega_x$)：
$$
\omega_x = \sin\psi \cdot \dot{a} - \cos\psi \cdot \dot{b} - \frac{a\sin\psi - b\cos\psi}{1+c} \dot{c}
$$

**Pitch轴角速度** ($\omega_y$)：
$$
\omega_y = \cos\psi \cdot \dot{a} + \sin\psi \cdot \dot{b} - \frac{a\cos\psi + b\sin\psi}{1+c} \dot{c}
$$

**Yaw轴角速度** ($\omega_z$)：
$$
\omega_z = \frac{b\dot{a} - a\dot{b}}{1+c} + \dot{\psi}
$$

**第一项**：由推力方向变化引起的绕Z轴旋转  
**第二项**：偏航角速度指令

### 6. **角速度记忆** (第872行)

```cpp
omg_old = omg;
```

保存当前角速度，用于下次遇到奇异性时作为备用值。

## 完整数据流

```
期望加速度 thr_acc
    ↓ 归一化
推力方向 z_B = [a, b, c]
    ↓ Hopf映射
倾斜四元数 q_tilt
    ↓ 组合偏航
完整姿态 q_total = q_tilt ⊗ q_yaw
    ↓
输出姿态 att

期望加加速度 jer
    ↓ 投影微分
方向导数 [ȧ, ḃ, ċ]
    ↓ 四元数微分
角速度 ω = [ωx, ωy, ωz]
    ↓
输出角速度 omg
```

## 优势总结

### 1. **避免奇异性**
传统欧拉角在 pitch = ±90° 时会出现万向锁，Hopf纤维化在除了正下方外都很稳定。

### 2. **计算效率**
直接代数公式，无需迭代或数值优化。

### 3. **数学优雅**
四元数乘法天然符合旋转组合规则。

### 4. **鲁棒性**
有明确的奇异性检测和备用策略。

## 应用场景

- **轨迹跟踪**：已知期望位置/速度/加速度轨迹
- **路径规划**：从几何路径生成动力学可行轨迹
- **前馈控制**：为MPC提供参考输入
- **仿真验证**：生成理想控制指令

## 与MPC的关系

在MPC中，这个函数用于：
1. **设置参考轨迹**：将轨迹规划器输出的位置/速度转换为姿态/角速度
2. **线性化点**：提供动力学线性化所需的参考状态
3. **前馈项**：生成 $u_{\text{des}}$，MPC求解器计算增量 $\delta u$

这是一个精妙的数学-物理-工程结合的典范，体现了现代控制理论的优雅与实用！

# 2. `buildConstraintMatrix` 函数解释
我来详细解释 `buildConstraintMatrix` 函数，这是MPC优化问题中**约束矩阵构建**的核心函数。

## 函数概述

该函数构建线性约束矩阵 $A$，将MPC优化问题中的所有约束（动力学、初始条件、控制界）转换为标准形式：
$$
l \leq A z \leq u
$$

其中 $z$ 是决策变量向量。

## 约束类型总结

函数处理三类约束：

| 约束类型 | 数量 | 约束形式 |
|---------|------|---------|
| **初始条件** | nx = 9 | $\delta x_0 = x_{\text{init}}$ (等式) |
| **动力学** | nstep×nx = 180 | $-F_x\delta x_k - F_u\delta u_k + \delta x_{k+1} = 0$ (等式) |
| **控制界** | 2×nstep×nu = 160 | $u_{\min} \leq \delta u_k \leq u_{\max}$ (不等式) |

**总约束数**：9 + 180 + 160 = 349

## 详细实现分析

### 阶段1: 统计非零元素 (第413-468行)

#### 1.1 初始化 (第413-415行)

```cpp
A_nnz_ = nx;  // 初始条件贡献9个非零元素
```

#### 1.2 动力学约束统计 (第417-422行)

```cpp
for (int k = 0; k < nstep; ++k)
{
  A_nnz_ += Fx[k].nonZeros() + Fu[k].nonZeros() + nx;
}
```

每个时间步 $k$ 的约束方程：
$$
-F_x[k] \delta x_k - F_u[k] \delta u_k + I \delta x_{k+1} = 0
$$

贡献的非零元素：
- `Fx[k].nonZeros()`: $F_x$ 矩阵的非零元素（约27个）
- `Fu[k].nonZeros()`: $F_u$ 矩阵的非零元素（约12个）  
- `nx`: 单位矩阵 $I$ 的对角元素（9个）

**总计**：约 20×(27+12+9) = 960 个非零元素

#### 1.3 控制约束统计 (第424-425行)

```cpp
A_nnz_ += 2 * nstep * nu;  // 2×20×4 = 160
```

每个控制变量有上下界两个约束，每个约束1个非零元素。

#### 1.4 内存分配 (第427-432行)

```cpp
A_data_ = (c_float *)malloc(A_nnz_ * sizeof(c_float));     // 值
A_indices_ = (c_int *)malloc(A_nnz_ * sizeof(c_int));      // 行索引
A_indptr_ = (c_int *)malloc((total_vars_ + 1) * sizeof(c_int));  // 列指针
```

### 阶段2: 按列统计非零元素数 (第434-468行)

这是构建CSC格式的关键步骤！

```cpp
std::vector<int> col_nnz(total_vars_, 0);  // 每列的非零元素数
```

#### 2.1 初始条件 (第437-441行)

```cpp
for (int i = 0; i < nx; ++i)
{
  int col = i;
  col_nnz[col]++;  // 第i列（对应δx0的第i个分量）有1个非零元素
}
```

**约束方程**：$\delta x_0 = x_{\text{init}}$

**矩阵形式**：
$$
\begin{bmatrix} 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{bmatrix}_{9\times 9}
$$

#### 2.2 动力学约束 (第443-467行)

这是最复杂的部分！

**变量索引计算**：
```cpp
int xk_offset = k * (nx + nu);      // δxk的起始位置
int uk_offset = xk_offset + nx;     // δuk的起始位置  
int xkp1_offset = (k + 1) * (nx + nu);  // δx{k+1}的起始位置
```

**示例** (k=0)：
- `xk_offset = 0` → δx0占据索引[0,8]
- `uk_offset = 9` → δu0占据索引[9,12]
- `xkp1_offset = 13` → δx1占据索引[13,21]

**约束矩阵的三个部分**：

##### Part A: -Fx[k] (对应δxk)

```cpp
for (int j = 0; j < Fx[k].outerSize(); ++j)
{
  for (Eigen::SparseMatrix<double>::InnerIterator it(Fx[k], j); it; ++it)
  {
    int col = xk_offset + it.col();
    col_nnz[col]++;
  }
}
```

遍历 $F_x[k]$ 的每个非零元素，为对应的 $\delta x_k$ 列增加计数。

##### Part B: -Fu[k] (对应δuk)

```cpp
for (int j = 0; j < Fu[k].outerSize(); ++j)
{
  for (Eigen::SparseMatrix<double>::InnerIterator it(Fu[k], j); it; ++it)
  {
    int col = uk_offset + it.col();
    col_nnz[col]++;
  }
}
```

遍历 $F_u[k]$ 的每个非零元素，为对应的 $\delta u_k$ 列增加计数。

##### Part C: I (对应δx{k+1})

```cpp
for (int i = 0; i < nx; ++i)
{
  int col = xkp1_offset + i;
  col_nnz[col]++;
}
```

单位矩阵的对角元素，每个 $\delta x_{k+1}$ 的分量都有1个非零元素。

#### 2.3 控制约束 (第469-485行)

```cpp
for (int k = 0; k < nstep; ++k)
{
  int uk_offset = k * (nx + nu) + nx;
  
  // 下界约束: δuk ≥ u_min
  for (int i = 0; i < nu; ++i)
  {
    int col = uk_offset + i;
    col_nnz[col]++;
  }
  
  // 上界约束: δuk ≤ u_max
  for (int i = 0; i < nu; ++i)
  {
    int col = uk_offset + i;
    col_nnz[col]++;
  }
}
```

每个控制变量在其对应列有2个非零元素（上界+下界）。

#### 2.4 设置列指针 (第487-491行)

```cpp
A_indptr_[0] = 0;
for (int col = 0; col < total_vars_; ++col)
{
  A_indptr_[col + 1] = A_indptr_[col] + col_nnz[col];
}
```

**CSC格式的列指针**：`A_indptr_[col]` 表示第col列的非零元素在数据数组中的起始位置。

**示例**：
```
col_nnz = [2, 3, 1, 2, ...]
A_indptr = [0, 2, 5, 6, 8, ...]
          第0列  第1列 第2列 第3列
```

### 阶段3: 填充数据 (第493-569行)

#### 3.1 临时数据结构 (第493-496行)

```cpp
std::vector<int> col_pos(total_vars_, 0);  // 每列当前填充位置
std::vector<c_float> temp_data(A_nnz_);    // 临时存储值
std::vector<c_int> temp_indices(A_nnz_);   // 临时存储行索引
```

`col_pos` 跟踪每列已填充的非零元素数量。

#### 3.2 初始条件约束填充 (第501-508行)

```cpp
for (int i = 0; i < nx; ++i)
{
  int col = i;
  int pos = A_indptr_[col] + col_pos[col];  // 计算存储位置
  temp_data[pos] = 1.0;                      // 值为1
  temp_indices[pos] = constraint_idx++;      // 行索引
  col_pos[col]++;                            // 更新列位置
}
```

**存储逻辑**：
- `pos`：在数据数组中的绝对位置 = 列起始位置 + 列内偏移
- `temp_data[pos]`：非零元素的值
- `temp_indices[pos]`：该元素所在的行号（约束编号）

#### 3.3 动力学约束填充 (第510-553行)

```cpp
for (int k = 0; k < nstep; ++k)
{
  // 填充 -Fx[k]
  for (int j = 0; j < Fx[k].outerSize(); ++j)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Fx[k], j); it; ++it)
    {
      int col = xk_offset + it.col();
      int pos = A_indptr_[col] + col_pos[col];
      temp_data[pos] = -it.value();          // 注意负号！
      temp_indices[pos] = constraint_idx + it.row();
      col_pos[col]++;
    }
  }
  
  // 填充 -Fu[k] (类似逻辑)
  
  // 填充 I (单位矩阵)
  for (int i = 0; i < nx; ++i)
  {
    int col = xkp1_offset + i;
    int pos = A_indptr_[col] + col_pos[col];
    temp_data[pos] = 1.0;
    temp_indices[pos] = constraint_idx + i;
    col_pos[col]++;
  }
  
  constraint_idx += nx;  // 移动到下一组约束
}
```

**关键点**：
- $F_x$ 和 $F_u$ 的元素要**取负**（因为约束是 $-F_x \delta x_k - F_u \delta u_k + \delta x_{k+1} = 0$）
- 行索引 = `constraint_idx + it.row()`（约束组起始位置 + 矩阵内行号）

#### 3.4 控制约束填充 (第555-579行)

```cpp
for (int k = 0; k < nstep; ++k)
{
  int uk_offset = k * (nx + nu) + nx;
  
  // 下界约束
  for (int i = 0; i < nu; ++i)
  {
    int col = uk_offset + i;
    int pos = A_indptr_[col] + col_pos[col];
    temp_data[pos] = 1.0;
    temp_indices[pos] = constraint_idx++;
    col_pos[col]++;
  }
  
  // 上界约束（类似）
}
```

每个控制变量的约束：$\delta u_k \geq u_{\min}$ 和 $\delta u_k \leq u_{\max}$

### 阶段4: 数据拷贝 (第571-573行)

```cpp
memcpy(A_data_, temp_data.data(), A_nnz_ * sizeof(c_float));
memcpy(A_indices_, temp_indices.data(), A_nnz_ * sizeof(c_int));
```

将临时数组拷贝到OSQP要求的C风格数组。

## 约束矩阵结构可视化

对于N=2的简化案例：

```
变量顺序: [δx0, δu0, δx1, δu1, δx2]
          9维  4维  9维  4维  9维

约束矩阵A的行:
┌─────────────────────────────────────┐
│ I  0  0  0  0   │ 初始条件 (9行)      │
├─────────────────────────────────────┤
│-Fx -Fu I  0  0  │ 动力学k=0 (9行)    │
│ 0  0 -Fx-Fu I   │ 动力学k=1 (9行)    │
├─────────────────────────────────────┤
│ 0  I  0  0  0   │ 控制下界k=0 (4行)  │
│ 0  I  0  0  0   │ 控制上界k=0 (4行)  │
│ 0  0  0  I  0   │ 控制下界k=1 (4行)  │
│ 0  0  0  I  0   │ 控制上界k=1 (4行)  │
└─────────────────────────────────────┘
```

## CSC格式示例

假设某列（对应 $\delta u_0$ 的第一个分量）有3个非零元素：

```
A_indptr_[9] = 45      (该列起始位置)
A_indptr_[10] = 48     (下一列起始位置)

A_data_[45] = -0.123   (动力学约束中的值)
A_data_[46] = 1.0      (下界约束)
A_data_[47] = 1.0      (上界约束)

A_indices_[45] = 12    (第12行约束)
A_indices_[46] = 189   (第189行约束)
A_indices_[47] = 193   (第193行约束)
```

## 性能优化技巧

### 1. **两遍扫描算法**
- 第一遍：统计非零元素数量
- 第二遍：填充数据
- 避免动态内存重分配

### 2. **稀疏迭代器**
只遍历 $F_x$ 和 $F_u$ 的非零元素，跳过大量零元素。

### 3. **列优先存储**
CSC格式天然适配OSQP的列优先算法。

## 与其他组件的配合

```
setStateMatricesandBounds()
    ↓ 计算 Fx, Fu
buildConstraintMatrix()
    ↓ 构建 A
buildConstraintVectors()
    ↓ 构建 l, u
solve()
    ↓ OSQP求解: min 0.5*z'*P*z  s.t. l ≤ Az ≤ u
```

这个函数是MPC约束表示的核心，将复杂的物理约束转换为求解器可直接使用的稀疏矩阵格式，展现了数值优化与控制理论的完美融合！

已汇总会话历史记录## `setStateMatricesandBounds` 函数详解

这是 MPC 控制器的**核心动力学线性化函数**，负责在每个预测时间步构建误差状态空间模型的系数矩阵和控制约束。

---

# 3. 流形线性化函数
### 一、函数作用

```cpp
void setStateMatricesandBounds(
    const int i,                        // 预测步索引 (0 ~ nstep-1)
    const Eigen::Quaterniond &q,        // 参考姿态四元数
    const Eigen::Vector3d &omg,         // 参考角速度 (body frame)
    const double t_step,                // 离散时间步长
    const double thracc)                // 参考推力加速度
```

**输出**：
- `Fx[i]`：9×9 稀疏矩阵，状态误差转移矩阵（$\frac{\partial f}{\partial \delta x}$）
- `Fu[i]`：9×4 稀疏矩阵，控制输入矩阵（$\frac{\partial f}{\partial \delta u}$）
- `u_lb[i]`, `u_ub[i]`：控制输入的上下界约束

---

### 二、线性化动力学模型

#### 误差状态空间定义

$$
\delta x = \begin{bmatrix} \delta p \\ \delta v \\ \delta R \end{bmatrix} \in \mathbb{R}^9, \quad
\delta u = \begin{bmatrix} \delta T \\ \delta \omega \end{bmatrix} \in \mathbb{R}^4
$$

其中：
- $\delta p$：位置误差（世界系）
- $\delta v$：速度误差（世界系）
- $\delta R$：SO(3) 李代数姿态误差
- $\delta T$：推力加速度误差
- $\delta \omega$：角速度误差（机体系）

#### 离散时间线性化模型

$$
\delta x_{k+1} = F_x[k] \cdot \delta x_k + F_u[k] \cdot \delta u_k
$$

---

### 三、代码逐块解析

#### **1. 构建 $F_x$ 矩阵（状态转移矩阵）**

$F_x$ 是 9×9 分块矩阵，结构如下：

$$
F_x = \begin{bmatrix}
I_3 & \Delta t \cdot I_3 & 0_{3\times3} \\
0_{3\times3} & I_3 & \Delta t \cdot R \cdot [e_z \times T] \\
0_{3\times3} & 0_{3\times3} & \exp(-[\omega]_\times \Delta t)
\end{bmatrix}
$$

##### **(1) 位置-位置块 (0,0)：单位矩阵**
```cpp
// (0,0): I_3
for (int k = 0; k < 3; ++k) {
    tripletList.push_back(Eigen::Triplet<double>(k, k, 1.0));
}
```
**物理含义**：位置误差保持项（惯性）

---

##### **(2) 速度-速度块 (3,3)：单位矩阵**
```cpp
// (3,3): I_3
for (int k = 0; k < 3; ++k) {
    tripletList.push_back(Eigen::Triplet<double>(3+k, 3+k, 1.0));
}
```
**物理含义**：速度误差保持项

---

##### **(3) 位置-速度耦合 (0,3)：$\Delta t \cdot I_3$**
```cpp
// (0,3): Δt * I_3
for (int k = 0; k < 3; ++k) {
    tripletList.push_back(Eigen::Triplet<double>(k, 3+k, t_step));
}
```
**物理含义**：运动学积分 $p_{k+1} = p_k + v_k \cdot \Delta t$

---

##### **(4) 姿态-姿态块 (6,6)：SO(3) 指数映射**
```cpp
// (6,6): exp(-[ω]_× * Δt)
Eigen::Matrix3d exp_mat = SO3::exp(-omg * t_step);
for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
        tripletList.push_back(Eigen::Triplet<double>(6+row, 6+col, exp_mat(row, col)));
    }
}
```

**数学原理**：
- SO(3) 上的误差传播：$\delta R_{k+1} = \exp(-[\omega]_\times \Delta t) \cdot \delta R_k$
- 负号来自左乘/右乘约定

**物理含义**：姿态误差在参考角速度作用下的旋转演化

---

##### **(5) 速度-姿态耦合 (3,6)：推力方向扰动**
```cpp
// (3,6): Δt * R * [e_z × (-T)]
Eigen::Matrix3d mat_3_6 = t_step * q.toRotationMatrix() * SO3::hat(Eigen::Vector3d(0, 0, -thracc));
```

**推导**：
- 机体推力：$F_{body} = [0, 0, T]^T$
- 世界系推力：$F_{world} = R \cdot F_{body}$
- 姿态扰动影响：$\delta F \approx R \cdot [F_{body}]_\times \cdot \delta R = R \cdot [e_z \times T] \cdot \delta R$

**物理含义**：姿态误差导致推力方向偏移，进而产生速度误差

---

#### **2. 构建 $F_u$ 矩阵（控制输入矩阵）**

$F_u$ 是 9×4 分块矩阵：

$$
F_u = \begin{bmatrix}
0_{3\times4} \\
\Delta t \cdot R \cdot e_z & 0_{3\times3} \\
0_{3\times1} & J_L(\omega \Delta t)^T \Delta t
\end{bmatrix}
$$

##### **(1) 速度-推力耦合 (3,0)**
```cpp
// (3,0): Δt * R * e_z
Eigen::Vector3d vec_3_0 = t_step * q.toRotationMatrix() * Eigen::Vector3d(0, 0, 1);
for (int k = 0; k < 3; ++k) {
    tripletList.push_back(Eigen::Triplet<double>(3+k, 0, vec_3_0(k)));
}
```

**物理含义**：推力扰动 $\delta T$ 沿机体 Z 轴方向产生加速度，转换到世界系后影响速度

---

##### **(2) 姿态-角速度耦合 (6,1:3)**
```cpp
// (6,1:3): J_L(ω*Δt)^T * Δt
Eigen::Matrix3d mat_6_1 = SO3::leftJacobian(omg * t_step).transpose() * t_step;
```

**数学原理**：
- SO(3) 左雅可比矩阵：$J_L(\theta) = \frac{\sin(\|\theta\|)}{\|\theta\|} I + \frac{1 - \cos(\|\theta\|)}{\|\theta\|^2} [\theta]_\times + \frac{\|\theta\| - \sin(\|\theta\|)}{\|\theta\|^3} \theta \theta^T$
- 连接角速度扰动与姿态误差：$\delta R_{k+1} \approx J_L(\omega \Delta t)^T \Delta t \cdot \delta \omega$

**物理含义**：角速度扰动导致的姿态变化（考虑了李群流形几何）

---

#### **3. 设置控制约束**

```cpp
u_ub[i] << thracc - param_.min_thrust,              // 推力上界
            param_.max_bodyrate_xy + omg(0),        // 角速度 x 上界
            param_.max_bodyrate_xy + omg(1),        // 角速度 y 上界
            param_.max_bodyrate_z + omg(2);         // 角速度 z 上界

u_lb[i] << -(param_.max_thrust - thracc),           // 推力下界
            omg(0) - param_.max_bodyrate_xy,        // 角速度 x 下界
            omg(1) - param_.max_bodyrate_xy,        // 角速度 y 下界
            omg(2) - param_.max_bodyrate_z;         // 角速度 z 下界
```

**约束形式**：
$$
u_{lb}[i] \leq \delta u[i] \leq u_{ub}[i]
$$

**实际控制量**：
$$
u_{actual} = u_{reference} + \delta u
$$

因此：
- 推力约束：$T_{min} \leq T_{ref} + \delta T \leq T_{max}$  
  → $\delta T \in [T_{min} - T_{ref}, T_{max} - T_{ref}]$
- 角速度约束：$-\omega_{max} \leq \omega_{ref} + \delta\omega \leq \omega_{max}$  
  → $\delta\omega \in [-\omega_{max} - \omega_{ref}, \omega_{max} - \omega_{ref}]$

---

### 四、矩阵结构可视化

#### $F_x$ 矩阵（9×9，约 27 个非零元素）

```
       δpx δpy δpz δvx δvy δvz δRx δRy δRz
δpx  [  1   0   0  Δt   0   0   0   0   0 ]
δpy  [  0   1   0   0  Δt   0   0   0   0 ]
δpz  [  0   0   1   0   0  Δt   0   0   0 ]
δvx  [  0   0   0   1   0   0  c11 c12 c13]  ← (3,6) 块
δvy  [  0   0   0   0   1   0  c21 c22 c23]
δvz  [  0   0   0   0   0   1  c31 c32 c33]
δRx  [  0   0   0   0   0   0  e11 e12 e13]  ← (6,6) 块
δRy  [  0   0   0   0   0   0  e21 e22 e23]
δRz  [  0   0   0   0   0   0  e31 e32 e33]
```

#### $F_u$ 矩阵（9×4，约 12 个非零元素）

```
       δT  δωx δωy δωz
δpx  [  0   0   0   0 ]
δpy  [  0   0   0   0 ]
δpz  [  0   0   0   0 ]
δvx  [ v1   0   0   0 ]  ← (3,0) 块
δvy  [ v2   0   0   0 ]
δvz  [ v3   0   0   0 ]
δRx  [  0  j11 j12 j13]  ← (6,1:3) 块
δRy  [  0  j21 j22 j23]
δRz  [  0  j31 j32 j33]
```

---

### 五、关键设计要点

1. **Triplet 格式构建稀疏矩阵**  
   - `tripletList.reserve(27)` 预分配内存避免重复扩容
   - `setFromTriplets()` + `makeCompressed()` 高效构建 CSC 格式

2. **SO(3) 李群/李代数表示**  
   - 姿态误差用李代数 $\mathfrak{so}(3)$ 而非欧拉角，避免奇异性
   - 使用指数映射和左雅可比矩阵处理非线性旋转

3. **误差状态 MPC**  
   - 求解误差 $\delta u$ 而非绝对控制量，提高数值稳定性
   - 参考轨迹 $(q, \omega, T)$ 在外层预计算

4. **控制约束对称化**  
   - 约束相对于参考值定义，适应不同飞行状态
   - 偏航轴独立限制（通常 $\omega_{max,z} < \omega_{max,xy}$）

---

### 六、调用流程

此函数在三个参考设置函数中被调用：

1. **悬停参考** (setHoverReference)：单次调用，所有时间步使用相同矩阵
2. **文本轨迹** (setTextReference)：循环 nstep 次，每步计算不同矩阵
3. **多项式轨迹** (setTrajectoryReference)：同样循环构建时变矩阵

---

### 七、性能优化

- **稀疏矩阵**：27/81 ≈ 33% 填充率，节省存储和计算
- **预留空间**：`reserve(27)` 避免 vector 动态扩容
- **列压缩格式**：与 OSQP 求解器内部格式一致，零拷贝传递
- **每控制周期重建**：虽然矩阵稀疏，但仍需每 10ms 重新计算（时变线性化）

---

**总结**：这是 OMMPC 的"物理引擎"，将四旋翼复杂的非线性动力学在参考轨迹附近线性化为标准的状态空间形式，使凸优化求解器能够在实时约束下高效求解最优控制序列。