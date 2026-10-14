# 垂直圆翻转的姿态参考与角速度前馈改进报告

日期：2026-10-14

项目：OMMPC Controller

范围：轨迹到参考姿态、推力和角速度的转换逻辑。

文中源码链接指向仓库文件；标为“本地归档”的测试快照、日志和报告不随本仓库上传。

本次改进的核心是：沿着参考推力方向的变化连续延续姿态，并使角速度前馈与该姿态变化一致。修改主要位于参考生成层，MPC 仍根据实际状态与参考状态的误差进行优化，没有增加人工触发的 180° 翻转动作。

需要明确：原来的 Hopf 方法已经包含角速度前馈公式，多项式轨迹也已经使用 jerk。原始文本轨迹把传入的 jerk 置为零，而原来的倒置附近备用处理会破坏参考姿态与角速度的一致性。这是两个不同的问题。

## 1. 从轨迹计算推力与姿态

项目使用世界坐标系中的期望推力加速度向量：

$$
\mathbf f_d=\ddot{\mathbf p}_d+g\mathbf e_3.
$$

这里的 $\mathbf f_d$ 单位为 m/s²，相当于推力除以质量。对于推力沿机体正 Z 轴的模型，期望状态应满足：

$$
R_d\mathbf e_3=\mathbf z_d
=\frac{\mathbf f_d}{\|\mathbf f_d\|},
\qquad
\frac{T_d}{m}=\|\mathbf f_d\|.
$$

因此：

- 推力大小由 $\|\mathbf f_d\|$ 决定。
- 机体 Z 轴应对齐 $\mathbf f_d$ 的方向。
- 绕该方向如何旋转，需要额外的航向约定。

垂直圆顶点若要求推力向下，应将机体 Z 轴转到向下，并保持正的集体推力。

## 2. 原始逻辑：Hopf 映射与倒置附近的备用处理

原代码将推力方向写成：

$$
\mathbf z_d=[a,b,c]^T.
$$

然后按四元数 $(w,x,y,z)$ 顺序计算：

$$
q_H=
\frac{1}{\sqrt{2(1+c)}}
\begin{bmatrix}
1+c\\-b\\a\\0
\end{bmatrix},
\qquad
q_d=q_Hq_{\mathrm{yaw}}.
$$

正常情况下，这个构造能把机体 Z 轴转到期望推力方向，再通过右乘航向四元数选择绕推力轴的姿态。

当推力完全向下时：

$$
\mathbf z_d=[0,0,-1]^T,\qquad c=-1.
$$

公式中的 $1+c$ 为零。原角速度计算同样包含 $1/(1+c)$，因此这张方向映射在向下方向不可用。该奇异性来自选用的映射，不表示四元数或 SO(3) 无法表达倒置。

原来的处理如下：

```cpp
if (1.0 + c > 1e-3 && thr_acc.norm() > 0.1) {
    // 使用 Hopf 公式计算姿态和角速度
} else {
    omg = omg_old;
    att = att_est;
}
```

进入备用分支后，参考姿态采用传入的备用姿态，角速度采用上一次保存的值。在轨迹预测循环内部，备用姿态通常已更新为上一预测步的参考姿态，这相当于暂停参考姿态更新；角速度却可能仍然非零，二者不再描述同一个运动。

原始参考循环每次从实际估计姿态 `odom.q` 开始，因此第一个预测步进入备用分支时，也可能直接依赖实际估计姿态。

原始代码快照：Hopf 构造（本地归档：`test_results/vertical_circle_fixed_20261013/legacy/ommpc_controller.hpp`）

### 2.1 倒置附近的独立检查

之前的独立测试让推力方向连续经过倒置点。以下角度表示沿同一俯仰路径延续的角度，不是限制在主值范围内的 Euler 角。

| 期望连续俯仰角 | 原逻辑参考角 | 推力轴方向误差 |
|---:|---:|---:|
| 177° | 177° | 约 0° |
| 178° | 177° | 约 1° |
| 180° | 177° | 约 3° |
| 182° | 177° | 约 5° |
| 183° | 等效的 183° | 约 0° |

参考在倒置附近冻结，离开备用区域后重新追上期望方向。增大 MPC 权重不能修复这个参考构造过程中的不一致。

检查源代码（本地归档：`test_results/aggressive_reference_audit_20261013/reference_chart_audit.cpp`） 与 检查输出（本地归档：`test_results/aggressive_reference_audit_20261013/chart_audit.csv`） 保留在测试目录中。

## 3. 原始文本路径：计算 jerk 后仍传入零

原始文本参考通过速度差分估计加速度和 jerk：

$$
\mathbf a_i\approx
\frac{\mathbf v_i-\mathbf v_{i-1}}{\Delta t},
\qquad
\mathbf j_i\approx
\frac{\mathbf a_i-\mathbf a_{i-1}}{\Delta t}.
$$

但实际调用是：

```cpp
computeFlatInputwithHopfFibration(
    des_acc_in_world,
    Eigen::Vector3d::Zero(),  // 实际传入的 jerk 为零
    yaw, yaw_dot, last_q, q, omg);
```

原代码担心文本差分和轨迹端点不连续造成过大的 jerk 噪声，因此禁用了这部分前馈。在固定 yaw 的正常分支中，会得到：

$$
\boldsymbol\omega_{\mathrm{ref}}=0.
$$

这时参考姿态要求无人机持续翻转，角速度前馈却要求不转动。实际翻转主要依赖姿态误差出现之后的 MPC 纠偏。误差动态中的旋转演化也使用了角速度参考，所以问题不只是少了一项主动转动指令。

多项式轨迹原本传入真实 jerk；上述“传零”问题主要存在于原始文本路径。

原始文本路径调用（本地归档：`test_results/vertical_circle_fixed_20261013/legacy/ommpc_controller.hpp`）

## 4. 新姿态构造：从上一参考姿态做连续运输

新函数 `continuousFlatnessReference` 从上一参考姿态出发，用最小旋转将旧推力方向对齐新推力方向，而不是每次从固定的向上方向重新构造姿态。

设上一参考姿态为 $q_{\mathrm{prev}}$：

$$
\mathbf z_{\mathrm{prev}}
=R(q_{\mathrm{prev}})\mathbf e_3.
$$

计算将 $\mathbf z_{\mathrm{prev}}$ 转到 $\mathbf z_d$ 的最小旋转 $\Delta q$，再左乘到上一参考姿态：

$$
q_d=\Delta q\,q_{\mathrm{prev}}.
$$

对应代码：

```cpp
transport = Eigen::Quaterniond::FromTwoVectors(old_z, z);
attitude = (transport * seed).normalized();
```

当参考经过 179°、180°、181° 时，相邻推力方向之间仍然只有小角度变化。该构造没有 $1+c$ 分母，因此不会仅因推力向下而进入冻结分支。

### 4.1 参考缓存与预测窗口

原始参考每次从 `odom.q` 开始，新实现从成员变量 `reference_attitude_` 开始。这使实际跟踪滞后不会反过来改变参考的连续延续过程。

预测窗口内，每个未来参考继续使用前一个参考姿态作为起点；跨控制周期时，只缓存窗口第 0 步的参考，不缓存最后一个未来预测点。否则下一周期可能错误地从未来姿态开始。

实际姿态仍进入 MPC 的反馈误差，参考缓存并没有替代真实状态反馈。

### 4.2 四元数符号连续性

$q$ 与 $-q$ 表示同一物理姿态。新代码使用以下处理选择与上一四元数相邻的表示：

```cpp
if (attitude.dot(seed) < 0.0)
    attitude.coeffs() *= -1.0;
```

这保护的是数值表示的连续性，便于插值、差分和记录。四元数符号变化本身并不等于物理姿态跳变；原 MPC 使用旋转矩阵计算姿态误差，也不会仅因为四元数整体变号而产生不同的旋转误差。

### 4.3 相邻方向完全相反时的备用处理

若相邻两次推力方向完全相反，最小旋转轴不唯一。代码选择上一参考姿态的机体 X 轴，构造 180° 旋转。

这是确定几何解的备用处理，不表示 10 ms 内突然旋转 180° 具有动态可行性。轨迹仍应提供足够连续、足够密集的方向变化。

[连续姿态函数](../include/continuous_flatness.h) 与 [参考缓存使用位置](../include/ommpc_controller.hpp)。

## 5. 新角速度前馈：从 jerk 推导推力方向变化率

重力保持不变，因此：

$$
\dot{\mathbf f}_d=\mathbf j_d.
$$

单位推力方向的导数是：

$$
\dot{\mathbf z}_d=
\frac{(I-\mathbf z_d\mathbf z_d^T)\mathbf j_d}
{\|\mathbf f_d\|}.
$$

该式去除了 jerk 沿推力轴的分量：沿轴分量只改变推力大小，垂直轴的分量才改变推力方向。

由刚体运动学：

$$
\dot{\mathbf z}_d
=\boldsymbol\omega_W\times\mathbf z_d.
$$

这个方程不能唯一确定绕推力轴的自转。一般可写成：

$$
\boldsymbol\omega_W
=\mathbf z_d\times\dot{\mathbf z}_d
+\lambda\mathbf z_d.
$$

当前实现选择不额外绕推力轴自转，即 $\lambda=0$，并将世界系角速度转换为机体系：

$$
\boldsymbol\omega_{\mathrm{ref}}
=R_d^T\left(\mathbf z_d\times\dot{\mathbf z}_d\right).
$$

对应代码：

```cpp
z_dot = (jerk - z * z.dot(jerk)) / force.norm();
body_rate = attitude.conjugate() * z.cross(z_dot);
```

### 5.1 纯俯仰翻转示例

设：

$$
\mathbf z_d=[\sin\theta,0,\cos\theta]^T.
$$

代入上式可得：

$$
\boldsymbol\omega_W=[0,\dot\theta,0]^T.
$$

经过 $\theta=180^\circ$ 时，角速度仍然正常，不会仅因为推力向下而冻结。对纯俯仰姿态，转到机体系后该俯仰角速度同样保持一致。

[方向导数与角速度计算](../include/continuous_flatness.h)

## 6. 前馈如何进入 MPC

本次没有把 MPC 改成独立的翻转控制器。实际姿态 $\hat R$ 与参考姿态 $R_d$ 的误差仍使用：

$$
\delta\boldsymbol\theta=\log(\hat R^TR_d).
$$

优化变量的输入误差定义是“参考减实际”，因此最终角速度命令为：

$$
\boldsymbol\omega_{\mathrm{cmd}}
=\boldsymbol\omega_{\mathrm{ref}}
-\delta\boldsymbol\omega^*.
$$

新的角速度前馈负责提供轨迹本身所需的转动，MPC 负责修正位置、速度和姿态误差。参考角速度也参与预测误差动态中的旋转演化。

[实际姿态误差计算](../include/ommpc_controller.hpp) 与 [最终角速度命令](../include/ommpc_controller.hpp)。

## 7. 配套修改：导数来源与正推力大小

### 7.1 一致的轨迹导数

新文本路径接收一致的 $p,v,a,j$。新翻转轨迹提供显式加速度和 jerk，使用 C3 插值连接入圈、圆周和出圈段，避免把有噪声的重复差分直接用于角速度前馈。

新的 14 列格式为：

```text
px py pz vx vy vz yaw yawdot ax ay az jx jy jz
```

原 8 列格式仍受支持，其导数通过平滑重建生成。当前连续模式不使用任意给定的 yaw/yawdot 轨迹，见下一节的适用边界。

[文本轨迹读取与导数处理](../include/text_trajectory.h)

### 7.2 推力大小与参考姿态一致

推力参考从原来的：

```cpp
thracc = des_acc_in_world.dot(body_z);
```

改为：

```cpp
thracc = des_acc_in_world.norm();
```

原投影使用实际或上一参考的机体方向。当期望推力已向下、实际飞机还向上时，投影可能为负，与正在要求翻转的参考姿态不一致。

现在使用正推力大小与指向期望力方向的参考姿态共同描述期望加速度。向下的推力通过倒置机体实现，不需要负的电机推力。实际控制仍受输入限制和飞机执行器动态约束。

[当前推力参考与姿态构造](../include/ommpc_controller.hpp)

## 8. 当前模式的选择与适用边界

高机动配置启用：

```yaml
use_fix_yaw: true
continuous_attitude_reference: true
```

旧 Hopf 正常分支仍然保留，便于保持其他示例的航向约定；当前代码也把旧模式的近极点备用分支改为调用连续构造函数。高机动模式则在整个参考生成过程中使用连续构造。

当前连续模式采用绕推力轴零自转的约定。这不等于任意三维运动中 Euler yaw 始终固定，也没有实现任意指定的 yaw 轨迹。若要跟踪独立指定的航向，还需要加入绕推力轴的旋转，并使相应角速度前馈与之匹配。

此外，该实现不能消除零推力的物理退化。当 $\|\mathbf f_d\|$ 接近零时，推力方向本来就不唯一，方向导数也可能变得很大。当前函数对非有限输入或 $\|\mathbf f_d\|<10^{-6}$ 返回失败；轨迹规划仍应避开不具有稳定推力方向和可行角速度的区段。

[高机动参数](../config/params_flip.yaml) 与 [模式参数读取](../src/ommpc_example.cpp)。

## 9. 验证结果与结论范围

连续姿态相关的三项单元测试已通过：

1. 连续 720° 的推力方向变化中，参考机体 Z 轴持续对齐推力方向，四元数保持相邻表示，角速度与指定的 4 rad/s 一致。
2. 三维方向路径中，解析角速度与相邻姿态的 SO(3) 差分一致。
3. 完全相反的方向能够得到确定姿态，零推力参考被拒绝。

[测试代码](../test/continuous_flatness_test.cpp) 与 测试结果（本地归档：`build/test_results/ommpc_controller/gtest-continuous_flatness_test.xml`） 可复查。

这部分修改解决的是倒置附近参考构造和前馈一致性的问题。之前原始变速圆约 20 s 的失稳发生在首次需要向下推力之前，不能全部归因于 Hopf 极点处理。成功完成实际翻转还依赖轨迹可行性、动力余量、内外环匹配和反馈时序。

相关飞行验证记录：

- 连续翻转与原有轨迹回归报告（本地归档：`test_results/flip_tracking_20261013/report.md`）
- PX4 估计反馈测试（本地归档：`test_results/flip_px4_estimate_20261013/report.md`）
- GPS 延迟与反馈时序诊断（本地归档：`test_results/flip_px4_estimate_20261013/diagnostics/report.md`）
- EKF GPS 延迟 120 ms 对照测试（本地归档：`test_results/flip_px4_delay120_20261014/report.md`）
