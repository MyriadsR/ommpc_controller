# 近期代码修改说明与验证

日期：2026-10-14。以下说明本次上传的已有修改及新增文件；飞行日志、工作空间构建产物和测试结果目录不纳入仓库。

## 1. 轨迹播放、参考构造与推力约束

| 文件 | 原有问题或行为 | 当前处理 |
|---|---|---|
| `include/text_trajectory.h`（新增） | 文本轨迹依赖行计数；速度差分与舍弃 jerk 会损失前馈一致性 | 校验有限数值与 8/14 列格式；8 列局部拟合导数，14 列保留显式 a/j；七次 Hermite 保持 C3；末端悬停 |
| `src/ommpc_example.cpp` | 回调计数播放；发布饱和与估计器输入不一致；高处降落固定下降距离不足 | 按 ROS 已过时间采样预测窗口；控制周期由 step_T 决定；记录最终发布推力；降落参考到 local origin 以下；PX4 自动解除武装后清理状态 |
| `include/continuous_flatness.h`（新增） | 向下推力附近 Hopf 单图姿态参考有奇异处理 | 平行运输参考机体轴、保持四元数符号连续；以 `z × dz/dt` 构造角速度前馈；拒绝零推力方向 |
| `include/ommpc_controller.hpp` | 推力大小投影到实际机体轴；预测末端姿态可能污染当前参考；QP 与实际限幅不一致 | 采用推力向量范数；仅缓存当前时刻参考；QP 取实际归一化限幅与物理加速度限制交集；轨迹专用最小推力不影响悬停/起降；保留默认线性模型，支持可选非线性校准 |
| `include/so3_math.h` | acos 域误差及接近 π 时的 sin 分母不稳定 | 用 Eigen AngleAxis 计算 SO(3) 对数，并保留小角度展开 |
| `config/params.yaml`、`params_vertical.yaml` | 发布限幅不显式进入 QP 配置 | 显式声明原归一化推力界 0.04–0.9 |
| `test/text_trajectory_test.cpp` | 新轨迹读取与导数需要验证 | 格式/非法行、恒速、圆运动导数和末端保持检查 |
| `test/continuous_flatness_test.cpp` | 连续翻转需验证几何一致性 | 完整俯仰回转、三维路径角速度一致性、反向与零推力检查 |
| `test/thrust_model_test.cpp` | 推力映射、限幅及起降最低推力需要验证 | 线性兼容、非线性互逆与实际限幅、轨迹保留推力不影响悬停/起降 |
| `CMakeLists.txt` | 没有上述回归检查入口 | 注册 gtest 和 Python nose 测试 |

连续姿态模式默认关闭，只在实验翻转配置启用。该模式的固定 yaw 约定是沿推力轴不自旋，不能等同于任意三维路径的绝对 Euler yaw 固定。原 Hopf 路径仍可使用。

## 2. 高机动垂直圆仿真

| 文件 | 用途与关键约定 |
|---|---|
| `config/params_flip.yaml` | 实验控制权重、正推力/角速度范围与二次推力校准；在线校准只在该配置关闭 |
| `traj/generate_vertical_flip.py`、`vertical_flip.traj`、`vertical_flip.json` | 1.2 m 半径、1.68 s 周期、三圈 XZ 圆，C3 进入退出，提供一致 p/v/a/j 和可行性元数据 |
| `sim/generate_agile_model.py`、`iris_agile.sdf` | 从 stock Iris 单独生成电机推力常数 ×2.4 的实验模型；质量/惯量保留，stock PX4 文件不改动 |
| `sim/ground_truth_odometry.py` | 显式理想反馈对照；Gazebo 世界速度转成控制器需要的 body FLU 速度 |
| `launch/flip_controller.launch` | 载入实验配置，可选择反馈话题及参考文件 |
| `launch/aggressive_flip_sitl.launch` | 实验 SITL 启动；区分理想/估计反馈；EKF 启动前设置仿真 GPS 延迟 120 ms，可覆盖 |
| `test/analyze_flip.py` | 逐圈展开俯仰旋转与圆顶倒置检查、误差及 IMU/控制量图 |

`iris_agile` 是显式增强动力的实验仿真模型，并非论文实机复刻。预测步数仍为现有 N=20，并非论文的 N=8。`analyze_flip.py` 的估计反馈坐标评价只做初始平移、假定轴对齐；精细 ENU yaw/平移对齐诊断保存在本地测试归档，不把此简化指标冒称为完整坐标校正结果。

## 3. PX4 100 Hz 同步估计反馈

| 文件 | 改动 |
|---|---|
| `sim/px4_synchronized_odometry.py` | 取单份 PX4 ODOMETRY 中的时间戳、位置、姿态和速度，不做预测/插值或真值注入；从原始 MAVLink 检查 frame ID |
| `sim/px4_odometry_math.py` | 当前 PX4 输出世界 NED 速度，而 MAVROS 1.20.x odometry 插件按 body FRD 转换；撤销该映射，转成 ENU 再依据同一姿态转成 body FLU；已为 body FRD 的输入保留 |
| `test/px4_odometry_math_test.py` | 检查 yaw、90° pitch、倒飞、原生 body frame、非法输入和未知坐标约定 |
| `test/run_tracking_sitl.py` | 自动起飞/跟踪/降落及日志采集；可请求三个 MAVLink 流 100 Hz；记录实际频率所需数据、原 PX4 参数及推力模型；解除武装后恢复实验内环参数 |
| `package.xml` | 补充仿真 Python 工具的 ROS/数值/绘图依赖和 nose 测试依赖；同步节点另需 pymavlink |

消息频率请求不能提高 GNSS 原始观测频率，也不会改变 EKF 算法。改变 stream rate 后应重启 SITL 再跑 baseline。升级 MAVROS 时需重新核对速度转换约定。

## 4. 五圈竖直螺旋翻转

| 文件 | 改动 |
|---|---|
| `traj/generate_vertical_helix.py`、`vertical_helix.traj`、`vertical_helix.json` | 论文 Fig. 7(a) 同类 YZ 竖直圆，沿 X 轴移动；五圈横滚翻转，半径 1.2 m、周期 1.68 s、圆段轴向位移 -2.4 m，总长 12.4 s；连接与输入可行性检查 |
| `test/analyze_helix.py` | 起飞前固定 yaw/平移对齐；Gazebo 真值仅用于离线评价；逐圈连续横滚、圆顶向下推力检查，SO(3) 姿态误差与论文式图 |
| `launch/*flip*.launch` | `ref_filename` 可切换三圈俯仰圆与五圈横滚螺旋；默认仍为三圈圆 |

螺旋参考为自行构造的同类轨迹，不是从论文图片恢复的原始数据。进入/退出段高度与论文示意图不同。

## 5. 原有加速圆脚本与文档

- `traj/circle_vertical.py`：保留此前 100 s、0.5→13 m/s 的加速压力测试设置。
- `traj/var_circle.py`：保留此前 80 s、0.5→10 m/s、起始高度 1.5 m；输出固定到脚本目录下的 `circle_horizon.txt`。
- 上述压力测试后段可能超出可用推力，和经过可行性检查的固定速翻转参考分开说明，不能把它们当作已通过的实机轨迹。
- `doc/mpc.md`：修正控制边界约束行数：每个输入一行同时表示上下界；20 步时 80 行控制界，总计 269 行约束。
- `doc/attitude_reference_improvement_report.md`：详细解释原 Hopf 逻辑、连续姿态和角速度前馈；仓库源码链接改为相对路径，本地归档标明未上传。
- `README.md`：参数兼容性、实验假设、依赖、启动和分析入口。
- `.gitignore`：忽略 Python 缓存，保留显式 `.traj/.json` 参考数据。

## 验证与已知限制

发布前重新运行编译和 15 个实际单元测试（4 轨迹、3 连续姿态、3 推力、5 坐标转换）。catkin 汇总 XML 的计数显示 25 tests；按独立测试用例为 15。另检查 Python/XML、两套生成器与预生成参考一致、C3 和密集采样输入限制，以及 git whitespace。

本地已有实际 SITL 验证：

| 测试 | 真值位置 RMSE | 最大误差 | 翻转检查 |
|---|---:|---:|---|
| 原约 30 Hz PX4 反馈 | 0.407 m | 0.999 m | 2/3 圈 |
| 100 Hz 原 MAVROS 组合反馈 | 0.407 m | 0.747 m | 3/3 圈 |
| 100 Hz 单消息同步反馈，两次 | 0.346 / 0.306 m | 0.720 / 0.752 m | 两次 3/3 圈 |
| 五圈螺旋，两次 | 0.241 / 0.318 m | 0.601 / 0.615 m | 两次 5/5 圈 |

此处不同轨迹不是同一受控精度比较。飞行重复数有限；首圈相位误差和位置估计/执行链滞后仍存在，不能据此声称达到论文平均 0.08 m 的实机精度。飞行原始日志已按用户要求清理；本地报告/结果图/汇总保留，但不上传到本仓库。
