# MATLAB 12-bit 梯度自适应权重

输入、输出均以 `uint16` 容器承载 u12 数值。模块从原始 Y12 生成 Sobel L1 梯度，经线性映射和 9×9 均值滤波得到权重，数学权重为 `W12/4095`。

## 运行

使用 MATLAB 和 Image Processing Toolbox。在 `Fixed_GW` 目录执行以下命令，Demo 会先运行测试，再生成整图比较及结果文件。

```matlab
report = demo_gradientWeight_fx;
```

单独运行自动化测试可使用以下命令。

```matlab
results = runtests('tests/test_gradientWeight_fx.m');
assertSuccess(results);
```

Demo 顶部的 `saveResults=true` 控制保存。图表保存到 `result/gradientWeightFixed`。原始整数矩阵、浮点参考及可复现参数保存在 `WLI_1_debug.mat`，PNG 用于显示。

当前受限执行环境将进程级 `MATLAB_PREFDIR`、`TEMP` 和 `TMP` 指向本目录的 `.matlab_runtime` 子目录。常规 MATLAB 会话可直接运行上述命令。

## 接口与文件

```matlab
fixedDir = pwd;
addpath(fullfile(fixedDir,'functions'));
addpath(fullfile(fileparts(fixedDir),'functions'));
params = prepareGradientWeightFxParams(GLow, GHigh, Wmin);
[W12, debugFx] = computeGradientWeight_fx(Y12, params);
[Wfloat, debugFloat] = computeGradientWeight_floatRef(Y12, GLow, GHigh, Wmin);
```

- `prepareGradientWeightFxParams` 接收三个有限实数标量，要求 `0 <= GLow < GHigh <= 1`、`0 <= Wmin <= 1`，量化后仍须高阈值大于低阈值。准备成功后返回一组完整配置。
- `computeGradientWeight_fx` 接收至少 3×3 的二维 `uint16` Y12，每个值为 0…4095，使用参数准备函数生成的配置。返回同尺寸 `uint16` 权重。
- `computeGradientWeight_floatRef` 复用上级项目的浮点梯度函数，明确设置原图梯度、线性映射和 `inputMax=4095`。
- `mapGradientWeight_fx` 和 `divide81Round_fx` 分别承载整数映射与常数除法，供主函数调用及直接穷举验证。
- `tests/test_gradientWeight_fx.m` 提供独立整数黄金模型与自动化测试。

第二输出 `debugFx` 包含 `Gx`、`Gy`、`G`、`G_Q4`、`W0_12`、`S9`、`S81`、`W12`。`debugFloat` 包含 `Gx`、`Gy`、`G`、`W0`、`W_smooth`。仅请求第一输出时省去调试结构体组装。

## 数值约定

`uN` 表示无符号总 N 位，`sN` 表示包含符号位的有符号总 N 位。下表分别列出逻辑位宽与 MATLAB 容器。

| 数据 | 数值范围或安全界 | 逻辑位宽 | MATLAB 容器 |
|---|---|---|---|
| Y12、Wmin12、W0、W12 | 0…4095 | u12 | uint16 |
| Sobel 差分 | −4095…4095 | s13 | int32 |
| Sobel 二倍差分 | −8190…8190 | s14 | int32 |
| Sobel 两项部分和 | −12285…12285 | s15 | int32 |
| Gx、Gy | −16380…16380 | s15 | int32 |
| abs(Gx)、abs(Gy) | 0…16380 | u14 | uint32 |
| G | 0…24570 | u15 | uint16 |
| G_Q4 | 0…393120 | u19 | uint32 |
| Tlow_Q4、Thigh_Q4 | 0…524160 | u19 | uint32 |
| D_Q4 | 1…524160 | u19 | uint32 |
| 4095−Wmin12 | 0…4095 | u12 | uint64 |
| 斜率预计算分子 | 0…4293918720 | u32 | uint64 |
| floor(D_Q4/2) | 0…262080 | u18 | uint64 |
| 斜率分子加舍入偏置 | 0…4294180800 | u32 | uint64 |
| K_Q20 | 0…4293918720 | u32 | uint32 |
| 中间段 delta_Q4 | 1…D_Q4−1 | u19 | uint64 |
| 插值乘积及加舍入偏置分子 | ≤4294705088 | u33 安全预留 | uint64 |
| 量化增量 q | 0…4095 | u12 | uint64 |
| Wmin12+q | 0…8190 保守界 | u13 | uint64 |
| 单方向9点和 S9 | 0…36855 | u16 | uint32 |
| 81点和 S81 | 0…331695 | u19 | uint32 |
| S81+40 | 40…331735 | u19 | uint64 |
| 常数828505 | 828505 | u20 | uint64 |
| (S81+40)×828505 | 33140200…274844106175 | u38 | uint64 |

Sobel 采用与原项目一致的卷积符号，归一化梯度满足 `G_float=G/32760`。3×3 Sobel 和 9×9 均值滤波均采用 `replicate` 边界。设计文档中混用了两种均值边界，本次按用户确认统一为 `replicate`。最小尺寸沿用原梯度函数的 3×3 约束。

阈值以 `round(524160*GLow)` 和 `round(524160*GHigh)` 量化，权重下限以 `round(4095*Wmin)` 量化。正数的最近整数舍入在半值处向上。

```text
D_Q4 = Thigh_Q4 - Tlow_Q4
K_Q20 = round((4095-Wmin12)*2^20/D_Q4)
```

端点分别输出 Wmin12 和4095。中间段先计算非负差值，再按以下公式运算，并限幅至 `[Wmin12,4095]`。

```text
q = (delta_Q4*K_Q20 + 524288) >> 20
W0 = clamp(Wmin12+q, Wmin12, 4095)
W12 = clamp(((S81+40)*828505) >> 26, Wmin12, 4095)
```

所有乘法和舍入偏置先在足够宽的容器中计算，缩窄前显式限幅或证明范围安全。整数数据通路使用标准整数类型、加减、比较、乘法和位移。参数预处理中的除法按一次配置更新执行。

默认 GLow、GHigh、Wmin 分别为 `0.0024894783978878089`、`0.0071887207526713029`、`0.3`，寄存器值依次为 `1305`、`3768`、`1229`，派生 `D_Q4=2463`、`K_Q20=1220146`。

## 验证及输出

自动化测试将正式实现与独立 `int64` 标量模型比较，黄金模型通过边界索引、直接81点累加和精确 `idivide` 求得期望值。测试包含正常图案、小尺寸与窄图、参数极值、非法输入，以及梯度和累加值全域穷举。

Demo 从现有 `testImg/WLI_1.jpg` 经原项目 RGB→YCbCr 函数提取 Y，再以 `uint16(round(Y8_float*4095/255))` 构造 Y12。浮点与定点参考使用同一 Y12。这组输入属于8位 JPEG 构造的12位测试数据。

指标相对于浮点参考的最近整数值计算，覆盖全图、四像素边界带、内部和距任一量化阈值一个原始梯度单位以内的区域。CSV 中 G 的 LSB 表示原始整数梯度单位，W0、Wfinal 的 LSB 表示1/4095权重。归一化误差分别除以32760和4095。空区域的统计值记为 NaN，像素数为0。

| 结果文件 | 内容 |
|---|---|
| WLI_1_Y12.png | Y12满量程归一化显示 |
| WLI_1_G_compare.png | 浮点、整数梯度及绝对误差 |
| WLI_1_W0_compare.png | 初始权重及绝对误差 |
| WLI_1_Wfinal_compare.png | 最终权重及绝对误差 |
| WLI_1_mapping_compare.png | 默认映射曲线和逐梯度误差 |
| WLI_1_metrics.csv | 各阶段、各区域误差及像素比例 |
| WLI_1_params.csv | UI参数、寄存器值、格式与单位 |
| WLI_1_tests.csv | 自动化测试逐项状态和耗时 |
| WLI_1_debug.mat | Y12、整数阶段矩阵、浮点参考、参数、指标及测试摘要 |

当图像出现超标像素时，Demo 保存 `WLI_1_outliers.csv`，其中包含像素坐标和权重误差分解。整数一致性验收要求逐点相等，默认初始映射精度目标为1 LSB，平滑后目标为2 LSB。

### 实际运行记录

2026-11-12 使用 MATLAB R2025b 和 Image Processing Toolbox 完成验证。17项自动化测试全部通过，独立整数黄金模型逐阶段逐点一致。512个二值邻域、24571个梯度映射输入及331696个除81输入均已枚举验证。

WLI_1 测试输入为1080×1428，共1542240像素，Y12范围0…4095。全图指标如下。

| 阶段 | MAE / LSB | RMSE / LSB | 最大绝对误差 / LSB | 完全相等比例 |
|---|---:|---:|---:|---:|
| 原始整数梯度 G | 0 | 0 | 0 | 100% |
| 初始权重 W0 | 0.136220692 | 0.369080875 | 1 | 86.3779308% |
| 最终权重 Wfinal | 0.236087120 | 0.485887971 | 1 | 76.3912880% |

默认映射全域最大误差为1 LSB，31个梯度输入相对于浮点量化值高1 LSB，其余相等。整图初始及最终权重超过1 LSB的像素数均为0，整数验收与默认浮点精度目标均通过。

整图 G 范围为0…12032，W0、Wfinal 范围均为1229…4095。7个 MATLAB 文件的 Code Analyzer 消息数为0。已重载生成的 MAT，核对矩阵类型、尺寸和全图精确除81结果，并验证参数 CSV 的17位有效数字回读与 MAT 中的 UI 参数完全一致。

已导出并检查梯度、初始权重、最终权重和映射曲线四张比较图。受限沙箱读取到空屏幕尺寸，创建图窗时发生图形库崩溃，最终验证使用可访问屏幕的沙箱外 MATLAB 进程完成，运行缓存和所有输出仍位于本目录。

后续 FPGA 验证需核对实际流水线边界、寄存器帧间原子更新、像素坐标及部署位宽。
