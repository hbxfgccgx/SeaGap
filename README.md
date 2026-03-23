# SeaGap

<p align="center">
　<img src="https://user-images.githubusercontent.com/39943988/213900070-0ff33731-443c-4666-b940-e2c100b433de.png" width="200">
</p>

Welcome to SeaGap!

SeaGap is a Software of the enhanced analysis for GNSS-acoustic positioning.

See [online manual](https://f-tommy.github.io/SeaGapDocs/) in detail.

## 长基线水声定位（LBL/GNSS-A）程序位置与结构

SeaGap 中进行长基线水声定位（GNSS-A/LBL）的核心程序位于 `src/` 目录，按功能可分为以下层次：

### 1) 传播时间计算层（物理模型）
- `src/traveltime.jl`  
  负责声学传播时间计算（精确与近似），是定位反演的基础。

### 2) 定位反演层（核心定位算法）
- `src/kinematic_array.jl`：动态阵列定位（逐时/逐组估计阵列位移）
- `src/kinematic_array_3d.jl`：三维动态阵列定位扩展
- `src/static_array.jl`：静态阵列定位（整期数据联合反演）
- `src/static_array_s.jl`, `src/static_array_AICBIC.jl`, `src/static_array_s_ABIC.jl`
- `src/static_array_grad.jl`, `src/static_array_TR.jl`
- `src/static_array_mcmcgrad.jl`, `src/static_array_mcmcgradc.jl`, `src/static_array_mcmcgradv.jl`  
  上述文件提供静态定位的不同约束、模型选择与 MCMC 扩展方案。

### 3) 反演支撑层（基函数与约束）
- `src/ntdbasis.jl`：时间 B 样条基函数（用于声速时变建模）
- `src/inv_func.jl`：平滑约束/反演辅助函数
- `src/ttres.jl`：传播时间残差计算

### 4) 数据与流程层（输入输出和预处理）
- `src/read_gnssa.jl`, `src/obsdata_format.jl`：观测数据读取与格式化
- `src/anttena2tr.jl`：GNSS 天线到换能器坐标换算

### 5) 结果分析与可视化层
- `src/plot_*.jl`, `src/position_kinematic.jl`, `src/convert_displacement.jl` 等  
  用于定位结果后处理、时序分析与绘图。

主入口在 `src/SeaGap.jl`：该文件通过 `include(...)` 组织上述模块，并导出定位与分析函数。
