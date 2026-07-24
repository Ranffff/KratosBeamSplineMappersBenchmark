# BeamSplineMapper With Rotational Recovery Short Case

这是第一个 Mok FSI beam-mapper case，目标是先跑通真实 co-simulation 数据链路：

- Fluid surface: `FluidModelPart.interface`
- Beam structure: `Structure.Parts_Beam`
- Mapper DTO: `kratos_beam_mapping`
- Mapper: `beam_spline_mapper_with_recovery_of_rotations`
- Structure element: `CrBeamElement3D2N`
- Beam model: Euler/co-rotational beam approximation, not Timoshenko

当前 `end_time = 2.0`，这是 short comparison 设置，不是正式 benchmark。正式对比前需要把
nearest-neighbor baseline 完整跑到 `t = 25.0 s` 并固化保存，后续脚本只读取这个 baseline
结果，不重复运行。

## 重要记录

普通 `beam_mapper` linear case 暂时取消，只保留 co-rotation 版本。原始 solid mesh 和 beam
版结构单元类型不同，这是预期差异：

- Baseline solid: `TotalLagrangianElement2D4N`
- Beam case: `CrBeamElement3D2N`

这里暂时没有使用 Timoshenko beam。原因是当前要优先保持 beam mapper / beam spline mapper
比较链路一致；普通 beam spline mapper 与 Timoshenko 结构模型的兼容性仍需后续确认。
