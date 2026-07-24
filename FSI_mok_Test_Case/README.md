# FSI Mok beam-spline mapper validation

该目录只维护当前两个 beam spline mapper 的理论、analytical verification 和
Mok FSI stability evidence。官方物理参考案例仍是：

`/root/dev/Kratos/Examples/co_simulation/validation/fsi_mok`

## 核心目录

- `CoSimulation_Cases/`：NearestNeighbor、BeamMapper_CoRotation、
  BeamSplineMapper 和 rotational-recovery mapper 的可运行 case。
- `Scripts/`：唯一维护的运行、验证和诊断脚本。
- `Notes/Beam_Abstraction_Short_Note.tex`：理论、实现、benchmark 设计和结果分析。
- `Notes/Beam_Abstraction_Short_Note.pdf`：上述笔记的编译版本。
- `TestCase_Output/MapperVerification/`：当前 mapper 二进制生成的 analytical
  forward/adjoint 结果。
- `TestCase_Output/StabilityTrials_s50/`：当前 post-fix FSI gate 结果。

## Analytical verification

所有 accuracy reference 都是 prescribed analytical field，不以其他 mapper
作为参考。脚本拒绝覆盖已有输出。

```bash
cd /root/dev/Kratos/FSI_mok_Test_Case
export PYTHONPATH=/root/dev/Kratos/bin/Release

python3 Scripts/validate_beam_spline_forward.py
python3 Scripts/validate_beam_spline_adjoint.py
python3 Scripts/validate_recovery_forward.py
python3 Scripts/validate_recovery_adjoint.py
```

Recovery 验证脚本可用 `--mode small`、`--mode finite` 或默认的 `--mode both`。
Adjoint 测试比较 analytical tangent transpose 与 centered finite-difference
directional work；finite difference 只用于验证。

## FSI benchmark

通用 runner 会对输入做快照，保持 structure/load 同步 scaling，监测
NaN、coupling failure 和 repeated structural nonconvergence，并在退出时校验输入
恢复。

示例：

```bash
python3 Scripts/run_fsi_benchmark.py \
  BeamSplineMapper_WithRotationalRecovery \
  --tag recovery_small_s50_dt005_t045 \
  --end-time 4.5 \
  --dt 0.05 \
  --alpha 0.03 \
  --iterations 60 \
  --scale 50 \
  --kernel-radius 0.50 \
  --regularization 1e-8 \
  --polynomial-level 4 \
  --rotation-recovery-mode small
```

FSI 只评价 stability 和 failure mechanism；rel L2 accuracy 由 analytical
verification 给出。推荐 gate 顺序是 `0.5 → 4.5 → 10 → 25 s`，未通过前一级时
不进入更长运行。

## 结果诊断

- `analyze_structural_newton.py`：按 FSI step/coupling iteration 汇总 Newton 历史。
- `analyze_vtk_history.py`：位移、转角、载荷和 step jump。
- `check_load_consistency.py`：fluid reaction 与 beam force/moment balance。
- `analyze_recovery_conditioning.py`：recovery saddle/polynomial conditioning。

当前结论、保留结果的精确路径以及 small/finite 理论区别见主笔记。
