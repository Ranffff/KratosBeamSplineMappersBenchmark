from pathlib import Path
import math
import os


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = SCRIPT_DIR.parent
PLOT_SUBDIR = "Summary_Plots"

DISPLACEMENT_CASE_DIR = OUTPUT_ROOT / "Displacement"
ROTATION_CASE_DIR = OUTPUT_ROOT / "Rotation_1x5"
TORSION_CASE_DIR = OUTPUT_ROOT / "Torsion"

DISPLACEMENT_SUMMARY = DISPLACEMENT_CASE_DIR / "sweep_summary.txt"
ROTATION_SUMMARY = ROTATION_CASE_DIR / "sweep_summary.txt"
TORSION_SUMMARY = TORSION_CASE_DIR / "sweep_summary.txt"

DISPLACEMENT_PLOT_OUTPUT = DISPLACEMENT_CASE_DIR / PLOT_SUBDIR
ROTATION_PLOT_OUTPUT = ROTATION_CASE_DIR / PLOT_SUBDIR
TORSION_PLOT_OUTPUT = TORSION_CASE_DIR / PLOT_SUBDIR

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def parse_summary(summary_path):
    records = []
    if not summary_path.exists():
        return records

    current_header = None
    for raw_line in summary_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("thresholds_") or "%" in line:
            continue

        parts = line.split()
        if parts[0].startswith("tip_"):
            current_header = parts
            continue
        if current_header is None or len(parts) < len(current_header):
            continue

        record = {}
        for name, value in zip(current_header, parts):
            record[name] = parse_value(value)
        records.append(record)

    return records


def parse_value(value):
    try:
        if value.lower() in {"none", "not_reached"}:
            return math.nan
        parsed = float(value)
        if parsed.is_integer() and "e" not in value.lower() and "." not in value:
            return int(parsed)
        return parsed
    except ValueError:
        return value


def consolidate_records(records, parameter_name):
    consolidated = {}
    for record in records:
        if parameter_name not in record:
            continue
        parameter_value = float(record[parameter_name])
        current = consolidated.get(parameter_value, {})
        current.update(record)
        consolidated[parameter_value] = current
    return [consolidated[key] for key in sorted(consolidated)]


def get_series(records, x_name, y_name):
    x_values = []
    y_values = []
    for record in records:
        if x_name not in record or y_name not in record:
            continue
        y_value = record[y_name]
        if isinstance(y_value, float) and math.isnan(y_value):
            continue
        x_values.append(float(record[x_name]))
        y_values.append(float(y_value))
    return x_values, y_values


def plot_absolute_errors(records, x_name, x_label, title, output_path):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    plotted = False

    for y_name, label, marker in [
        ("spline_error", "beam_spline_mapper vs co-rotation", "o"),
        ("linear_error", "linear beam mapper vs co-rotation", "s"),
    ]:
        x_values, y_values = get_series(records, x_name, y_name)
        if not x_values:
            continue
        ax.plot(x_values, y_values, marker=marker, linewidth=1.8, label=label)
        plotted = True

    if not plotted:
        plt.close(fig)
        return

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("selected displacement error to co-rotation")
    ax.set_yscale("log")
    configure_x_axis(ax, x_name, records)
    configure_y_axis(ax, log_scale=True)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_percent_errors(records, x_name, x_label, title, output_path):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    plotted = False

    for y_name, label, marker in [
        ("spline_error_percent", "beam_spline_mapper vs co-rotation", "o"),
        ("linear_error_percent", "linear beam mapper vs co-rotation", "s"),
    ]:
        x_values, y_values = get_series(records, x_name, y_name)
        if not x_values:
            continue
        ax.plot(x_values, y_values, marker=marker, linewidth=1.8, label=label)
        plotted = True

    if not plotted:
        plt.close(fig)
        return

    for threshold in (1.0, 10.0):
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0, alpha=0.45)
        ax.text(
            0.99,
            threshold,
            f"{threshold:g}%",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=9,
        )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("selected normalized error to co-rotation [%]")
    configure_x_axis(ax, x_name, records)
    configure_y_axis(ax, log_scale=False)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_threshold_report(rotation_records, output_path):
    lines = []
    for mapper_name, field_name in [
        ("beam_spline_mapper", "spline_error_percent"),
        ("linear beam mapper", "linear_error_percent"),
    ]:
        lines.append(mapper_name)
        for threshold in (1.0, 10.0):
            rotation = first_parameter_at_threshold(rotation_records, "tip_rotation", field_name, threshold)
            if rotation is None:
                lines.append(f"  {threshold:g}%: not reached")
            else:
                lines.append(f"  {threshold:g}%: theta = {rotation:g}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")


def first_parameter_at_threshold(records, x_name, y_name, threshold):
    for record in records:
        if x_name not in record or y_name not in record:
            continue
        if float(record[y_name]) >= threshold:
            return float(record[x_name])
    return None


def configure_x_axis(ax, x_name, records):
    if x_name in {"tip_rotation", "tip_torsion"}:
        values = [float(record[x_name]) for record in records if x_name in record]
        if not values:
            return
        max_value = max(values)
        ax.set_xlim(0.0, max_value * 1.02)
        pi_ticks = [
            0.0,
            0.25 * math.pi,
            0.5 * math.pi,
            0.75 * math.pi,
            1.0 * math.pi,
            1.25 * math.pi,
            1.5 * math.pi,
            1.75 * math.pi,
            2.0 * math.pi,
        ]
        visible_ticks = [tick for tick in pi_ticks if tick <= max_value * 1.02]
        ax.xaxis.set_major_locator(mticker.FixedLocator(visible_ticks))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_pi_tick))
        return

    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))


def configure_y_axis(ax, log_scale):
    if log_scale:
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
        ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=tuple(range(2, 10))))
        ax.yaxis.set_major_formatter(mticker.LogFormatterSciNotation(base=10.0))
        return

    _, y_max = ax.get_ylim()
    ax.set_ylim(bottom=0.0, top=y_max * 1.03)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    if abs(y_max) < 1.0e-3:
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1e"))
    else:
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))


def format_pi_tick(value, _position):
    multiple = value / math.pi
    if math.isclose(multiple, 0.0, abs_tol=1.0e-12):
        return "0"
    if math.isclose(multiple, round(multiple), abs_tol=1.0e-12):
        integer = int(round(multiple))
        return "π" if integer == 1 else f"{integer}π"
    return f"{multiple:g}π"


def main():
    displacement_records = consolidate_records(
        parse_summary(DISPLACEMENT_SUMMARY),
        "tip_displacement",
    )
    rotation_records = consolidate_records(
        parse_summary(ROTATION_SUMMARY),
        "tip_rotation",
    )
    torsion_records = consolidate_records(
        parse_summary(TORSION_SUMMARY),
        "tip_torsion",
    )

    plot_absolute_errors(
        displacement_records,
        "tip_displacement",
        "prescribed tip displacement",
        "Mapper Error vs Prescribed Tip Displacement",
        DISPLACEMENT_PLOT_OUTPUT / "displacement_absolute_error_to_corotation.png",
    )
    plot_absolute_errors(
        rotation_records,
        "tip_rotation",
        "prescribed tip rotation [rad]",
        "Mapper Absolute Error vs Prescribed Tip Rotation",
        ROTATION_PLOT_OUTPUT / "rotation_absolute_error_to_corotation.png",
    )
    plot_percent_errors(
        rotation_records,
        "tip_rotation",
        "prescribed tip rotation [rad]",
        "Mapper Relative Error vs Prescribed Tip Rotation",
        ROTATION_PLOT_OUTPUT / "rotation_relative_error_to_corotation.png",
    )
    plot_absolute_errors(
        torsion_records,
        "tip_torsion",
        "prescribed tip torsion [rad]",
        "Mapper Absolute Error vs Prescribed Tip Torsion",
        TORSION_PLOT_OUTPUT / "torsion_absolute_error_to_corotation.png",
    )
    plot_percent_errors(
        torsion_records,
        "tip_torsion",
        "prescribed tip torsion [rad]",
        "Mapper Relative Error vs Prescribed Tip Torsion",
        TORSION_PLOT_OUTPUT / "torsion_relative_error_to_corotation.png",
    )
    write_threshold_report(rotation_records, ROTATION_PLOT_OUTPUT / "rotation_threshold_report.txt")

    print("PNG plots written to:")
    for output_dir in (DISPLACEMENT_PLOT_OUTPUT, ROTATION_PLOT_OUTPUT, TORSION_PLOT_OUTPUT):
        if output_dir.exists():
            print(" ", output_dir)


if __name__ == "__main__":
    main()
