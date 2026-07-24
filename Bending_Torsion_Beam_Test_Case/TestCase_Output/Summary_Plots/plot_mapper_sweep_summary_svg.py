from pathlib import Path
import html
import math


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = SCRIPT_DIR.parent
PLOT_SUBDIR = "Summary_Plots"

DISPLACEMENT_CASE_DIR = OUTPUT_ROOT / "Displacement"
ROTATION_CASE_DIR = OUTPUT_ROOT / "Rotation_10x100"
TORSION_CASE_DIR = OUTPUT_ROOT / "Torsion"

DISPLACEMENT_SUMMARY = DISPLACEMENT_CASE_DIR / "sweep_summary.txt"
ROTATION_SUMMARY = ROTATION_CASE_DIR / "sweep_summary.txt"
TORSION_SUMMARY = TORSION_CASE_DIR / "sweep_summary.txt"

DISPLACEMENT_PLOT_OUTPUT = DISPLACEMENT_CASE_DIR / PLOT_SUBDIR
ROTATION_PLOT_OUTPUT = ROTATION_CASE_DIR / PLOT_SUBDIR
TORSION_PLOT_OUTPUT = TORSION_CASE_DIR / PLOT_SUBDIR


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
    series = collect_series(records, x_name, [
        ("spline_error", "beam_spline_mapper vs co-rotation", "o"),
        ("linear_error", "linear beam mapper vs co-rotation", "s"),
    ])
    if not series:
        return

    write_svg_plot(
        output_path,
        title,
        x_label,
        "selected displacement error to co-rotation",
        series,
        log_y=True,
        angle_axis=x_name in {"tip_rotation", "tip_torsion"},
    )


def plot_percent_errors(records, x_name, x_label, title, output_path):
    series = collect_series(records, x_name, [
        ("spline_error_percent", "beam_spline_mapper vs co-rotation", "o"),
        ("linear_error_percent", "linear beam mapper vs co-rotation", "s"),
    ])
    if not series:
        return

    write_svg_plot(
        output_path,
        title,
        x_label,
        "selected normalized error to co-rotation [%]",
        series,
        thresholds=[1.0, 10.0],
        log_y=False,
        angle_axis=x_name in {"tip_rotation", "tip_torsion"},
    )


def collect_series(records, x_name, specs):
    collected = []
    used_labels = set()
    for y_name, label, marker in specs:
        if label in used_labels:
            continue
        x_values, y_values = get_series(records, x_name, y_name)
        if not x_values:
            continue
        collected.append({
            "label": label,
            "marker": marker,
            "x": x_values,
            "y": y_values,
        })
        used_labels.add(label)
    return collected


def write_svg_plot(output_path, title, x_label, y_label, series, thresholds=None, log_y=False, angle_axis=False):
    width = 960
    height = 600
    margin_left = 90
    margin_right = 230
    margin_top = 70
    margin_bottom = 80
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom

    all_x = [x for item in series for x in item["x"]]
    all_y = [y for item in series for y in item["y"] if y > 0.0]
    if thresholds:
        all_y.extend(value for value in thresholds if value > 0.0)
    if not all_x or not all_y:
        return

    x_min, x_max = min(all_x), max(all_x)
    if angle_axis:
        x_min = 0.0
        x_max *= 1.02
    if math.isclose(x_min, x_max):
        x_min -= 1.0
        x_max += 1.0

    if log_y:
        y_min = 10.0 ** math.floor(math.log10(min(all_y)))
        y_max = 10.0 ** math.ceil(math.log10(max(all_y)))
    else:
        y_min = 0.0
        y_max = max(all_y) * 1.08
        if math.isclose(y_min, y_max):
            y_max = 1.0

    def x_to_px(x_value):
        return margin_left + (x_value - x_min) / (x_max - x_min) * plot_width

    def y_to_px(y_value):
        if log_y:
            safe_value = max(y_value, y_min)
            t = (math.log10(safe_value) - math.log10(y_min)) / (math.log10(y_max) - math.log10(y_min))
        else:
            t = (y_value - y_min) / (y_max - y_min)
        return margin_top + (1.0 - t) * plot_height

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="{width / 2}" y="32" text-anchor="middle" font-family="Arial" font-size="20">{escape(title)}</text>')

    x0, y0 = margin_left, margin_top + plot_height
    svg.append(f'<line x1="{x0}" y1="{margin_top}" x2="{x0}" y2="{y0}" stroke="#222" stroke-width="1.2"/>')
    svg.append(f'<line x1="{x0}" y1="{y0}" x2="{margin_left + plot_width}" y2="{y0}" stroke="#222" stroke-width="1.2"/>')

    for tick in make_x_ticks(x_min, x_max, angle_axis):
        px = x_to_px(tick)
        svg.append(f'<line x1="{px:.2f}" y1="{y0}" x2="{px:.2f}" y2="{y0 + 5}" stroke="#222"/>')
        svg.append(f'<line x1="{px:.2f}" y1="{margin_top}" x2="{px:.2f}" y2="{y0}" stroke="#ddd"/>')
        svg.append(f'<text x="{px:.2f}" y="{y0 + 24}" text-anchor="middle" font-family="Arial" font-size="12">{format_x_tick(tick, angle_axis)}</text>')

    y_ticks = make_log_ticks(y_min, y_max) if log_y else make_nice_ticks(y_min, y_max, 6)
    for tick in y_ticks:
        py = y_to_px(tick)
        svg.append(f'<line x1="{x0 - 5}" y1="{py:.2f}" x2="{x0}" y2="{py:.2f}" stroke="#222"/>')
        svg.append(f'<line x1="{x0}" y1="{py:.2f}" x2="{margin_left + plot_width}" y2="{py:.2f}" stroke="#ddd"/>')
        svg.append(f'<text x="{x0 - 10}" y="{py + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12">{format_number(tick)}</text>')

    if thresholds:
        for threshold in thresholds:
            if threshold < y_min or threshold > y_max:
                continue
            py = y_to_px(threshold)
            svg.append(f'<line x1="{x0}" y1="{py:.2f}" x2="{margin_left + plot_width}" y2="{py:.2f}" stroke="#333" stroke-dasharray="6 5" opacity="0.65"/>')
            svg.append(f'<text x="{margin_left + plot_width - 8}" y="{py - 5:.2f}" text-anchor="end" font-family="Arial" font-size="12">{threshold:g}%</text>')

    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
    for index, item in enumerate(series):
        color = colors[index % len(colors)]
        points = " ".join(f'{x_to_px(x):.2f},{y_to_px(y):.2f}' for x, y in zip(item["x"], item["y"]))
        svg.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.2" points="{points}"/>')
        for x, y in zip(item["x"], item["y"]):
            px, py = x_to_px(x), y_to_px(y)
            if item["marker"] == "s":
                svg.append(f'<rect x="{px - 3:.2f}" y="{py - 3:.2f}" width="6" height="6" fill="{color}"/>')
            else:
                svg.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="3.4" fill="{color}"/>')

        legend_x = margin_left + plot_width + 26
        legend_y = margin_top + 24 + index * 28
        svg.append(f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 24}" y2="{legend_y}" stroke="{color}" stroke-width="2.2"/>')
        svg.append(f'<text x="{legend_x + 32}" y="{legend_y + 4}" font-family="Arial" font-size="12">{escape(item["label"])}</text>')

    svg.append(f'<text x="{margin_left + plot_width / 2}" y="{height - 25}" text-anchor="middle" font-family="Arial" font-size="14">{escape(x_label)}</text>')
    svg.append(
        f'<text x="22" y="{margin_top + plot_height / 2}" text-anchor="middle" '
        f'transform="rotate(-90 22 {margin_top + plot_height / 2})" font-family="Arial" font-size="14">{escape(y_label)}</text>'
    )
    svg.append("</svg>")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(svg) + "\n")


def make_x_ticks(minimum, maximum, angle_axis):
    if not angle_axis:
        return make_nice_ticks(minimum, maximum, 7)

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
    return [tick for tick in pi_ticks if minimum <= tick <= maximum]


def make_ticks(minimum, maximum, count):
    if count <= 1:
        return [minimum]
    step = (maximum - minimum) / (count - 1)
    return [minimum + i * step for i in range(count)]


def make_nice_ticks(minimum, maximum, count):
    if count <= 1 or math.isclose(minimum, maximum):
        return [minimum]

    raw_step = (maximum - minimum) / (count - 1)
    magnitude = 10.0 ** math.floor(math.log10(abs(raw_step)))
    normalized = raw_step / magnitude
    if normalized <= 1.0:
        nice_normalized = 1.0
    elif normalized <= 2.0:
        nice_normalized = 2.0
    elif normalized <= 2.5:
        nice_normalized = 2.5
    elif normalized <= 5.0:
        nice_normalized = 5.0
    else:
        nice_normalized = 10.0

    step = nice_normalized * magnitude
    start = math.floor(minimum / step) * step
    end = math.ceil(maximum / step) * step

    ticks = []
    value = start
    while value <= end + 0.5 * step:
        if value >= minimum - 0.5 * step and value <= maximum + 0.5 * step:
            ticks.append(0.0 if math.isclose(value, 0.0, abs_tol=1.0e-15) else value)
        value += step
    return ticks


def make_log_ticks(minimum, maximum):
    min_power = math.floor(math.log10(minimum))
    max_power = math.ceil(math.log10(maximum))
    return [10.0 ** power for power in range(min_power, max_power + 1)]


def format_number(value):
    if value == 0.0:
        return "0"
    if abs(value) >= 1.0e4 or abs(value) < 1.0e-3:
        return f"{value:.1e}"
    return f"{value:g}"


def format_x_tick(value, angle_axis):
    if not angle_axis:
        return format_number(value)

    multiple = value / math.pi
    if math.isclose(multiple, 0.0, abs_tol=1.0e-12):
        return "0"
    if math.isclose(multiple, round(multiple), abs_tol=1.0e-12):
        integer = int(round(multiple))
        return "π" if integer == 1 else f"{integer}π"
    return f"{multiple:g}π"


def escape(value):
    return html.escape(str(value))


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
        DISPLACEMENT_PLOT_OUTPUT / "displacement_absolute_error_to_corotation.svg",
    )
    plot_absolute_errors(
        rotation_records,
        "tip_rotation",
        "prescribed tip rotation [rad]",
        "Mapper Absolute Error vs Prescribed Tip Rotation",
        ROTATION_PLOT_OUTPUT / "rotation_absolute_error_to_corotation.svg",
    )
    plot_percent_errors(
        rotation_records,
        "tip_rotation",
        "prescribed tip rotation [rad]",
        "Mapper Relative Error vs Prescribed Tip Rotation",
        ROTATION_PLOT_OUTPUT / "rotation_relative_error_to_corotation.svg",
    )
    plot_absolute_errors(
        torsion_records,
        "tip_torsion",
        "prescribed tip torsion [rad]",
        "Mapper Absolute Error vs Prescribed Tip Torsion",
        TORSION_PLOT_OUTPUT / "torsion_absolute_error_to_corotation.svg",
    )
    plot_percent_errors(
        torsion_records,
        "tip_torsion",
        "prescribed tip torsion [rad]",
        "Mapper Relative Error vs Prescribed Tip Torsion",
        TORSION_PLOT_OUTPUT / "torsion_relative_error_to_corotation.svg",
    )

    print("SVG plots written to:")
    for output_dir in (DISPLACEMENT_PLOT_OUTPUT, ROTATION_PLOT_OUTPUT, TORSION_PLOT_OUTPUT):
        if output_dir.exists():
            print(" ", output_dir)


if __name__ == "__main__":
    main()
