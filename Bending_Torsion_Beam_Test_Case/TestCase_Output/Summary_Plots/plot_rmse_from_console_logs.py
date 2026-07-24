from argparse import ArgumentParser
from pathlib import Path
import csv
import math
import os
import re


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_DIR = SCRIPT_DIR.parent / "Rotation_1x5"
OUTPUT_SUBDIR = "summary_plots"

os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PARAMETER_NAMES = (
    "tip_rotation",
    "tip_torsion",
    "tip_displacement",
    "tip_rotation_v",
    "tip_displacement_v",
    "tip_rotation_w",
    "tip_displacement_w",
)


def parse_arguments():
    parser = ArgumentParser(
        description="RMSE"
    )
    parser.add_argument(
        "result_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_RESULT_DIR,
        help="Simulation result directory containing case folders with console_log.txt.",
    )
    parser.add_argument(
        "--output-subdir",
        default=OUTPUT_SUBDIR,
        help="Output folder created below the simulation result directory.",
    )
    return parser.parse_args()


def parse_console_log(log_path):
    text = log_path.read_text(errors="replace")
    error_type = extract_text_value(text, "error_type")
    selected_error = extract_text_value(text, "selected_error")
    if error_type and error_type != "normalized_rmse":
        return None
    if selected_error and selected_error.lower() != "normalized rmse":
        return None

    record = {
        "case_name": extract_text_value(text, "case_name") or log_path.parent.name,
        "log_path": str(log_path),
    }
    for name in PARAMETER_NAMES:
        value = extract_float_value(text, name)
        if value is not None:
            record[name] = value

    reference_rmse = extract_float_value(text, "reference_scale")
    if reference_rmse is not None:
        record["reference_rmse"] = reference_rmse

    table_values = parse_error_summary_table(text)
    if table_values:
        record.update(table_values)
    else:
        record.update(parse_legacy_error_summary(text))

    required = {
        "spline_rmse",
        "spline_nrmse_percent",
        "linear_rmse",
        "linear_nrmse_percent",
    }
    if not required.issubset(record):
        return None
    return record


def extract_text_value(text, name):
    match = re.search(
        rf"^\s*{re.escape(name)}\s*:\s*(.+?)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return match.group(1).strip() if match else None


def extract_float_value(text, name):
    match = re.search(
        rf"^\s*{re.escape(name)}\s*:\s*({FLOAT_PATTERN})\s*$",
        text,
        flags=re.MULTILINE,
    )
    return float(match.group(1)) if match else None


def parse_error_summary_table(text):
    mapper_pattern = re.compile(
        rf"^\s*(beam_spline_mapper|beam_mapper_linear)\s+"
        rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})%\s*$",
        flags=re.MULTILINE,
    )
    values = {}
    for mapper_name, absolute_error, relative_error in mapper_pattern.findall(text):
        prefix = "spline" if mapper_name == "beam_spline_mapper" else "linear"
        values[f"{prefix}_rmse"] = float(absolute_error)
        values[f"{prefix}_nrmse_percent"] = float(relative_error)
    return values


def parse_legacy_error_summary(text):
    values = {}
    for prefix, label in (
        ("spline", "spline_error_to_corotation"),
        ("linear", "linear_error_to_corotation"),
    ):
        match = re.search(
            rf"^\s*{label}\s*:\s*({FLOAT_PATTERN})\s+"
            rf"percent\s*:\s*({FLOAT_PATTERN})\s*$",
            text,
            flags=re.MULTILINE,
        )
        if match:
            values[f"{prefix}_rmse"] = float(match.group(1))
            values[f"{prefix}_nrmse_percent"] = float(match.group(2))
    return values


def collect_records(result_dir):
    records = []
    for log_path in sorted(result_dir.rglob("console_log.txt")):
        if OUTPUT_SUBDIR in log_path.parts:
            continue
        record = parse_console_log(log_path)
        if record:
            records.append(record)
    return records


def select_x_parameter(records):
    for name in ("tip_rotation", "tip_torsion", "tip_displacement"):
        if all(name in record for record in records) and len(
            {record[name] for record in records}
        ) > 1:
            return name
    return None


def write_csv(records, output_path):
    fieldnames = [
        "case_name",
        *PARAMETER_NAMES,
        "reference_rmse",
        "spline_rmse",
        "spline_nrmse_percent",
        "linear_rmse",
        "linear_nrmse_percent",
        "log_path",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def plot_sweep(records, x_name, result_name, output_dir):
    records = sorted(records, key=lambda record: record[x_name])
    x_values = [record[x_name] for record in records]

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.plot(
        x_values,
        [record["spline_rmse"] for record in records],
        marker="o",
        linewidth=1.8,
        label="beam_spline_mapper vs co-rotation",
    )
    ax.plot(
        x_values,
        [record["linear_rmse"] for record in records],
        marker="s",
        linewidth=1.8,
        label="linear beam mapper vs co-rotation",
    )
    ax.set_title("RMSE")
    ax.set_xlabel(parameter_label(x_name))
    ax.set_ylabel("displacement RMSE")
    ax.set_yscale("log")
    configure_x_axis(ax, x_name, x_values)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "rmse_from_console_logs.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.plot(
        x_values,
        [record["spline_nrmse_percent"] for record in records],
        marker="o",
        linewidth=1.8,
        label="beam_spline_mapper vs co-rotation",
    )
    ax.plot(
        x_values,
        [record["linear_nrmse_percent"] for record in records],
        marker="s",
        linewidth=1.8,
        label="linear beam mapper vs co-rotation",
    )
    for threshold in (1.0, 10.0):
        ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0, alpha=0.4)
    ax.set_title("Normalized RMSE")
    ax.set_xlabel(parameter_label(x_name))
    ax.set_ylabel("normalized RMSE [%]")
    ax.set_ylim(bottom=0.0)
    configure_x_axis(ax, x_name, x_values)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "normalized_rmse_from_console_logs.png", dpi=220)
    plt.close(fig)


def plot_categorical(records, result_name, output_dir):
    labels = [record["case_name"].replace("_", " ") for record in records]
    positions = list(range(len(records)))
    width = 0.36

    fig, ax = plt.subplots(figsize=(max(9.0, 1.8 * len(records)), 5.8))
    ax.bar(
        [position - width / 2 for position in positions],
        [record["spline_rmse"] for record in records],
        width,
        label="beam_spline_mapper vs co-rotation",
    )
    ax.bar(
        [position + width / 2 for position in positions],
        [record["linear_rmse"] for record in records],
        width,
        label="linear beam mapper vs co-rotation",
    )
    ax.set_title("RMSE")
    ax.set_ylabel("displacement RMSE")
    ax.set_xticks(positions, labels, rotation=20, ha="right")
    ax.set_ylim(bottom=0.0)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "rmse_from_console_logs.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(9.0, 1.8 * len(records)), 5.8))
    ax.bar(
        [position - width / 2 for position in positions],
        [record["spline_nrmse_percent"] for record in records],
        width,
        label="beam_spline_mapper vs co-rotation",
    )
    ax.bar(
        [position + width / 2 for position in positions],
        [record["linear_nrmse_percent"] for record in records],
        width,
        label="linear beam mapper vs co-rotation",
    )
    ax.set_title("Normalized RMSE")
    ax.set_ylabel("normalized RMSE [%]")
    ax.set_xticks(positions, labels, rotation=20, ha="right")
    ax.set_ylim(bottom=0.0)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "normalized_rmse_from_console_logs.png", dpi=220)
    plt.close(fig)


def configure_x_axis(ax, x_name, x_values):
    if x_name in {"tip_rotation", "tip_torsion"}:
        ax.set_xlim(0.0, max(x_values) * 1.02)
        ticks = [
            multiple * math.pi
            for multiple in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
            if multiple * math.pi <= max(x_values) * 1.02
        ]
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_pi_tick))
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=7))


def parameter_label(name):
    return {
        "tip_rotation": "prescribed tip rotation [rad]",
        "tip_torsion": "prescribed tip torsion [rad]",
        "tip_displacement": "prescribed tip displacement",
    }[name]


def format_pi_tick(value, _position):
    multiple = value / math.pi
    if math.isclose(multiple, 0.0, abs_tol=1.0e-12):
        return "0"
    if math.isclose(multiple, 1.0, abs_tol=1.0e-12):
        return "π"
    if math.isclose(multiple, round(multiple), abs_tol=1.0e-12):
        return f"{int(round(multiple))}π"
    return f"{multiple:g}π"


def main():
    arguments = parse_arguments()
    result_dir = arguments.result_dir.expanduser().resolve()
    if not result_dir.is_dir():
        raise RuntimeError(f"Result directory does not exist: {result_dir}")

    records = collect_records(result_dir)
    if not records:
        raise RuntimeError(
            f"No normalized RMSE results were found in console logs below: {result_dir}"
        )

    output_dir = result_dir / arguments.output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    x_name = select_x_parameter(records)
    if x_name:
        records.sort(key=lambda record: record[x_name])
        plot_sweep(records, x_name, result_dir.name, output_dir)
    else:
        plot_categorical(records, result_dir.name, output_dir)
    write_csv(records, output_dir / "rmse_from_console_logs.csv")

    print(f"Parsed {len(records)} console logs")
    print(f"Plots and extracted data written to: {output_dir}")


if __name__ == "__main__":
    main()
