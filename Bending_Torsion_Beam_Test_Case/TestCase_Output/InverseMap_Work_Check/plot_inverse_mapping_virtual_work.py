from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


WORK_CHECK_DIR = Path(__file__).resolve().parent
TESTCASE_OUTPUT_DIR = WORK_CHECK_DIR.parent
CASE_ROOT = TESTCASE_OUTPUT_DIR.parent
REPOSITORY_ROOT = CASE_ROOT.parent.parent

WORK_SUMMARY = WORK_CHECK_DIR / "summary.txt"
COMPARISON_SUMMARY = (
    TESTCASE_OUTPUT_DIR / "InverseMap_Mapper_Comparison" / "summary.txt"
)
OUT = (
    REPOSITORY_ROOT
    / "BeamSplineDocumentation"
    / "Thesis"
    / "MasterThesis"
    / "images"
    / "verification"
    / "inverse_mapping"
)

COLORS = {
    "beam_spline_mapper": "#3B5BA7",
    "beam_mapper_linear": "#D9772A",
    "beam_mapper_corotation": "#2A9D8F",
}
LABELS = {
    "beam_spline_mapper": "BeamSplineMapper",
    "beam_mapper_linear": "BeamMapper (linear)",
    "beam_mapper_corotation": "BeamMapper (corotation)",
}


def read_table(path):
    metadata = {}
    header = None
    rows = []
    with path.open() as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                content = line[1:].strip()
                if ":" in content:
                    key, value = content.split(":", 1)
                    metadata[key.strip()] = value.strip()
                continue
            tokens = line.split()
            if header is None:
                header = tokens
                continue
            if len(tokens) != len(header):
                break
            rows.append(dict(zip(header, tokens)))
    return metadata, rows


def configure_plot_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 13,
            "axes.labelsize": 15,
            "axes.titlesize": 16,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 0.8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
        }
    )


def annotate_log_bars(ax, bars, values):
    for bar, value in zip(bars, values):
        ax.annotate(
            f"{value:.1e}",
            xy=(bar.get_x() + bar.get_width() / 2.0, value),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=90,
            color="#252525",
        )


def plot_work_check(metadata, rows):
    case_labels = {
        "bending_y": r"Bending $y$" + "\n" + r"$\theta_z=0.25\pi$",
        "bending_z": r"Bending $z$" + "\n" + r"$\theta_y=0.20\pi$",
        "torsion": "Torsion\n" + r"$\theta_x=0.30\pi$",
        "displacement_y": r"Displacement $y$" + "\n" + r"$u_y=1$",
    }
    labels = [case_labels[row["kinematics_mode"]] for row in rows]
    errors = np.array([float(row["relative_error"]) for row in rows])
    tolerance = float(metadata["relative_tolerance"])

    fig, ax = plt.subplots(figsize=(10.8, 6.1))
    x = np.arange(len(rows))
    bars = ax.bar(
        x,
        errors,
        width=0.58,
        color=COLORS["beam_spline_mapper"],
        edgecolor="#222222",
        linewidth=0.7,
        zorder=3,
    )
    ax.axhline(
        tolerance,
        color="#B02A37",
        linestyle="--",
        linewidth=1.7,
        label=rf"Relative tolerance $={tolerance:.0e}$",
        zorder=4,
    )
    ax.set_yscale("log")
    ax.set_ylabel("Relative virtual-work error")
    ax.set_xticks(x, labels)
    ax.grid(axis="y", which="both", color="#D9D9D9", linewidth=0.7, zorder=0)
    ax.legend(frameon=False, loc="upper right")
    annotate_log_bars(ax, bars, errors)
    ax.set_ylim(min(errors) / 8.0, tolerance * 20.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(
        OUT
        / (
            "inverse_map_work_check__mesh_coarse_beam_surface_1x5"
            "__epsilon_1e-7__mapper_beam_spline_mapper"
            "__relative_virtual_work_error_by_case.png"
        ),
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def plot_mapper_comparison(rows):
    case_order = [
        ("rotation", "7.8539816339744828e-01"),
        ("rotation", "1.5707963267948966e+00"),
        ("displacement", "1.0000000000000000e+00"),
        ("torsion", "9.4247779607693793e-01"),
    ]
    case_labels = [
        r"Rotation" + "\n" + r"$\theta_z=0.25\pi$",
        r"Rotation" + "\n" + r"$\theta_z=0.50\pi$",
        r"Displacement" + "\n" + r"$u_y=1$",
        r"Torsion" + "\n" + r"$\theta_x=0.30\pi$",
    ]
    mapper_order = [
        "beam_spline_mapper",
        "beam_mapper_linear",
        "beam_mapper_corotation",
    ]

    indexed = {
        (row["kinematics_mode"], row["parameter_value"], row["mapper_type"]): row
        for row in rows
    }

    fig, ax = plt.subplots(figsize=(11.8, 6.5))
    x = np.arange(len(case_order))
    width = 0.23

    all_errors = []
    for mapper, offset in zip(mapper_order, [-width, 0.0, width]):
        errors = np.array(
            [
                float(indexed[(mode, parameter, mapper)]["relative_work_error"])
                for mode, parameter in case_order
            ]
        )
        all_errors.extend(errors)
        bars = ax.bar(
            x + offset,
            errors,
            width=width,
            color=COLORS[mapper],
            edgecolor="#222222",
            linewidth=0.6,
            label=LABELS[mapper],
            zorder=3,
        )
        annotate_log_bars(ax, bars, errors)

    ax.set_yscale("log")
    ax.set_ylabel("Relative virtual-work error")
    ax.set_xticks(x, case_labels)
    ax.grid(axis="y", which="both", color="#D9D9D9", linewidth=0.7, zorder=0)
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.set_ylim(min(all_errors) / 15.0, max(all_errors) * 25.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(
        OUT
        / (
            "inverse_map_mapper_comparison__mesh_coarse_beam_surface_1x5"
            "__epsilon_1e-7__beam_spline_vs_linear_vs_corotation"
            "__relative_virtual_work_error_by_case.png"
        ),
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def write_compact_csv(work_rows, comparison_rows):
    output = OUT / (
        "inverse_mapping__mesh_coarse_beam_surface_1x5__epsilon_1e-7"
        "__virtual_work_consistency_error_data.csv"
    )
    fields = [
        "study",
        "case",
        "mapper",
        "surface_virtual_work",
        "beam_virtual_work",
        "absolute_virtual_work_error",
        "relative_virtual_work_error",
    ]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in work_rows:
            writer.writerow(
                {
                    "study": "inverse_map_work_check",
                    "case": row["kinematics_mode"],
                    "mapper": "beam_spline_mapper",
                    "surface_virtual_work": row["surface_directional_work"],
                    "beam_virtual_work": row["beam_generalized_work"],
                    "absolute_virtual_work_error": row["absolute_error"],
                    "relative_virtual_work_error": row["relative_error"],
                }
            )
        for row in comparison_rows:
            writer.writerow(
                {
                    "study": "inverse_map_mapper_comparison",
                    "case": f'{row["kinematics_mode"]}_{row["parameter_value"]}',
                    "mapper": row["mapper_type"],
                    "surface_virtual_work": row["surface_work"],
                    "beam_virtual_work": row["beam_work"],
                    "absolute_virtual_work_error": row["absolute_work_error"],
                    "relative_virtual_work_error": row["relative_work_error"],
                }
            )


def main():
    configure_plot_style()
    OUT.mkdir(parents=True, exist_ok=True)
    work_metadata, work_rows = read_table(WORK_SUMMARY)
    _, comparison_rows = read_table(COMPARISON_SUMMARY)
    plot_work_check(work_metadata, work_rows)
    plot_mapper_comparison(comparison_rows)
    write_compact_csv(work_rows, comparison_rows)
    print(f"Created inverse-mapping virtual-work figures and CSV in: {OUT}")


if __name__ == "__main__":
    main()
