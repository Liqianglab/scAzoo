#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
FIGURE_FILES = PACKAGE / "_current_figure_files"
OUT = PACKAGE / "figure9_source_data"
RAW = OUT / "raw_support"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)

    image_rows: list[dict[str, object]] = []
    for src in sorted(FIGURE_FILES.glob("Figure 9_20260602.*")):
        dst = RAW / src.name
        shutil.copy2(src, dst)
        image_rows.append(
            {
                "figure": "Figure 9",
                "source_file": str(src),
                "copied_to": str(dst.relative_to(PACKAGE)),
                "file_type": src.suffix.lstrip("."),
                "size_bytes": src.stat().st_size,
                "sha256": sha256(src),
                "role": "current composed figure file; not numeric source data",
            }
        )

    if not image_rows:
        image_rows.append(
            {
                "figure": "Figure 9",
                "source_file": "NOT_FOUND",
                "copied_to": "",
                "file_type": "",
                "size_bytes": "",
                "sha256": "",
                "role": "current Figure 9 file was not found in _current_figure_files",
            }
        )

    write_csv(
        OUT / "Fig9_image_file_registry.csv",
        image_rows,
        ["figure", "source_file", "copied_to", "file_type", "size_bytes", "sha256", "role"],
    )

    panel_rows = [
        {
            "panel": "a",
            "displayed_content": "Animal experimental design: Vehicle, EDS, EDS + testosterone rescue time course",
            "required_raw_source": "animal metadata, treatment dates, dose, route, sampling day, exclusion criteria",
            "current_workspace_status": "figure_image_available_only",
            "notes": "Protocol schematic is available in the composed figure, but raw animal metadata table was not located.",
        },
        {
            "panel": "b",
            "displayed_content": "CYP11A1 immunofluorescence and CYP11A1+DAPI interstitial cell density",
            "required_raw_source": "raw microscopy files, ROI annotations, cell-count table, animal/sample IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Do not infer values from the rendered bar plot; add measured count data.",
        },
        {
            "panel": "c",
            "displayed_content": "Serum testosterone concentration",
            "required_raw_source": "hormone assay output table with animal IDs, group, value, unit, batch",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Needed to support androgen rescue claim.",
        },
        {
            "panel": "d",
            "displayed_content": "WT1/DHH and CYP11A1/GAS1 immunofluorescence representative images",
            "required_raw_source": "raw microscopy files and channel metadata",
            "current_workspace_status": "raw_images_not_located",
            "notes": "Representative images are embedded in the composed figure only.",
        },
        {
            "panel": "e",
            "displayed_content": "GAS1 relative fluorescence intensity",
            "required_raw_source": "per-image or per-ROI GAS1 intensity measurements with animal/sample IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Needed for HH receptor mismatch claim.",
        },
        {
            "panel": "f",
            "displayed_content": "DHH relative fluorescence intensity",
            "required_raw_source": "per-image or per-ROI DHH intensity measurements with animal/sample IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Needed for Sertoli-side DHH retained/increased claim.",
        },
        {
            "panel": "g",
            "displayed_content": "CLDN11 and AMH immunofluorescence representative images",
            "required_raw_source": "raw microscopy files and channel metadata",
            "current_workspace_status": "raw_images_not_located",
            "notes": "Representative images are embedded in the composed figure only.",
        },
        {
            "panel": "h",
            "displayed_content": "CLDN11 continuity index",
            "required_raw_source": "junction-continuity measurements with image/ROI/animal IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Needed for BTB damage/rescue claim.",
        },
        {
            "panel": "i",
            "displayed_content": "AMH relative fluorescence intensity",
            "required_raw_source": "per-image or per-ROI AMH intensity measurements with animal/sample IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Needed for Sertoli dysfunction phenotype.",
        },
        {
            "panel": "j",
            "displayed_content": "SYCP3 immunofluorescence representative images and relative intensity",
            "required_raw_source": "raw microscopy files and SYCP3 intensity measurements with animal/sample IDs",
            "current_workspace_status": "raw_quantification_not_located",
            "notes": "Figure label should be checked for SYCP3 spelling before submission.",
        },
        {
            "panel": "k",
            "displayed_content": "Mechanistic summary model",
            "required_raw_source": "derived from panels b-j and prior figures",
            "current_workspace_status": "schematic_only",
            "notes": "No numeric source data expected, but model claims should cite supporting panels.",
        },
    ]

    write_csv(
        OUT / "Fig9_panels_data_availability.csv",
        panel_rows,
        ["panel", "displayed_content", "required_raw_source", "current_workspace_status", "notes"],
    )

    template_rows: list[dict[str, object]] = []
    assays = [
        ("b", "CYP11A1+DAPI interstitial cells", "cells/mm2"),
        ("c", "Serum testosterone", "ng/mL"),
        ("e", "GAS1 relative fluorescence intensity", "relative_intensity"),
        ("f", "DHH relative fluorescence intensity", "relative_intensity"),
        ("h", "CLDN11 continuity index", "percent_or_index"),
        ("i", "AMH relative fluorescence intensity", "relative_intensity"),
        ("j", "SYCP3 relative fluorescence intensity", "relative_intensity"),
    ]
    for panel, assay, unit in assays:
        template_rows.append(
            {
                "panel": panel,
                "assay": assay,
                "group": "Vehicle|EDS|EDS+T",
                "animal_id": "",
                "sample_id": "",
                "image_or_assay_id": "",
                "roi_id": "",
                "value": "",
                "unit": unit,
                "biological_replicate": "",
                "technical_replicate": "",
                "source_file": "",
                "notes": "",
            }
        )

    write_csv(
        OUT / "Fig9_raw_quantification_template.csv",
        template_rows,
        [
            "panel",
            "assay",
            "group",
            "animal_id",
            "sample_id",
            "image_or_assay_id",
            "roi_id",
            "value",
            "unit",
            "biological_replicate",
            "technical_replicate",
            "source_file",
            "notes",
        ],
    )

    mapping_rows = [
        {
            "figure_panel": f"Fig9{row['panel']}",
            "source_data_file": "Fig9_raw_quantification_template.csv"
            if row["current_workspace_status"].startswith("raw_quantification")
            else "Fig9_image_file_registry.csv",
            "code_file": "figure9_code/build_figure9_source_data.py",
            "status": row["current_workspace_status"],
            "comment": row["notes"],
        }
        for row in panel_rows
    ]
    write_csv(
        OUT / "Fig9_file_mapping.csv",
        mapping_rows,
        ["figure_panel", "source_data_file", "code_file", "status", "comment"],
    )

    readme = OUT / "README.txt"
    readme.write_text(
        "Figure 9 source data package\n\n"
        "This directory was created from the files currently present in the workspace. "
        "Only the composed Figure 9 PNG was located. Numeric animal validation raw data "
        "and raw microscopy files were not present in the current workspace, so they are "
        "not reconstructed from the rendered figure.\n\n"
        "Generated files:\n"
        "- Fig9_image_file_registry.csv: current composed figure file registry and checksum.\n"
        "- Fig9_panels_data_availability.csv: panel-level source-data availability and required raw inputs.\n"
        "- Fig9_raw_quantification_template.csv: template for animal-level and ROI-level measurements.\n"
        "- Fig9_file_mapping.csv: panel-to-source-data/code mapping.\n"
        "- raw_support/: copied current Figure 9 composed image.\n\n"
        "To complete this package, add raw hormone assay tables, microscopy image metadata, "
        "ROI-level fluorescence measurements, and cell-count/continuity-index tables to raw_support/ "
        "and replace the template rows with measured values.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
