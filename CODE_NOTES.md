# Code notes

The code folders preserve scripts used during figure-level analysis and source-data preparation.

The repository should be cited as supporting figure-level reproducibility rather than full end-to-end reprocessing from raw sequencing reads. Several scripts depend on already exported source tables or historical workstation paths. These path references are not hidden; they are audited in `metadata/local_path_audit.csv`.

Suggested use:

1. Inspect `metadata/latest_panel_manifest.csv` to identify source files and code for each panel.
2. Open the relevant `source_data/<figure>_source_data/` directory.
3. Review the matching `code/<figure>_code/` scripts.
4. Use `environment/requirements.txt` and `environment/R_packages.txt` as dependency hints rather than a locked environment.
