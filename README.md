# Figure source data and analysis code

This repository contains figure-level source data, supplementary tables, and analysis scripts for the manuscript. It is intended to support manuscript review and source-data inspection.

## Contents

- `source_data/`: figure-level and supplementary-figure source tables.
- `code/`: scripts used to generate or assemble figure-level source data and selected plots.
- `figures/final_png/`: final manuscript figure PNGs.
- `supplementary_tables/`: supplementary tables, including the updated Figure 9 animal-validation source data.
- `metadata/`: figure-to-source mapping, panel manifest, file manifest, excluded-file manifest, and path audit.
- `environment/`: detected Python and R package lists.

## Data availability

Processed single-cell analysis data have been deposited in OMIX under accession number `OMIX017477`.

This repository does not include raw FASTQ files. It also does not claim full end-to-end reprocessing from raw sequencing reads. The repository is organized for figure-level reproducibility and reviewer inspection of source tables, selected code, and supporting metadata.

## Figure 9 animal validation

Figure 9 animal-validation source data were updated on 2026-06-06.

Key files:

- `source_data/figure9_source_data/Fig9_animal_level_source_data_20260606.csv`
- `source_data/figure9_source_data/animal_validation_20260606/Fig9_animal_level_source_data_20260606.csv`
- `source_data/figure9_source_data/animal_validation_20260606/Fig9_prism_source_files_manifest_20260606.csv`
- `source_data/figure9_source_data/animal_validation_20260606/Fig9_raw_image_file_registry_20260606.csv`
- `supplementary_tables/Supplementary_Table_11_EDS_Animal_Validation_Source_Data.xlsx`
- `supplementary_tables/Supplementary_Tables_11_12_13.xlsx`

The Figure 9 values were extracted from selected GraphPad Prism `pzfx` tables. Original animal IDs and ROI IDs were not present in the Prism files, so deidentified replicate IDs were generated from group and row order. Remaining limitations are listed in `source_data/figure9_source_data/animal_validation_20260606/Fig9_remaining_gaps_20260606.csv`.

## Large files

Standard GitHub upload rejects individual files larger than 100 MB. Therefore, AI files, zip archives, and other large files were not included in this clean repository. They are listed in:

- `metadata/large_assets_for_lfs_or_release_manifest.csv`
- `metadata/excluded_files_manifest.csv`

Those files remain in the full local source package and should be handled by Git LFS, a GitHub release asset, Zenodo, or the institutional data repository if they need to be distributed.

## Code notes

Some historical scripts still contain local path references because they were originally used on the analysis workstation. The audit is provided in `metadata/local_path_audit.csv`. Source tables are included so reviewers can inspect the data even where a historical script is not directly runnable after cloning.

Recommended read order:

1. `metadata/figure_index.csv`
2. `metadata/latest_panel_manifest.csv`
3. `metadata/current_figure_to_source_data_mapping.csv`
4. `metadata/source_data_gaps_and_actions.csv`
5. `metadata/repository_manifest.csv`

## License

No reuse license has been selected yet. Add a license before making this repository fully public if broad reuse is intended.
