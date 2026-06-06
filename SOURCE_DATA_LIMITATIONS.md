# Source data limitations

This repository is suitable for figure-level source-data review, but it is not a complete raw-data archive.

Known limitations:

- Raw FASTQ files are not included.
- Complete Seurat RDS, h5ad, h5, loom, and mtx count-matrix files were not found in the audited source package.
- Some source tables are exported figure-level analysis tables rather than original upstream objects.
- Figure 9 Prism files provide replicate-level values, but original animal IDs, ROI IDs, and full treatment metadata were not present.
- Some historical scripts contain workstation-local paths. See `metadata/local_path_audit.csv`.
- Large AI and zip assets were excluded from ordinary GitHub packaging and are listed in `metadata/large_assets_for_lfs_or_release_manifest.csv`.

Current Figure 9 status:

- CYP11A1, serum testosterone, DHH, GAS1, AMH, CLDN11, and SYCP3 values are now included.
- Raw Figure 9 image-file registry is included.
- Round-spermatids-per-tubule quantification remains missing if that claim is retained.
