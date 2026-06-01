# VDJ ↔ GEX QC Dashboard

Interactive Shiny dashboard for B-cell receptor (BCR) VDJ ↔ GEX barcode matching, clonotype analysis, barcode recovery assessment, and VDJ quality control.

<p align="center">
  <img src="Screenshot.png" width="1200">
</p>

<p align="center">
  <i>Figure 1. Example dashboard showing VDJ ↔ GEX barcode recovery, QC metrics, clonotype statistics, and Sankey visualization.</i>
</p>

---

## Overview

The **VDJ ↔ GEX QC Dashboard** is a `{golem}`-based Shiny application designed to evaluate barcode recovery between single-cell V(D)J sequencing and gene expression (GEX) datasets.

The application helps researchers understand:

* How many VDJ cells successfully map to GEX
* How many productive BCRs were recovered
* Where barcode losses occur during processing
* Clonotype expansion patterns
* V gene and J gene usage
* Sample-specific and lane-specific QC metrics

The dashboard was originally developed for large-scale immune repertoire studies involving PBMC, vaccine-response, infection, and antigen-specific B-cell datasets.

---

## Features

### VDJ Quality Control

* Total VDJ contigs
* Productive full-length contigs
* Unique VDJ barcodes
* Paired heavy/light chain recovery
* IGH, IGK, and IGL chain statistics
* Duplicate heavy chain detection
* Duplicate light chain detection
* Cells with more than two productive chains

### Clonotype Analysis

* Unique clonotypes
* Expanded clonotypes
* Singleton clonotypes
* Largest clonotype size
* Clone size distributions
* Top clonotype visualization

### VDJ ↔ GEX Integration

* Barcode matching between VDJ and GEX datasets
* Config-aware barcode matching
* Recovery statistics
* Missing barcode identification
* Mapping percentages

### Barcode Loss Tracking

When both merged and filtered Seurat objects are provided, the dashboard can track:

```text
Paired VDJ Cells
        ↓
Merged Seurat Object
        ↓
Filtered Seurat Object
```

This allows users to determine where cells are lost during processing and quality control.

### V/J Gene Usage

* V gene usage frequencies
* J gene usage frequencies
* Chain-specific usage summaries

### Visualization

Interactive Sankey diagram showing:

```text
VDJ Contigs
    ↓
VDJ Barcodes
    ↓
Paired BCR Cells
    ↓
Mapped to GEX
```

### Exportable Reports

* Flow summary CSV
* Barcode loss summary CSV
* VDJ QC summary CSV
* Per-config summary CSV
* Clone size summary CSV
* V gene usage CSV
* J gene usage CSV
* Interactive Sankey HTML

---

## Installation

### Install Dependencies

```r
install.packages(c(
  "shiny",
  "golem",
  "dplyr",
  "readr",
  "purrr",
  "stringr",
  "tibble",
  "DT",
  "networkD3",
  "htmlwidgets",
  "ggplot2",
  "plotly",
  "Seurat"
))
```

### Install Package

```r
devtools::install_github("foocheung/vdjgexqc")
```

Or install locally:

```r
devtools::install_local("vdjgexqc")
```

---

## Running the Application

```r
library(vdjgexqc)

run_app()
```

For large Seurat objects:

```r
options(shiny.maxRequestSize = 5000 * 1024^2)

run_app()
```

This increases the upload limit to approximately 5 GB.

---

## Input Files

### Required

#### Filtered Seurat Object

Upload a filtered/final Seurat object:

```r
saveRDS(pbmc_filt, "pbmc_filt.rds")
```

#### VDJ Files

Upload one or more Cell Ranger VDJ annotation files:

```text
multi_config_1_filtered_contig_annotations.csv
multi_config_2_filtered_contig_annotations.csv
multi_config_3_filtered_contig_annotations.csv
multi_config_4_filtered_contig_annotations.csv
multi_config_5_filtered_contig_annotations.csv
```

---

### Optional

#### Merged or Pre-QC Seurat Object

```r
saveRDS(pbmc_merged, "pbmc_merged.rds")
```

Providing both merged and filtered objects enables barcode loss tracking.

---

## Expected Cell Naming Convention

By default, the application assumes Seurat cell names follow:

```text
multi_config_1_AAACCTGAGTAACTTG-1
```

and VDJ files follow:

```text
multi_config_1_filtered_contig_annotations.csv
```

Custom barcode parsing rules can be supplied through the dashboard sidebar.

---

## Workflow

### Step 1

Upload a filtered Seurat object.

### Step 2

Optionally upload a merged/pre-QC Seurat object.

### Step 3

Upload one or more VDJ annotation files.

### Step 4

Review barcode parsing settings.

### Step 5

Click:

```text
Run VDJ ↔ GEX QC
```

### Step 6

Explore:

* Sankey visualization
* Flow summary
* Barcode loss statistics
* VDJ QC metrics
* Per-config summaries
* Clone size distributions
* V gene usage
* J gene usage
* Matched barcode tables

### Step 7

Export tables and reports.

---

## Example Outputs

### Flow Summary

| Metric              | Example |
| ------------------- | ------: |
| Total VDJ Contigs   |  25,872 |
| Unique VDJ Barcodes |  12,038 |
| Paired BCR Cells    |  11,409 |
| Mapped to GEX       |   7,987 |
| Percent Mapped      |   70.0% |

### VDJ QC Summary

| Metric            | Example |
| ----------------- | ------: |
| IGH Contigs       |  12,518 |
| IGK Contigs       |   7,683 |
| IGL Contigs       |   5,671 |
| Unique Clonotypes |  10,694 |
| Largest Clone     |      34 |

---

## Applications

This dashboard is useful for:

* Single-cell immune repertoire studies
* B-cell receptor recovery QC
* PBMC datasets with VDJ sequencing
* Antigen-specific B-cell studies
* Vaccine-response studies
* Dengue virus studies
* HIV studies
* COVID-19 immune profiling
* Longitudinal immune monitoring

---

## Future Development

Planned additions include:

* B-cell subset integration
* Clonotype overlap analysis
* Subject-level summaries
* UMAP overlays for expanded clones
* Repertoire diversity metrics
* Differential clonotype analysis
* Interactive clonotype exploration

---

## Citation

If you use this dashboard in a publication, please cite the GitHub repository and include the version used for analysis.

---

## License

MIT License

---

## Contact

Questions, bug reports, and feature requests are welcome through GitHub Issues.
