# Q1K

Preprocessing pipeline for EEG and eye-tracking data from the **Quebec 1000 Families** project. Processes raw EEG/ET recordings through BIDS conversion, artifact rejection, synchronization, segmentation, and epoch cleaning.

## Pipeline Overview

```
Raw EEG/ET data
     │
     ▼
┌─────────────┐
│ 1. q1k-init │  BIDS conversion + HTML reports
└─────┬───────┘
      ▼
┌──────────────────┐
│ 2. q1k-pylossless│  Artifact detection (PyLossless, via Slurm)
└─────┬────────────┘
      ▼
┌──────────────────┐
│ 3. q1k-sync-loss │  EEG/ET synchronization + lossless cleaning
└─────┬────────────┘
      ▼
┌────────────────┐
│ 4. q1k-segment │  Task-specific epoch segmentation
└─────┬──────────┘
      ▼
┌────────────────┐
│ 5. q1k-autorej │  AutoReject epoch cleaning (via Slurm)
└─────┬──────────┘
      ▼
  Grand averages (standalone marimo notebooks in grands/)
```

## Supported Tasks

| Code   | Task | Eye-Tracking |
|--------|------|:---:|
| RS     | Resting State | No |
| RSRio  | Resting State (Rio movie) | No |
| VEP    | Visual Evoked Potential | Yes |
| AEP  | Auditory Evoked Potential | No |
| GO   | Gap Overlap | Yes |
| PLR  | Pupillary Light Reflex | Yes |
| VS   | Visual Search | Yes |
| NSP  | Naturalistic Social Preference | Yes |
| TO   | Tone Oddball | No |

## Installation

```bash
git clone https://github.com/lina-usc/q1k.git
cd q1k
pip install -e .
```

On Narval (HPC), use the shared virtual environment:

```bash
source /project/def-emayada/q1k/experimental/q1k_env/bin/activate
```

### Dependencies

- [MNE-Python](https://mne.tools/) (>= 1.0)
- [MNE-BIDS](https://mne.tools/mne-bids/)
- [PyLossless](https://github.com/lina-usc/pylossless)
- [AutoReject](https://autoreject.github.io/)
- [marimo](https://marimo.io/) (for interactive notebooks)
- numpy, pandas, matplotlib, plotly

## Usage

Each pipeline stage has a CLI command. All commands accept `--subject` for a single subject or `--all` for batch processing.

### Stage 1: BIDS Initialization (`q1k-init`)

Converts raw EEG (.mff) and eye-tracking (.asc) files into BIDS format, extracts events, and generates per-subject report notebooks.

```bash
# Single subject
q1k-init --project-path /project/def-emayada/q1k/experimental/ \
         --task AEP --subject Q1K_HSJ_100123_F1 --site HSJ

# All unprocessed subjects
q1k-init --project-path /project/def-emayada/q1k/experimental/ \
         --task AEP --all --site HSJ
```

**Input:** `sourcedata/eeg/Q1K*/` (raw .mff files)
**Output:** BIDS-formatted files in the project root + per-subject marimo report notebooks.

### Stage 2: PyLossless Artifact Rejection (`q1k-pylossless`)

Runs the PyLossless pipeline for automated artifact detection. On Narval, use `--slurm` to submit cluster jobs.

```bash
# Submit Slurm jobs for all unprocessed subjects
q1k-pylossless --project-path /project/def-emayada/q1k/experimental/ \
               --task RS --all --slurm

# Run locally for a single subject
q1k-pylossless --project-path /project/def-emayada/q1k/experimental/ \
               --task RS --subject 0042P
```

**Input:** BIDS EEG files from Stage 1
**Output:** `derivatives/pylossless/` (ICA components, flagged channels/epochs)

### Stage 3: EEG/ET Sync + Lossless Cleaning (`q1k-sync-loss`)

Synchronizes eye-tracking data with EEG recordings (for ET-enabled tasks), applies PyLossless annotations (ICA cleaning, bad channel interpolation, average re-reference), and filters the data (1-90 Hz bandpass, 60 Hz notch).

```bash
q1k-sync-loss --project-path /project/def-emayada/q1k/experimental/ \
              --task VEP --all
```

**Prerequisites:** For ET-enabled tasks (VEP, GO, NSP, PLR, VS), `.ascii` eye-tracking files must be present in sourcedata (converted from Eyelink .edf using `edf2ascii`).

**Input:** `derivatives/pylossless/` + raw BIDS files
**Output:** `derivatives/pylossless/derivatives/sync_loss/`

### Stage 4: Epoch Segmentation (`q1k-segment`)

Segments continuous data into task-specific epochs. Each task has its own event selection logic and epoch time window.

```bash
q1k-segment --project-path /project/def-emayada/q1k/experimental/ \
            --task RS --all

# Use legacy postproc path instead of sync_loss
q1k-segment --project-path /project/def-emayada/q1k/experimental/ \
            --task RS --all --derivative-base postproc
```

| Task | Epoch Window | Conditions |
|------|-------------|------------|
| RS   | -0.2 to 0.8s | rest |
| VEP  | -1.0 to 2.0s | sv06_d (6Hz), sv15_d (15Hz) |
| AEP  | -1.0 to 2.0s | ae06_d (6Hz), ae40_d (40Hz) |
| GO   | -1.0 to 1.0s | dtbc_d (baseline), dtoc_d (overlap), dtgc_d (gap) |
| PLR  | -2.0 to 4.0s | plro_d |
| VS   | -1.0 to 1.0s | 5-item, 9-item, 13-item arrays |
| TO   | -1.0 to 1.0s | standard, deviant |

**Input:** `derivatives/pylossless/derivatives/sync_loss/`
**Output:** `derivatives/.../segment/epoch_fif_files/{TASK}/` (.fif epoch files)

### Stage 5: AutoReject (`q1k-autorej`)

Applies the AutoReject algorithm to automatically repair or reject bad epochs. On Narval, use `--slurm`.

```bash
q1k-autorej --project-path /project/def-emayada/q1k/experimental/ \
            --task RS --all --slurm
```

**Input:** Epoch files from Stage 4
**Output:** `derivatives/.../autorej/epoch_fif_files/{TASK}/`

### Pipeline Tracking (`q1k-tracking`)

Optional monitoring tool that cross-references REDCap demographic/task-completion exports with pipeline outputs at every stage. Generates per-task tracking CSVs, Sankey diagrams showing data flow, and data loss comparison reports against a manually curated SharePoint sheet.

This is **not a processing stage** — it can be run at any time and none of the processing stages depend on it.

```bash
# Generate tracking CSVs for all tasks
q1k-tracking --project-path /project/def-emayada/q1k/experimental/ \
             --redcap-dir /path/to/demographics_redcap/2025_12_12/

# Generate data loss comparison report against SharePoint
q1k-tracking --project-path /project/def-emayada/q1k/experimental/ \
             --redcap-dir /path/to/demographics_redcap/2025_12_12/ \
             --sharepoint /path/to/Sharepoint_dataloss_report.xlsx \
             --mni-upload-date 2025-12-08 --hsj-upload-date 2025-12-08
```

**Input:** REDCap CSV exports (Demographics + TaskCompletion) + BIDS filesystem
**Output:** Per-task tracking CSVs, multi-sheet Excel data loss report

Interactive notebooks are also available:

```bash
marimo edit q1k/notebooks/tracking_pipeline.py    # Pipeline tracking dashboard
marimo edit q1k/notebooks/tracking_data_loss.py   # Data loss comparison
```

### Grand Averages

Standalone marimo notebooks for group-level analysis are in `grands/`. Open interactively:

```bash
marimo edit grands/group_erp_tf_VEP.py
```

These compute grand average ERPs, time-frequency representations (Morlet wavelets), and cluster permutation statistics across participants.

## Per-Subject Report Notebooks

Each CLI stage generates a **marimo notebook per subject** as a processing log. These are saved alongside HTML exports:

```
reports/
├── init/
│   └── AEP/
│       ├── Q1K_HSJ_100123_F1_AEP_init.py    # Re-openable in marimo
│       └── Q1K_HSJ_100123_F1_AEP_init.html   # Quick HTML view
├── sync_loss/
│   └── VEP/
│       ├── 0042P_VEP_sync_loss.py
│       └── 0042P_VEP_sync_loss.html
└── segment/
    └── RS/
        ├── 0042P_RS_segment.py
        └── 0042P_RS_segment.html
```

To re-inspect a subject's results interactively:

```bash
marimo edit reports/segment/RS/0042P_RS_segment.py
```

## Package Structure

```
q1k/
├── config.py       # Constants: tasks, EOG channels, frequency bands
├── io.py           # Path utilities for the project directory structure
├── bids.py         # BIDS filename parsing and data writing
├── slurm.py        # Slurm job submission utilities
├── init/           # Stage 1: BIDS conversion
│   ├── tools.py    #   Event processing, EEG/ET combination
│   └── cli.py      #   q1k-init entry point
├── pylossless/     # Stage 2: Artifact rejection
│   ├── pipeline.py #   PyLossless execution
│   ├── config.yaml #   PyLossless configuration
│   └── cli.py      #   q1k-pylossless entry point
├── sync_loss/      # Stage 3: Sync + cleaning
│   ├── tools.py    #   apply_ll, eeg_et_combine
│   └── cli.py      #   q1k-sync-loss entry point
├── segment/        # Stage 4: Segmentation
│   ├── tasks.py    #   Task-specific segmentation functions
│   └── cli.py      #   q1k-segment entry point
├── autorej/        # Stage 5: AutoReject
│   ├── pipeline.py #   AutoReject execution
│   └── cli.py      #   q1k-autorej entry point
├── tracking/       # Pipeline tracking (optional)
│   ├── tools.py    #   REDCap loading, stage scanning, comparison
│   └── cli.py      #   q1k-tracking entry point
├── notebooks/      # Marimo template notebooks
│   ├── init_report.py
│   ├── sync_loss_report.py
│   ├── segment_*.py
│   ├── tracking_pipeline.py
│   └── tracking_data_loss.py
└── slurm/          # Slurm batch scripts
    ├── pylossless_job.sh
    └── autorej_job.sh
```

## Site Codes

| Code | Site |
|------|------|
| HSJ  | CHU Sainte-Justine (Hospital) |
| MHC  | McGill Health Centre (The Neuro) |
| NIM  | Neuroimaging |

## Authors

- James Desjardins (<jim.a.desjardins@gmail.com>)
- Gabriel Blanco Gomez (<gabriel.blancogomez@mail.mcgill.ca>)
- Christian O'Reilly (<christian.oreilly@sc.edu>)

Q1kpipeline extension, Integration, Validation and release
- Sweety Ramnani (sramnani@email.sc.edu)
## License

MIT
