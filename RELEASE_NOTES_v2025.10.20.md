# Release v2025.10.20 — Reproducibility & Resumability improvements (ARDaP)

This release focuses on reproducibility, resumability and reliable report generation. It stages helper scripts in Nextflow so script edits invalidate process caches, makes helper scripts produce per-sample outputs (no cross-sample appends), fixes the HTML/CSV mismatch, and adds small monitoring/collection helpers used during testing.

## Highlights

- Make pipeline outputs configurable:
  - Add support for `--outdir` (set output publish directory) in `main.nf`.
  - Default remains `Outputs/` when `--outdir` is not provided.

- Resume-safety and reproducibility:
  - Helper scripts (reporting, SQL-driven steps) are staged as inputs to Nextflow processes so changes to helpers change process hashes.
  - Refactor reporting helpers to write per-sample files (overwrite) instead of appending a shared file.

- Deterministic reports and CSVs:
  - `bin/AbR_reports.sh` now writes a deterministic canonical CSV: header `ID,Class,Drug,Status,Level,Details`.
  - Ensures `Status`/`Level`/`Details` are populated (placeholders when upstream data missing) so `Report_html.sh` reliably marks resistant rows.
  - Deduplication and sanitization of AbR outputs to avoid duplication/noisy tokens.

- Propagated fixes:
  - Same per-sample, overwrite and CSV-generator behavior applied to `AbR_reports copy.sh`, `AbR_reports_mix.sh`, and `AbR_reports_mix copy.sh`.

- Operational/QA additions (non-invasive):
  - Small watcher script used during testing: `/tmp/ardap_report_watcher.sh` (polls `Outputs/AbR_reports` and `Reports/data/patientDrugSusceptibilityData.csv`, logs to `/tmp/ardap_report_events.log`).
  - Snapshot/collection routine used to create a timestamped archive of produced reports (e.g., `Outputs/collected_YYYYMMDD_HHMMSS.tar.gz`).

- Run-time/resume behavior validated:
  - Smoke-tested on sample `SRR2102060` — produced a working `patientDrugSusceptibilityData.csv` with populated `Status` fields and a rendered `_report.html` showing resistant rows.

## Files changed (representative)

- `main.nf` — add `params.outdir` support; stage helper scripts as process inputs; wire `publishDir` to `params.outdir`.
- `bin/AbR_reports.sh` — refactored to per-sample generation, deterministic CSV header, dedupe AbR outputs, placeholders where upstream files missing, write `${seq}.drug.table.txt` (overwrite).
- `bin/AbR_reports_mix.sh`, `bin/AbR_reports copy.sh`, `bin/AbR_reports_mix copy.sh` — aligned behavior with the main reporting script.
- `bin/Report_html.sh` — works with deterministic CSV; no API/semantic changes besides compatibility with generated CSV.

## Upgrade / migration notes

- Nextflow version: this pipeline is DSL1; we recommend running with Nextflow v20.x (the environment used during testing was `nextflow 20.04.1`). If you upgrade to a newer Nextflow, test cautiously.

- How to run with custom output directory:

  ```bash
  nextflow run main.nf \
    --fastq '/path/to/fastq/*_{R1,R2}*.fastq.gz' \
    --outdir '/path/to/my_outputs' \
    --size 1000000 \
    -resume -with-trace
  ```

  If you do not set `--outdir`, the pipeline uses `Outputs/` in the repo root.

  To silence the Nextflow WARN about `outdir`, add a default in `nextflow.config` or `main.nf`:

  ```groovy
  params.outdir = params.outdir ?: "${baseDir}/Outputs"
  ```

- Resume behavior:
  - Because helpers are now staged, re-running with `-resume` will pick up processed tasks unless helper scripts changed (which will cause affected processes to re-run — intentional).
  - Per-sample files are used and overwritten, so resume is safe and avoids cross-sample contamination.

## Breaking changes / compatibility

- Any external scripts that relied on a shared global `drug.table.txt` being appended-to across samples should be updated. This release switches to per-sample `*.drug.table.txt` and deterministic CSV generation. (A repo search during the work found ~20 legacy references; they may require manual review.)

- The CSV header/order is now deterministic — if you have downstream parsers that expect other orders/columns, update them to accept `ID,Class,Drug,Status,Level,Details`.

## Known issues & TODOs

- Nextflow WARN: If you don't set `params.outdir` in `nextflow.config`, Nextflow warns about undefined `outdir`. This is harmless; you can silence it by adding a default as shown above.

- Remaining scripts referencing global `drug.table.txt` may need to be migrated (optional but recommended). I can patch these consecutively if you’d like.

- I added a temporary watcher and collection helper for testing/monitoring; these are not required for production runs and can be removed if you prefer no extra files under /tmp.

## Release QA summary

- Smoke tests passed for sample `SRR2102060`:
  - Report: `SRR2102060_report_correct.html` rendered expected resistant rows.
  - CSV: `Reports/data/patientDrugSusceptibilityData.csv` contains populated `Status/Level/Details`.

- Per-sample outputs were produced and archived during testing:
  - Example archive: `Outputs/collected_20251020_213743.tar.gz` (contains the above reports and CSV).

- Stale Nextflow session lock handling: during the run a LevelDB `LOCK` was found and remedied by terminating the process that held the lock; documented as an operational note (not a code change).

## Release checklist (mapping requirements → status)

- [x] Make `publishDir` specifiable via `--outdir` — Done (`main.nf`).
- [x] Stage helper scripts so process hashes change when scripts change — Done.
- [x] Refactor helpers to be per-sample for safe resume — Done (`bin/AbR_reports.sh` + variants).
- [x] Fix HTML/CSV mismatch (Status missing) — Done (deterministic CSV generator).
- [x] Propagate fixes across report script variants — Done for main copies and mixes.
- [x] Smoke test using sample `SRR2102060` — Done; working HTML & CSV observed.

## Next steps / suggestions

- (Optional) Patch remaining legacy scripts referencing shared `drug.table.txt` (I can do this).
- (Optional) Extract the deterministic CSV generator into a shared helper script to avoid duplication across `AbR_reports*` scripts.
- (Optional) Create a formal GitHub release tag (I can draft the tag and push the release if you want).
- Update `README.md` to document `--outdir` and resume behavior for users.

If you want, I can:
- Create the Git tag and open a GitHub release using this release note as the body, or
- Produce a short PR that lists exact file diffs and test artifacts, or
- Implement the remaining cleanup (patch legacy script references and extract the CSV generator).
