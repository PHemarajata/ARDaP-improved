## ARDaP release v2025.10.24

Date: 2025-10-24

Summary
- Fix snpEff invocation in helper scripts so snpEff runs from the dedicated Java-21 conda env (`ardap_snpeff21`) by default. This resolves class-version issues observed with snpEff and Java versions.
- Prefer `mamba run -n ardap_snpeff21 snpEff` (with a `conda run` fallback and final `snpEff` fallback) in `bin/Masked_alignment.sh`, `bin/SNP_matrix.sh`, and `bin/VariantCalling.sh`.
- Added `mosdepth` to the GATK env manifest earlier (environment files) to fix deduplicate coverage computation.
- Cleaned repository by adding `.gitignore` entries for common pipeline outputs and transient files (Outputs/, SmokeData/, smoke reports, work caches).

Notes
- Smoke tests completed successfully after the fixes; annotated VCFs and reports were produced under a smoke output folder locally during testing but these folders are intentionally ignored from the repo.
- This release creates a local annotated tag; push it to the remote and create a GitHub release if you want a published release artifact.

How to publish
1. Push the branch and the tag to the remote:

   git push origin main
   git push origin v2025.10.24

2. Create a GitHub release (web UI) for tag `v2025.10.24` and paste these notes, or use the GitHub CLI:

   gh release create v2025.10.24 -t "v2025.10.24" -F RELEASE_NOTES_v2025.10.24.md

If you'd like, I can push the branch and tag and create the GitHub release for you (requires network access and write permissions).
