# BugSigDB Access Benchmark

This file records a single-run benchmark of the non-R BugSigDB data access methods.

Generated at: `2026-04-20T20:15:10+00:00`

## Methodology

- one end-to-end run per access method
- files downloaded sequentially
- bytes read fully but not persisted to disk
- timings include manifest resolution for methods that require API discovery
- host environment: local workstation run from `bsdbng`, this repository's former
  name at the time of the run

## Results

| Method | Files | Declared size | Downloaded size | Total seconds | Throughput |
| --- | ---: | ---: | ---: | ---: | ---: |
| raw_csv_exports | 3 | 24.03 MiB | 24.03 MiB | 0.47 | 50.92 MiB/s |
| bugsigdbexports_full_dump | 1 | 26.52 MiB | 26.52 MiB | 0.29 | 91.30 MiB/s |
| bugsigdbexports_gmt_bundle | 15 | 44.02 MiB | 44.02 MiB | 0.66 | 66.94 MiB/s |
| zenodo_release_15272273 | 17 | 49.20 MiB | 49.20 MiB | 5.60 | 8.78 MiB/s |

## Notes

- `raw_csv_exports` is the source-faithful path and remains the recommended first ingestion surface.
- `bugsigdbexports_full_dump` is the convenience single-file path.
- `bugsigdbexports_gmt_bundle` is the heaviest path and is optimized for set-based analysis, not canonical ingestion.
- `zenodo_release_15272273` is the reproducible pinned-release path.
