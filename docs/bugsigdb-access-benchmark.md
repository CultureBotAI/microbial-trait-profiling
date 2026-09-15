# BugSigDB Access Benchmark

This file records a single-run benchmark of the non-R BugSigDB data access methods.

Generated at: `2026-08-26T19:42:10+00:00`

## Methodology

- one end-to-end run per access method
- files downloaded sequentially
- bytes read fully but not persisted to disk
- timings include manifest resolution for methods that require API discovery
- host environment: local workstation run from `bsdbng`

## Results

| Method | Files | Declared size | Downloaded size | Total seconds | Throughput |
| --- | ---: | ---: | ---: | ---: | ---: |
| raw_csv_exports | 3 | 27.44 MiB | 27.44 MiB | 5.41 | 5.07 MiB/s |
| bugsigdbexports_v1_3_1_full_dump | 1 | 15.49 MiB | 15.49 MiB | 1.44 | 10.72 MiB/s |
| bugsigdbexports_v1_3_1_gmt_bundle | 15 | 26.03 MiB | 26.03 MiB | 7.38 | 3.53 MiB/s |
| zenodo_release_v1_3_1 | 17 | 41.53 MiB | 41.53 MiB | 49.52 | 0.84 MiB/s |

## Notes

- `raw_csv_exports` is the live source-faithful path and remains the preferred surface once separate-table ingestion is implemented.
- `bugsigdbexports_v1_3_1_full_dump` is the citable convenience single-file path currently used by bsdbng.
- `bugsigdbexports_v1_3_1_gmt_bundle` is the heaviest path and is optimized for set-based analysis, not canonical ingestion.
- `zenodo_release_v1_3_1` is the reproducible pinned-release path (10.5281/zenodo.19740674).
