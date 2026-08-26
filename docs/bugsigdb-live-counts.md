# BugSigDB Live Counts

This note records the current raw BugSigDB export counts used for comparison
with the pinned `v1.3.1` release workflow.

As of 2026-08-26, using the raw CSV tables linked from `Help:Export`:

- studies: 2,191 rows
  - 2,169 `Complete`
  - 22 `Incomplete`
- experiments: 9,140 rows
  - 7,210 `Homo sapiens`
  - 938 `Mus musculus`
  - 237 `Rattus norvegicus`
- signatures: 14,474 rows
  - 14,468 `Complete`
  - 4 `Incomplete`
  - 2 blank `State`

The citable `BugSigDBExports` `v1.3.1` full dump used by the regenerated
branch artifacts contains 7,425 complete rows and 1,246 parsed study IDs. The
bsdbng ingest writes 1,245 usable study YAML files because one parsed study has
no valid experiment/signature content after ingestion filtering.

## Source Mapping

- `studies` count: row count of `studies.csv`
- `experiments` count: row count of `experiments.csv`
- `signatures` count: row count of `signatures.csv`
- citable branch input: `waldronlab/BugSigDBExports` tag `v1.3.1`
  (`10.5281/zenodo.19740674`)

The live raw CSV exports are fresher than `v1.3.1`, but they are not the primary
ingestion source for this branch yet because bsdbng still lacks complete
separate-table join support.

## Sources

- BugSigDB `Help:Export`: https://bugsigdb.org/Help:Export
- BugSigDB main page: https://bugsigdb.org/Main_Page
- BugSigDBExports `v1.3.1`: https://github.com/waldronlab/BugSigDBExports/tree/v1.3.1
- Zenodo record: https://zenodo.org/records/19740674
