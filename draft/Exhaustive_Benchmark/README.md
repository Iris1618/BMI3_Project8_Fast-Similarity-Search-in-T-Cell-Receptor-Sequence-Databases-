# Exhaustive search and benchmark

This folder adds a length-filtered exhaustive-search baseline, detailed result validation and preliminary benchmarking to Project 8. It uses the repository's existing `verification_basic.py`, `verification.py`, and `Inverted_Index/inverted_index.py`; those files are not duplicated here.

## Repository placement

Place the entire `Exhaustive_Benchmark` folder at the repository root, next to `verification_basic.py`, `verification.py`, and `Inverted_Index/`. Run the commands below **from the repository root**. Python 3.10 or later is required by the current repository's optimized verifier.

```text
repository root/
  verification_basic.py             # already in repository
  verification.py                   # already in repository
  Inverted_Index/                   # already in repository
  Exhaustive_Benchmark/             # this folder
    __init__.py
    baseline.py
    compare.py
    test_comparison.py
    README.md
```

## Data and input

Data is not included in this folder. Obtain the team's prepared data package and use your local paths. TSV files must be UTF-8 tab-separated, with the following columns:

- Reference: `reference_id`, `sequence`.
- Queries: `query_id`, `sequence`; other columns may be present.
- Optional truth: the team's `truth_matches.tsv` and `truth_summary.tsv`, plus its `reference_set` name.

`d` is the inclusive maximum global Levenshtein distance. Matches have exact distances; insertion, deletion and substitution each cost 1. Annotations are not part of the benchmark or the ID/distance comparison. The baseline's `search_one` can attach annotations if a caller passes an ID-to-records dictionary.

`read_sequences` is the index module's existing TSV reader. It requires unique, nonempty IDs and nonempty uppercase sequences using the 20 standard amino acids. The reference and query inputs are loaded once and shared by all methods.

## 1. Run the self-contained tests

```powershell
python -m unittest -v Exhaustive_Benchmark.test_comparison
```

These small tests need no VDJdb package. They include the length-filtered baseline, both seed strategies, fallback, no matches, incorrect distances, a deliberately broken verifier, truth coverage and benchmark output.

## 2. Validate one threshold

The paths below are examples: replace `DATA` with the actual location of the prepared data package. Output must be a new file, and its parent folder must already exist.

```powershell
python -m Exhaustive_Benchmark.compare validate --reference "DATA/real/reference_n1000.tsv" --queries "DATA/real/queries.tsv" --d 1 --k 3 --truth "DATA/real/truth_matches.tsv" --truth-summary "DATA/real/truth_summary.tsv" --reference-set real_n1000 --output "validation_d1.csv"
```

For each query the CSV records:

1. `indexed_basic`: indexed search with basic verification compared with exhaustive + basic verification.
2. `indexed_optimized`: indexed search with optimized verification compared with exhaustive + basic verification.
3. `exhaustive_vs_supplied_truth`: exhaustive search compared with the team's supplied truth, only when all three truth options are supplied.

Columns include query ID and sequence, `d`, expected and observed match counts, `match`, missing IDs, unexpected IDs, wrong distances, candidate count, length-eligible count, and index/fallback mode. Empty results still produce rows. A distance discrepancy is shown as `[reference_id, expected_distance, observed_distance]` inside `wrong_distances`. The program returns a nonzero exit code if any comparison fails. CSV files are UTF-8 with BOM for convenient Excel opening.

The truth comparison checks query coverage and reference counts against the supplied summary and checks truth row counts, but it does not independently prove that the input file has the intended exact content. Keep the team's data version and checksums together with your experiment records. The validation compares reference IDs and distances; it does not validate biological annotations.

Repeat with `--d 0`, `--d 2`, and `--d 3`, using a new output path each time. `--limit N` selects the first N queries for a quick pilot; omit it for all queries. `--strategy prefix` selects the other index seed-placement strategy.

## 3. Run a preliminary benchmark

```powershell
python -m Exhaustive_Benchmark.compare benchmark --reference "DATA/real/reference_n1000.tsv" --queries "DATA/real/queries.tsv" --d 1 --k 3 --limit 10 --repeats 5 --output "benchmark_d1.csv" --metadata "benchmark_d1_metadata.json"
```

The benchmark first checks that the indexed methods agree with the exhaustive baseline. It then reports:

- `exhaustive_basic`: length-filtered exhaustive search with the basic verifier.
- `indexed_basic`: candidate generation and basic verification, using an index built before timing.
- `indexed_optimized`: candidate generation and optimized verification, using the same index.
- `candidate_basic` and `candidate_optimized`: verifier time per pair on the exact same precomputed candidate pairs. These measurements exclude candidate generation.

Whole-search timing is in milliseconds per query and includes result construction. Verifier-only timing is in milliseconds per candidate pair. Input loading, correctness checks and file output are excluded. Index building is timed separately and saved in the metadata JSON. If no candidate pairs exist, per-pair timing is blank, not zero.

The methods are warmed up and their timing order is rotated. The CSV retains each batch average, the mean, median and standard error across repeated batch averages. This standard error describes run-to-run timing variability for the same query batch, not biological variation or variation across independent queries. The metadata includes the exact query IDs, candidate and length-eligible totals, fallback count, index build samples, Python version and platform. Outputs never intentionally overwrite an existing file.

Index build cost should be considered separately and, if relevant, amortized over repeated queries. Do not treat this pilot as a final performance claim.

## 4. Measure Python allocation peaks separately

```powershell
python -m Exhaustive_Benchmark.compare memory --reference "DATA/real/reference_n1000.tsv" --queries "DATA/real/queries.tsv" --d 1 --k 3 --limit 10 --output "memory_d1.csv"
```

This optional run uses Python's `tracemalloc` for each method separately. Input reference and query data are loaded before tracing. Indexed methods then build an index and search the query batch; the exhaustive method searches the same query batch without building an index. The CSV records the highest **additional traced Python allocations** during each run. This is not full process RSS, does not include preloaded input data or all native allocations, and should not be reported as total program memory. These runs are separate from timing because memory tracing adds overhead. For final results, supplement them with a reproducible process-level peak RSS measurement if required.

The supplied real query file contains exact-copy queries; discuss adding held-out real queries and use the team's synthetic edge-case and motif datasets for wider coverage. Simulated insertions/deletions are algorithm stress tests and do not automatically model biological V/D/J junctions.
