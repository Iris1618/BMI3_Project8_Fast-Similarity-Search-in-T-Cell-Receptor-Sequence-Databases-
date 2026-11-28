# TCR Similarity Search

Find T-cell receptor sequences within a given edit-distance threshold using an inverted index and dynamic programming.

Requires Python 3.10 or newer. No additional packages are needed.

## Project Structure

```text
TCR_Search/
├── main.py                  # Program entry point and search workflow
├── tcr_search/
│   ├── data_input.py        # Read and check input files
│   ├── inverted_index.py    # Build the index and generate candidates
│   ├── verification.py      # Verify sequence edit distances
│   ├── exhaustive_search.py # Search all eligible references
│   └── benchmark.py         # Compare the two search methods
├── data/
│   ├── real/
│   │   ├── reference_all.tsv
│   │   ├── reference_n1000.tsv
│   │   ├── reference_n5000.tsv
│   │   ├── reference_n10000.tsv
│   │   ├── queries.tsv
│   │   └── annotations.jsonl.gz
│   └── synthetic/
│       ├── length12_gqe_frequency/
│       │   ├── queries.tsv          # 61 general edit queries
│       │   ├── queries_gqe.tsv      # 10 common-seed probes
│       │   ├── reference_gqe_0pct.tsv
│       │   ├── reference_gqe_50pct.tsv
│       │   └── reference_gqe_90pct.tsv
│       └── sequence_length_comparison/
│           ├── length8/            # queries.tsv and reference.tsv
│           ├── length12/           # queries.tsv and reference.tsv
│           └── length18/           # queries.tsv and reference.tsv
├── results/                 # Generated output files
├── tests/                   # Input and search tests
└── README.md
```

## Usage

Run these commands from the `TCR_Search` directory. Relative paths are resolved from your current working directory.

`search`, `validate`, `benchmark`, and `memory` require `--reference` and `--queries`. `summarize` reads saved benchmark reports instead. Running `python main.py` without arguments displays help.

### 1. Search

Search the 1,000-sequence reference subset and include the full stored annotations:

```bash
python main.py search --reference data/real/reference_n1000.tsv --queries data/real/queries.tsv --annotations data/real/annotations.jsonl.gz --d 1 --k 3 --output results/real_n1000.jsonl
```

`--d` sets the inclusive maximum edit distance (default `1`), and `--k` sets the seed length (default `3`). Use `--method exhaustive` to search without the index. Replace the reference path with `data/real/reference_n5000.tsv`, `data/real/reference_n10000.tsv`, or `data/real/reference_all.tsv` to change the reference size; the full reference contains 111,414 unique sequences. All four reference files share the same annotation file.

`--annotations` accepts `.jsonl` or `.jsonl.gz`; gzip files are decompressed while reading, so manual extraction is unnecessary. If the option is omitted, search looks beside the reference TSV for `annotations.jsonl.gz`, then `annotations.jsonl`. If neither exists, results contain empty annotation lists. An explicitly supplied file must exist, and a loaded annotation file must contain records for every selected reference ID.

### 2. Benchmark

Each experiment compares indexed search with length-filtered exhaustive search using the same verifier. Matching IDs and exact distances are checked before timing. Both methods therefore use the same edit-distance calculation; these experiments measure indexing and candidate filtering, not DP optimisation.

The commands below are for PowerShell, run from `TCR_Search`. Each command creates its output directory, writes one CSV and one `.metadata.json` per run, and automatically generates `summary.csv`. Start with a missing or empty output directory. Choose a new `--output-dir` to repeat an experiment; existing results are not overwritten.

Use `--vary` to name the condition being changed. All other experimental conditions are kept fixed. `--d` and `--k` accept multiple values for their corresponding comparison. Seed positions are selected automatically. For length comparisons, query and reference files are paired in the order given, rather than searching every reference with every query file.

#### A. Real data: database size

Question: how do search time, speedup, index construction time and memory change as the database grows? Use the same 100 real queries with the nested reference subsets of 1,000, 5,000 and 10,000 sequences.

```powershell
python main.py benchmark `
  --vary size `
  --reference data/real/reference_n1000.tsv data/real/reference_n5000.tsv data/real/reference_n10000.tsv `
  --queries data/real/queries.tsv `
  --d 1 --k 3 --repeats 5 `
  --output-dir results/real_size
```

Outputs: `n1000.csv`, `n5000.csv`, `n10000.csv`, their metadata files, and `summary.csv`. Each row of the summary describes a different reference size. The supplied real queries are exact copies present in all three subsets; use the synthetic experiments below to cover edited queries.

#### B. Synthetic data: GQE frequency

Question: how does a common seed affect candidate counts and speedup? Keep reference length 12 and size 1,000; vary the proportion containing GQE from 0% to 50% to 90%.

Use `queries_gqe.tsv`, the 10 queries copied from the original data package's `synthetic/motif_probes/L12/queries.tsv`. With `k=3` and `d=2`, every probe selects GQE as one of its automatically chosen seeds. Distance and seed-length experiments use the 61 general edit queries in `sequence_length_comparison/length12/`.

```powershell
$data = "data/synthetic/length12_gqe_frequency"
python main.py benchmark `
  --vary gqe `
  --reference "$data/reference_gqe_0pct.tsv" "$data/reference_gqe_50pct.tsv" "$data/reference_gqe_90pct.tsv" `
  --queries "$data/queries_gqe.tsv" `
  --d 2 --k 3 --repeats 5 `
  --output-dir results/gqe_frequency
```

Outputs: `gqe_0pct.csv`, `gqe_50pct.csv`, `gqe_90pct.csv`, their metadata files, and `summary.csv`. Compare average candidate counts and search times; the proportion of references sharing a seed does not by itself determine the number of final matches.

#### C. Synthetic data: sequence length

Question: how does sequence length affect the benefit of indexing? Compare the length-8, length-12 and length-18 data groups, each with 1,000 reference sequences and 61 queries of the same edit-category composition. All reference databases have 0% GQE. Each database uses its own paired queries.

```powershell
$data = "data/synthetic/sequence_length_comparison"
python main.py benchmark `
  --vary length `
  --reference "$data/length8/reference.tsv" "$data/length12/reference.tsv" "$data/length18/reference.tsv" `
  --queries "$data/length8/queries.tsv" "$data/length12/queries.tsv" "$data/length18/queries.tsv" `
  --d 2 --k 3 --repeats 5 `
  --output-dir results/sequence_length
```

Outputs: `length8.csv`, `length12.csv`, `length18.csv`, their metadata files, and `summary.csv`. Compare speedup, average candidates and `fallback_count`. A query needs at least `k * (d + 1)` characters for indexed seed selection; with these parameters a length-8 query falls back to verifying all length-eligible references. The groups include insertion, deletion and short-query cases, so individual query lengths differ from the reference length. Different groups also contain different sequences; interpret this as a comparison under the same generation rules, rather than an isolated causal effect of length alone.

#### D. Distance threshold

Question: what happens when the permitted edit distance increases? Use the 0%-GQE length-12 reference and its 61 general edit queries from `sequence_length_comparison/length12/`, varying only `d`.

```powershell
$data = "data/synthetic/sequence_length_comparison/length12"
python main.py benchmark `
  --vary distance `
  --reference "$data/reference.tsv" --queries "$data/queries.tsv" `
  --d 0 1 2 3 --k 3 --repeats 5 `
  --output-dir results/distance_threshold
```

Outputs: `d0.csv` through `d3.csv`, their metadata files, and `summary.csv`. A larger `d` changes both the valid match set and how many seeds are required. Compare candidate counts, fallback counts and speedup while checking that both methods agree at each threshold. Do not confuse `d`, the search acceptance threshold, with `applied_edits`, which records how a query was generated.

#### E. Seed length

Vary seed length `k`, keeping the same length-12 data and `d`. Seed positions are selected automatically:

```powershell
$data = "data/synthetic/sequence_length_comparison/length12"
python main.py benchmark `
  --vary seed_length `
  --reference "$data/reference.tsv" --queries "$data/queries.tsv" `
  --d 2 --k 2 3 4 --repeats 5 `
  --output-dir results/seed_length
```

Outputs: `k2.csv`, `k3.csv`, `k4.csv`, their metadata files, and `summary.csv`. `k` is the length of each seed, not the number of seeds. Longer seeds may give smaller posting lists, but more queries may be too short to fit `d + 1` disjoint seeds.

#### Reading the reports

Each benchmark CSV contains an `exhaustive` row and an `indexed` row. Times are measured over five repeated query batches by default, with mean, median, standard error and samples recorded. A timing sample is the whole batch's elapsed time divided by its query count; it is not the median of individual-query latencies. Standard error is calculated across the repeated batch averages, not as an error bar for the median. Index construction time is measured separately; memory is measured in separate tracing runs after timing. Input loading, correctness checks, annotation attachment and file output are outside the timed search loops.

Each summary has one row per run, sorted by the changed condition:

| Summary column | Meaning |
|---|---|
| `comparison`, `reference` | Experiment type and reference name |
| `gqe_percentage` | Included only with `--vary gqe`; percentage label from the reference filename |
| `reference_count`, `query_count` | Input sizes |
| `average_reference_length`, `average_query_length` | Mean sequence lengths; a length-12 reference group can contain queries of other lengths |
| `d`, `k` | Maximum edit distance and seed length |
| `fallback_count` | Number of queries using the short-query fallback |
| `exhaustive_median_ms_per_query`, `indexed_median_ms_per_query` | Median search time per query for each method |
| `search_speedup` | Exhaustive median divided by indexed median; excludes index construction |
| `average_candidates_per_query` | Total indexed candidates divided by query count |
| `index_build_median_ms` | Median time to build one index |
| `exhaustive_incremental_peak_python_mib`, `indexed_incremental_peak_python_mib` | Incremental peak Python allocations in MiB; indexed includes index construction, not total process memory or RSS |
| `memory_ratio` | Indexed incremental peak memory divided by exhaustive incremental peak memory: greater than 1 means more memory, less than 1 means less. Blank if exhaustive peak is zero. This is an incremental-memory ratio, not a total-process-memory ratio. |

The summary checks the conditions that should remain fixed for the selected `--vary`. Query paths and IDs are compared except in length experiments, where different query sets are intentional. Use unchanged input files and the same execution environment throughout a comparison. Older reports may leave new length/fallback fields blank. For means and standard errors, consult the individual benchmark CSVs.

For one benchmark without a batch summary, the previous single-reference command with `--output` still works. `--output` and `--output-dir` cannot be combined; a custom `--metadata` path applies only to single-file output. To rebuild a summary from existing reports without rerunning benchmarks, use `python main.py summarize --vary size --input-dir results/real_size --output results/real_size/summary_copy.csv`, changing `--vary` and paths for the desired experiment.

## Input and Output

Input files use tab-separated columns:

- Reference file: `reference_id` and `sequence`.
- Query file: `query_id` and `sequence`.

IDs must be unique within each file. Sequences must contain uppercase standard amino-acid letters. Your own query file only needs `query_id` and `sequence`; extra query columns are retained as text in `query_metadata`.

The supplied real query file contains 100 `exact_copy` queries: each query was copied unchanged from its `source_id`, with `applied_edits = 0` and `source_distance = 0`. All 100 source sequences occur in each supplied reference subset. These columns describe query preparation and do not control search. With `d=0`, search returns exact matches; with a larger threshold, it also returns any other sequences within that distance. `source_distance` describes the query's distance to its source, while each result's `distance` is computed against that particular matched reference.

Search produces one JSON Lines record per query, including queries with no matches:

- `query_id` and `query_sequence`: the input query.
- `query_metadata`: all additional query columns.
- `matches`: all reference sequences within the threshold, or `[]` if none match. Each match contains `reference_id`, `sequence`, the exact `distance`, and `annotations`.

The annotation file stores each reference ID with a `records` list. Search returns this entire list as `annotations`, preserving each record's `source_line` and `fields`. A sequence can have several original records with different epitopes, MHC values, or references; these records remain separate. Completeness here means all fields stored in the supplied annotation file, which was prepared from the VDJdb slim dataset.

By default, search results appear in the terminal. Use `--output` to save them to a file. Existing output files are not overwritten; choose a new path when needed.

## Tests

```bash
python -B -m unittest discover -s tests
```
