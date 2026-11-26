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
│   └── synthetic/             # Empty for now
├── results/                 # Generated output files
├── tests/                   # Input and search tests
└── README.md
```

## Usage

Run these commands from the `TCR_Search` directory. Relative paths are resolved from your current working directory.

Every command (`search`, `validate`, `benchmark`, or `memory`) requires `--reference` and `--queries`. Running `python main.py` without arguments displays help.

### 1. Search

Search the 1,000-sequence reference subset and include the full stored annotations:

```bash
python main.py search --reference data/real/reference_n1000.tsv --queries data/real/queries.tsv --annotations data/real/annotations.jsonl.gz --d 1 --k 3 --output results/real_n1000.jsonl
```

`--d` sets the inclusive maximum edit distance (default `1`), and `--k` sets the seed length (default `3`). Use `--method exhaustive` to search without the index. Replace the reference path with `data/real/reference_n5000.tsv`, `data/real/reference_n10000.tsv`, or `data/real/reference_all.tsv` to change the reference size; the full reference contains 111,414 unique sequences. All four reference files share the same annotation file.

`--annotations` accepts `.jsonl` or `.jsonl.gz`; gzip files are decompressed while reading, so manual extraction is unnecessary. If the option is omitted, search looks beside the reference TSV for `annotations.jsonl.gz`, then `annotations.jsonl`. If neither exists, results contain empty annotation lists. An explicitly supplied file must exist, and a loaded annotation file must contain records for every selected reference ID.

### 2. Benchmark

Check that both search methods return the same IDs and distances, then compare their execution time and memory usage:

```bash
python main.py benchmark --reference data/real/reference_n1000.tsv --queries data/real/queries.tsv --d 1 --k 3 --output results/benchmark_real_n1000.csv
```

Both methods use the same verifier. Timing is repeated five times by default, and index construction time is reported separately. Memory is measured in a separate run after timing.

This command saves `results/benchmark_real_n1000.csv` and `results/benchmark_real_n1000.metadata.json`. Without `--output`, the default names are `results/benchmark.csv` and `results/benchmark.metadata.json`. Search timing excludes data loading, validation, annotation attachment, and file output. Memory figures cover incremental Python allocations after input loading, including index construction for indexed search; they are not total process memory. Benchmarking does not load annotations.

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
