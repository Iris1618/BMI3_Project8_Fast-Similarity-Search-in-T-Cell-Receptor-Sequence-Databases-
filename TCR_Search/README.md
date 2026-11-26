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
├── data/examples/
│   ├── reference.tsv
│   └── queries.tsv
├── results/                 # Generated output files
├── tests/                   # Empty for now
└── README.md
```

## Usage

Run these commands from the project directory:

### 1. Search

Run indexed search on the examples using the default settings (`d=1`, `k=3`):

```bash
python main.py
```

Set the parameters and save the results, or switch to exhaustive search:

```bash
# Indexed search
python main.py search --d 1 --k 2 --output results/indexed.jsonl

# Exhaustive search
python main.py search --method exhaustive --d 1 --output results/exhaustive.jsonl
```

`--d` sets the maximum edit distance, and `--k` sets the seed length.

### 2. Benchmark

Check that both search methods return the same IDs and distances, then compare their execution time and memory usage:

```bash
python main.py benchmark --d 1 --k 2
```

Both methods use the same verifier. Timing is repeated five times by default, and index construction time is reported separately. Memory is measured in a separate run after timing.

The reports are saved as `results/benchmark.csv` and `results/benchmark.metadata.json`. Memory figures cover Python allocations after input loading, including index construction for indexed search.

## Input and Output

Input files use tab-separated columns:

- Reference file: `reference_id` and `sequence`.
- Query file: `query_id` and `sequence`.

IDs must be unique within each file. Sequences must contain uppercase standard amino-acid letters. Use `--reference` and `--queries` to provide your own files.

The examples contain 14 references and 8 queries. The default search returns 16 matches.

Search results contain reference IDs, sequences, and exact distances. By default, search results appear in the terminal. Exported results and benchmark reports are saved in `results/`. Existing files are not overwritten; choose a new `--output` path when needed.
