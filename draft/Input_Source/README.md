# Input preparation

Run the commands below from the repository root. Python 3.9 or later and the
Python standard library are sufficient. The synthetic script imports the group's
`verification_basic` function to check each edited query against its source template.

## Real data

Download the [VDJdb 2026-06-03 release](https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03/vdjdb-2026-06-03.zip)
and extract it. Use `vdjdb.slim.txt` as the input. The release is described in the
[VDJdb paper](https://doi.org/10.1093/nar/gkx760).

```bash
python3 -m Input_Source.prepare_real \
  --input /path/to/vdjdb-2026-06-03/vdjdb.slim.txt \
  --out data/real \
  --sizes 1000 5000 10000 --queries 100 --seed 20261111
```

The script selects human TRB records, removes missing or invalid CDR3 sequences,
and deduplicates by CDR3 amino acid sequence. It trims outer whitespace and uses
uppercase letters. All original annotation rows remain linked to their reference ID.
Lengths outside 8-30 aa are flagged for review and retained.

The expected counts for this release are 120,575 valid human TRB records and
111,414 unique sequences. Nested subsets share one shuffled sequence order.
The queries are copied from the smallest subset, so every subset contains their
exact matches. `--species` and `--gene` can be set for another source selection.

The default command prints

```text
Kept 120575 records and 111414 unique sequences
```

The reference sizes and query count are configurable. We provide the current
settings as a starting dataset. If the group needs larger datasets, change the
parameters and regenerate the inputs.

The 1,000/5,000/10,000 sizes and 100 queries are practical initial benchmark
settings. Compare them using the same queries and record the chosen parameters.

| Output | Contents |
| --- | --- |
| `reference_all.tsv` | All unique reference IDs and sequences. |
| `reference_n*.tsv` | Nested reference subsets of the requested sizes. |
| `queries.tsv` | Exact-copy queries with their source IDs. |
| `annotations.jsonl.gz` | Each reference ID and its full original annotation records. |
| `excluded_target_rows.tsv` | Missing or invalid CDR3 rows among the selected species and chain. |
| `length_distribution.tsv` | Unique sequence counts by length. |
| `length_review_flags.tsv` | Sequences outside the review length range. |
| `cleaning_report.json` | Source filename and SHA-256, counts, parameters and Python version. |

## Synthetic data

```bash
python3 -m Input_Source.prepare_synthetic \
  --out data/synthetic \
  --size 1000 --lengths 8 12 18 --motif-percentages 0 50 90 \
  --max-edits 3 --anchors 5 --probe-count 10
```

Each library has unique sequences starting with C and ending with F.
`L12_motif050` contains length-12 references with 50% containing GQE. GQE is
placed centrally in the motif-positive pool. Percentages are rounded to a whole
number of sequences. The requested proportions must leave at least `--anchors`
background sequences in each library, and the motif-positive pool must fit the
available sequence combinations. Increase the length when generating larger pools.

With the defaults, each of nine libraries has 1,000 references and 61 queries.
Five common background templates produce 5 exact copies, 15 substitutions,
15 insertions and 15 deletions. Another 6 queries have lengths 1-3, and 5
random W/Y queries are labelled `candidate_no_match` for the benchmark member
to check against the whole reference library. Edit positions use zero-based
indices in the sequence at that step. `source_distance` applies to the source
template. Use that reference-query pair when checking a distance boundary.

These are amino acid stress tests with edits at interior positions. They test
search behaviour under controlled changes. V/D/J junction positions are not
inferred during generation.

The `motif_probes/L8`, `L12` and `L18` directories each contain 10 common GQE
queries. A probe changes one residue outside GQE in a template shared by all
positive-proportion libraries. Each probe is absent as an exact sequence from
every library. Use the same probe file against all three proportions of that length.
Set `--probe-count 0` to omit probes. The default synthetic seed is 20262111,
and the probe seed is 20261113. Both seeds are configurable.

| Output in each `L*_motif*` directory | Contents |
| --- | --- |
| `reference.tsv` | Reference IDs and sequences for this condition. |
| `queries.tsv` | Generated queries and source information. |
| `edit_history.jsonl` | The operations used to construct each template-derived query. |
| `generation_report.json` | Reference length, actual motif count and query categories. |

`generation_report.json` in the synthetic output root records the complete
generation settings. `motif_probes/L*/queries.tsv` supplies the common queries.
The short, hand-checkable files in `Input_Source/examples/` can be used to check
fallback behaviour with true matches.

For example, generate a smaller practice set with

```bash
python3 -m Input_Source.prepare_synthetic \
  --out data/practice --size 100 --lengths 12 \
  --motif-percentages 0 50 90 --probe-count 5
```

Use a fresh output directory for each run. Files are published after preparation
finishes. Fixed seeds and input files give reproducible sequence TSVs. These
scripts prepare inputs for the group's exhaustive search and benchmark stage.
Existing truth files must be checked against the exact query and reference files
used in each run.

## Data fields and handoff

Reference TSV files contain `reference_id` and `sequence`. Query TSV files contain
`query_id` and `sequence`, followed by `category`, `source_id`, `source_sequence`,
`applied_edits` and `source_distance`. Empty metadata means it does not apply.
IDs are stable within the files and reference IDs are derived from the sequence.

The group's index can read the generated files directly. For example

```bash
python3 Inverted_Index/inverted_index.py \
  --reference data/synthetic/L12_motif050/reference.tsv \
  --queries data/synthetic/motif_probes/L12/queries.tsv \
  --k 3 --d 2 --strategy spread --output candidates.jsonl
```

For the length-12 GQE experiment, `spread`, k=3 and d=2 select starts 0, 4 and 9.
The seed at start 4 is GQE. Record the actual seeds and candidate counts when
comparing the three motif proportions.

## Pseudocode

```text
READ the VDJdb column names
FOR each record
    SELECT the requested species and chain
    CLEAN and CHECK the CDR3 amino acid sequence
    GROUP identical sequences and preserve every annotation record
SORT the unique sequences and shuffle using a fixed seed
WRITE nested reference subsets
SAMPLE queries from the smallest subset
WRITE the preparation settings and counts
```

```text
FOR each reference length
    GENERATE unique background and GQE-containing sequence pools
    FOR each requested GQE percentage
        SELECT sequences from the two pools and write the reference library
        COPY common templates and apply substitutions, insertions and deletions
        CHECK the source-query distances using the group's verifier
        ADD short queries and candidate negative queries
        SAVE queries, edit histories and generation settings
    GENERATE one common GQE probe set from shared positive templates
    CHECK that each probe is absent as an exact sequence from every library
    WRITE the common probes
```

## Contribution and assistance

This directory covers the input/source preparation contribution. AI assistance
was used to draft the input preparation scripts and this README. The group should
record this assistance in its coursework disclosure and follow the module's AI rules.
