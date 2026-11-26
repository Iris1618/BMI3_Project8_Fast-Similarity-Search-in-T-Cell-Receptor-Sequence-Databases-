# Input preparation

This folder prepares inputs for the group's TCR search program.
Run the commands below from the **repository root**, using Python 3.10 or newer.
No third-party packages are needed.

## Read the code in this order

1. common.py: check amino-acid letters, assign IDs and save files.
2. prepare_real.py: read VDJdb, clean sequences and create nested subsets.
3. prepare_synthetic.py: generate references, edit templates and create GQE queries.

The scripts use ordinary loops, lists, dictionaries and sets. A set records values
already seen so duplicate checks stay fast. There are no classes or comprehensions.

| Python standard-library module | Why it is used |
| --- | --- |
| csv | Write tab-separated tables. |
| json | Save annotations, edit histories and settings. |
| random | Generate repeatable data using fixed seeds. |
| argparse | Read settings supplied on the command line. |
| pathlib | Work with file and directory paths. |
| sys | Record the Python version. |
| gzip | Save the existing compressed annotation format. |
| hashlib | Keep existing sequence IDs and record the source-file checksum. |

The synthetic script borrows the group's basic DP function from
draft/Verification/verification_basic.py to check each query against its source.
It does not calculate the complete retrieval truth.

## Real data

Download and extract the [VDJdb 2026-06-03 release](https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03/vdjdb-2026-06-03.zip).
Use vdjdb.slim.txt as the input. See the [VDJdb paper](https://doi.org/10.1093/nar/gkx760).

~~~bash
python3 -m draft.Input_Source.prepare_real \
  --input /path/to/vdjdb-2026-06-03/vdjdb.slim.txt \
  --out data/real \
  --sizes 1000 5000 10000 --queries 100 --seed 20261111
~~~

The defaults select human TRB records. Missing or invalid CDR3 sequences are
excluded; outer whitespace is removed and lowercase letters become uppercase.
Identical sequences are grouped while every original annotation row is retained.
Lengths outside 8-30 aa are flagged and retained.

For the pinned release, expect 120,575 valid human TRB records and 111,414 unique
sequences. The script shuffles once and takes prefixes to make nested reference
subsets. Queries are sampled from the smallest subset, so every subset contains
their exact matches. Species and chain can be changed with --species and --gene.

| Output | Contents |
| --- | --- |
| reference_all.tsv | All unique reference IDs and sequences. |
| reference_n*.tsv | Nested reference subsets of the requested sizes. |
| queries.tsv | Exact-copy queries with their source IDs. |
| annotations.jsonl.gz | All original annotation records grouped by reference ID. |
| excluded_target_rows.tsv | Missing or invalid target CDR3 records and their reasons. |
| length_distribution.tsv | Unique sequence counts by length. |
| length_review_flags.tsv | Sequences outside the review length range. |
| cleaning_report.json | Counts, settings, source checksum and Python version. |

The reference sizes and query count are parameters. We provide the current
settings as a starting dataset. If the group needs larger datasets, change the
parameters and regenerate the inputs. The sizes 1,000/5,000/10,000 and 100 queries
are practical initial settings; record the parameters used in each experiment.

## Synthetic data

~~~bash
python3 -m draft.Input_Source.prepare_synthetic \
  --out data/synthetic \
  --size 1000 --lengths 8 12 18 --motif-percentages 0 50 90 \
  --anchors 5 --max-edits 3 --probe-count 10
~~~

The defaults create nine reference libraries: three lengths multiplied by three
GQE proportions, with 1,000 unique references in each. References start with C
and end with F. GQE is placed centrally in motif-positive sequences. Proportions
are rounded to whole sequence counts and must leave enough background templates.
The sequence length must provide enough distinct motif-positive combinations.

Each library has 61 ordinary queries: 5 exact copies, 15 substitution queries,
15 insertion queries, 15 deletion queries, 6 short queries and 5 candidate negatives.
The five background templates are shared across proportions of the same length.
The default synthetic seed is 20262111.

- --anchors chooses the number of background templates.
- --max-edits sets the largest requested source-query distance during generation.
- --probe-count sets the number of common GQE queries per length; use 0 to omit them.

Each common GQE query changes one residue outside GQE in a shared positive
template. It is absent as an exact sequence from every reference library.
Use the same motif_probes/L12/queries.tsv file against all three length-12 libraries;
likewise use L8 and L18 against the corresponding lengths. The default probe seed
is 20261113. For length 12, k=3, d=2 and spread placement select GQE as one seed.

| Output | Contents |
| --- | --- |
| L*_motif*/reference.tsv | Reference IDs and sequences for one condition. |
| L*_motif*/queries.tsv | Ordinary queries and source information. |
| L*_motif*/edit_history.jsonl | Operations, zero-based positions, old and new residues. |
| L*_motif*/generation_report.json | Actual motif count and query categories. |
| motif_probes/L*/queries.tsv | Common GQE queries for one length. |
| generation_report.json | Settings for the complete synthetic generation run. |

The edits are amino-acid stress tests at interior positions. V/D/J junction
positions are not inferred. source_distance describes the distance to the source
template, not the distance to every reference. candidate_no_match queries still
need exhaustive validation at the chosen search threshold.

## Handoff to search and benchmark

Reference TSV columns: reference_id, sequence.
Query TSV columns: query_id, sequence, category, source_id, source_sequence,
applied_edits, source_distance. Empty metadata means it does not apply.
The search reader accepts these extra query columns.

~~~bash
python3 TCR_Search/main.py search \
  --reference data/synthetic/L12_motif050/reference.tsv \
  --queries data/synthetic/motif_probes/L12/queries.tsv \
  --d 2 --k 3 --strategy spread

python3 TCR_Search/main.py validate \
  --reference data/synthetic/L12_motif050/reference.tsv \
  --queries data/synthetic/motif_probes/L12/queries.tsv \
  --d 2 --k 3 --output results/check_L12_motif050.csv
~~~

Repeat with motif000 and motif090 while keeping the query file and search settings
the same. Benchmark candidate counts, query time and index costs.

Use a fresh output directory for each generation run. Files are written directly
to that directory. If generation fails after writing starts, the directory may
contain incomplete files; use a new directory for the next attempt.
Fixed seeds and the same source file reproduce the sequence TSVs and IDs.
Check supplied truth files against the exact references and queries used.

The current search application leaves annotations empty. The saved real-data
annotations can be joined to matches by reference_id after retrieval.

## English pseudocode

~~~text
READ VDJdb records
SELECT the requested species and chain
CLEAN and CHECK each CDR3 sequence
GROUP identical sequences and retain all annotations
SHUFFLE the unique sequences with a fixed seed
WRITE nested reference subsets
SAMPLE queries from the smallest subset
SAVE counts and preparation settings
~~~

~~~text
FOR each sequence length
    GENERATE pools with and without GQE
    FOR each requested GQE proportion
        MIX the pools and SAVE the reference library
        COPY templates and APPLY substitutions, insertions and deletions
        CHECK each query against its source using the group's verifier
        ADD short queries and candidate negatives
        SAVE queries, edit histories and settings
    CHANGE residues outside GQE in shared templates
    REJECT queries identical to any reference
    SAVE the common GQE queries
~~~

## Contribution and assistance

This directory covers the input/source preparation contribution. AI assistance
was used to draft and simplify these scripts and this README. Record that
assistance in the coursework disclosure and follow the module's AI rules.
