# Prepared Project 8 data

[Download Project8_Input_Data.zip](Project8_Input_Data.zip) (about 9.4 MiB).
This is the complete prepared benchmark data package, alongside the small
examples and generation scripts already in the repository.

## Unpack and run

From the repository root, extract the archive with Python:

```bash
python3 -m zipfile -e data/Project8_Input_Data.zip .
```

The archive already contains a data/ directory. After extraction, use:

```bash
python3 TCR_Search/main.py search \
  --reference data/real/reference_n1000.tsv \
  --queries data/real/queries.tsv --d 1 --k 3
```

For the common GQE queries:

```bash
python3 TCR_Search/main.py validate \
  --reference data/synthetic/L12_motif050/reference.tsv \
  --queries data/synthetic/motif_probes/L12/queries.tsv \
  --d 2 --k 3 --strategy spread --output results/check_L12_motif050.csv
```

Repeat with motif000 and motif090 using the same query file and settings.
Use a new output filename for each run.

## What is included

| Location after extraction | Contents |
| --- | --- |
| data/real/ | 111,414 unique human TRB sequences; nested 1,000/5,000/10,000-reference subsets; 100 common queries; annotations and cleaning reports. |
| data/synthetic/L*_motif*/ | Nine libraries: lengths 8/12/18 and GQE proportions 0/50/90%, with 1,000 references and 61 ordinary queries each. |
| data/synthetic/motif_probes/L*/ | Thirty common GQE queries, ten per length, shared across the three proportions. |
| data/synthetic/edge_cases/ | Small hand-checkable boundary examples. |
| truth_matches.tsv and truth_summary.tsv | Exhaustive expected results for d=0,1,2,3, beside the corresponding inputs. |
| data/STEP2_IO.txt | Input/output fields, no-match output and a real annotated example. |
| data/manifest.json and data/checksums.json | Provenance, generation settings and per-file SHA-256 checksums. |
| data/sources/VDJdb-LICENSE.txt | The license distributed with the VDJdb source release. |

## Checking truth

Select the rows for the exact reference_set, threshold and query_id. Compare
the full set of (reference_id, distance) pairs against the search output.
A truth_summary row with n_matches=0 records a query with no matches.

For example, the length-12 50% library uses reference_set=L12_motif050_n1000.
The real 1,000-reference subset uses reference_set=real_n1000.
Keep every query in the output, including those with an empty matches list.

The main.py validate command calculates its own exhaustive comparison; it does
not automatically read the bundled truth files. The benchmark member can use
the bundled truth as an additional check.

## Relation to the generation scripts

The package's reference/query IDs and sequences match the fixed-seed defaults
of draft/Input_Source/. Its negative queries have already been checked and are
labelled verified_no_match. Freshly generated queries start as candidate_no_match
until the benchmark member validates them.

The frozen package also has more detailed reports and an extra flag column in
length_review_flags.tsv. Annotation JSON formatting differs, but all 111,414
reference annotation objects agree with the current generator output.
Keep the frozen reports and checksums with this package.

## Source and integrity

Real data source: [VDJdb 2026-06-03 release](https://github.com/antigenomics/vdjdb-db/releases/tag/2026-06-03).
The archive contains prepared human TRB inputs; the original full VDJdb release
can be downloaded separately from that source. Retain its license and attribution.
Synthetic data are controlled amino-acid stress tests, not observed receptors.

The package contains 96 files, with 95 file checksums recorded inside it.
Archive SHA-256:

```text
f434be751925d3b20012edd3216634c0a71217a068c8705a7e36ff898597cf7b
```
