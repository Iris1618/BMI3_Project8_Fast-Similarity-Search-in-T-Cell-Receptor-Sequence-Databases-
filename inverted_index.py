"""Build a k-mer index of TCR sequences and find candidate reference IDs for each query."""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from types import MappingProxyType

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
STRATEGIES = ("spread", "prefix")


# 1. Build the index, select seeds, and find candidates.
class InvertedIndex:
    def __init__(self, references, k=3):
        require_integer(k, "k", 1)
        self._k = k
        self.references = MappingProxyType(dict(references))
        postings = defaultdict(set)
        by_length = defaultdict(set)

        for rid, sequence in self.references.items():
            check_id(rid)
            check_sequence(sequence, allow_empty=True)
            by_length[len(sequence)].add(rid)
            # With k=2, CASS adds this ID to the lists for CA, AS, and SS.
            for start in range(len(sequence) - k + 1):
                kmer = sequence[start:start + k]
                postings[kmer].add(rid)

        self.postings = _freeze_groups(postings)
        self.by_length = _freeze_groups(by_length)

    @property
    def k(self):
        return self._k

    def generate_candidates(self, query, d, strategy="spread"):
        check_sequence(query, allow_empty=True)
        m = len(query)
        starts = select_seed_positions(m, self.k, d, strategy)
        length_groups = [ids for length, ids in self.by_length.items() if abs(length - m) <= d]

        # Short query: keep every reference within d residues of the query length.
        if not starts:
            candidates = frozenset().union(*length_groups)
            return CandidateResult(candidates, "fallback", len(candidates))

        # At most d edits can disrupt d of the d+1 seeds. A hit on any seed is enough here.
        seeds = tuple(query[start:start + self.k] for start in starts)
        postings = [self.postings.get(seed, ()) for seed in seeds]
        hits = frozenset().union(*postings)
        eligible_count = sum(len(ids) for ids in length_groups)
        candidates = hits  # Reuse all hits when every reference passes the length filter.
        if eligible_count < len(self.references):
            candidates = frozenset(rid for rid in hits if abs(len(self.references[rid]) - m) <= d)
        return CandidateResult(
            candidates, "indexed", eligible_count, seed_positions=starts, seeds=seeds,
            posting_sizes=[len(ids) for ids in postings], union_count=len(hits),
        )


def select_seed_positions(m, k, d, strategy="spread"):
    """Return zero-based starts of d+1 disjoint seeds, or an empty tuple for fallback."""
    require_integer(m, "m", 0)
    require_integer(k, "k", 1)
    require_integer(d, "d", 0)
    if strategy not in STRATEGIES:
        raise ValueError("strategy must be spread or prefix")
    if m < k * (d + 1):
        return ()
    if d == 0:
        return (0,)
    if strategy == "prefix":
        return tuple(i * k for i in range(d + 1))
    # Spread the starts from 0 to m-k, rounding each position down.
    return tuple(i * (m - k) // d for i in range(d + 1))


class CandidateResult:
    """Candidate IDs for the distance verifier, along with the seed lookup details."""

    def __init__(self, candidate_ids, mode, length_eligible_count,
                 seed_positions=(), seeds=(), posting_sizes=(), union_count=None):
        self.candidate_ids = frozenset(candidate_ids)
        self.mode = mode
        self.length_eligible_count = length_eligible_count
        self.seed_positions = tuple(seed_positions)
        self.seeds = tuple(seeds)
        self.posting_sizes = tuple(posting_sizes)
        self.union_count = union_count
        self.fallback_reason = "query_length_lt_k_times_d_plus_1" if mode == "fallback" else None

    def to_dict(self):
        return dict(candidate_ids=sorted(self.candidate_ids), candidate_count=len(self.candidate_ids),
                    mode=self.mode, fallback_reason=self.fallback_reason,
                    seed_positions=list(self.seed_positions), seeds=list(self.seeds),
                    posting_sizes=list(self.posting_sizes), union_count=self.union_count,
                    length_eligible_count=self.length_eligible_count)


def _freeze_groups(groups):
    # Freeze each group in place so the old sets can be released one at a time.
    for key, ids in groups.items():
        groups[key] = frozenset(ids)
    return MappingProxyType(dict(groups))


def require_integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def check_id(identifier):
    if not isinstance(identifier, str) or not identifier or identifier.strip() != identifier:
        raise ValueError("ID must be a nonempty string without outer whitespace")


def check_sequence(sequence, allow_empty=False):
    if not isinstance(sequence, str) or (not sequence and not allow_empty):
        raise ValueError("expected an amino acid sequence")
    for letter in sequence:
        if letter not in AMINO_ACIDS:
            raise ValueError(f"invalid amino acid: {letter!r}; use uppercase letters")


# 2. Read input and write results.
def read_sequences(path, id_column):
    """Read a UTF-8 TSV file and map each record ID to its sequence."""
    records = {}
    with open(path, encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file, delimiter="\t", strict=True)
        try:
            header = next(reader, [])
            if id_column not in header or "sequence" not in header:
                raise ValueError(f"requires {id_column} and sequence columns")
            if len(header) != len(set(header)):
                raise ValueError("duplicate column names")
            id_col, seq_col = header.index(id_column), header.index("sequence")
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError("malformed TSV row")
                rid, sequence = row[id_col], row[seq_col]
                check_id(rid)
                if rid in records:
                    raise ValueError(f"duplicate ID {rid}")
                check_sequence(sequence)
                records[rid] = sequence
        except csv.Error as error:
            raise ValueError(f"{path}:{reader.line_num}: malformed TSV: {error}") from error
        except ValueError as error:
            raise ValueError(f"{path}:{reader.line_num}: {error}") from error
    return records


def run_search(args, references, queries, output):
    index = InvertedIndex(references, args.k)
    counts = {"indexed": 0, "fallback": 0}
    for query_id, sequence in queries.items():
        result = index.generate_candidates(sequence, args.d, args.strategy)
        counts[result.mode] += 1
        row = result.to_dict()
        row.update(query_id=query_id, k=index.k, d=args.d, strategy=args.strategy)
        print(json.dumps(row, ensure_ascii=False), file=output)
    summary = dict(k=index.k, reference_count=len(references), query_count=len(queries),
                   indexed_queries=counts["indexed"], fallback_queries=counts["fallback"])
    print(json.dumps(summary), file=sys.stderr)


# 3. Command-line interface.
def main(argv=None):
    parser = argparse.ArgumentParser(description="Find candidate IDs before DP verification.")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--k", type=int, default=3, help="seed length (default: 3)")
    parser.add_argument("--d", type=int, required=True, help="maximum edit distance, inclusive")
    parser.add_argument("--strategy", choices=STRATEGIES, default="spread")
    parser.add_argument("--output", type=Path, help="candidate JSONL file (default: stdout)")
    args = parser.parse_args(argv)
    try:
        require_integer(args.k, "k", 1)
        require_integer(args.d, "d", 0)
        refs = read_sequences(args.reference, "reference_id")
        queries = read_sequences(args.queries, "query_id")
        if args.output is None:
            run_search(args, refs, queries, sys.stdout)
        else:
            with open(args.output, "a", encoding="utf-8", newline="\n") as output:
                # Check the open file before truncating it to protect both input files.
                target = os.fstat(output.fileno())
                if any(os.path.samestat(target, p.stat()) for p in (args.reference, args.queries)):
                    raise ValueError("output must not overwrite an input file")
                if output.seekable():
                    output.truncate(0)
                run_search(args, refs, queries, output)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")
    except KeyboardInterrupt:
        print("\nCancelled.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
