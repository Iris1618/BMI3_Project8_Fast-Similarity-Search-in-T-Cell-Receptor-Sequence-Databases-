"""Build a k-mer index and retrieve candidates for full-sequence verification."""

from collections import defaultdict
from types import MappingProxyType

from .data_input import AMINO_ACIDS, STRATEGIES, check_id, check_sequence, require_integer


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
