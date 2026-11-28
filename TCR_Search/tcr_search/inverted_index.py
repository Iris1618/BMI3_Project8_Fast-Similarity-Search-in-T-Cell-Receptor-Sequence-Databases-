"""Build a k-mer index and find candidates for edit-distance verification."""
from collections import defaultdict
from types import MappingProxyType
from .data_input import STRATEGIES, check_id, check_sequence, require_integer

def select_seed_positions(m, k, d, strategy="spread"):
    """Return d+1 disjoint seed starts, or an empty tuple for fallback."""
    require_integer(m, "m", 0)
    require_integer(k, "k", 1)
    require_integer(d, "d", 0)
    if strategy not in STRATEGIES:
        raise ValueError("strategy must be spread or prefix")
    if m < k * (d + 1):
        return ()
    if d == 0 or strategy == "prefix":
        return tuple(i * k for i in range(d + 1))
    return tuple(i * (m - k) // d for i in range(d + 1))  # End at m - k.

class InvertedIndex:
    def __init__(self, references, k=3):
        require_integer(k, "k", 1)
        self._k = k
        self.references = MappingProxyType(dict(references))  # Keep a read-only copy.
        postings = defaultdict(set)
        by_length = defaultdict(set)
        for rid, sequence in self.references.items():
            check_id(rid)
            check_sequence(sequence, allow_empty=True)
            by_length[len(sequence)].add(rid)
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
        groups = [ids for n, ids in self.by_length.items() if abs(n - m) <= d]
        hits = set()
        if not starts:  # Too short for d+1 seeds. Use all length-eligible IDs.
            for ids in groups:
                hits.update(ids)
            return CandidateResult(hits, "fallback", len(hits))

        # Any reference within d edits must contain at least one of these seeds.
        seeds = [query[start:start + self.k] for start in starts]
        sizes = []
        for seed in seeds:
            ids = self.postings.get(seed, ())
            hits.update(ids)
            sizes.append(len(ids))
        union_count = len(hits)  # Unique seed hits before length filtering.
        eligible = sum(len(ids) for ids in groups)
        if eligible < len(self.references):
            hits = {rid for rid in hits if abs(len(self.references[rid]) - m) <= d}
        return CandidateResult(hits, "indexed", eligible, seed_positions=starts,
                               seeds=seeds, posting_sizes=sizes, union_count=union_count)


class CandidateResult:
    """Candidate IDs and the evidence used to select them."""
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
    # Replace each set in turn to limit peak memory during construction.
    for key in groups:
        groups[key] = frozenset(groups[key])
    return MappingProxyType(dict(groups))
