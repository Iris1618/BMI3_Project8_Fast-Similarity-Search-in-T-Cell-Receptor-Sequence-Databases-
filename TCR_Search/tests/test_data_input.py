"""Check lossless metadata loading and reference/annotation consistency."""

import gzip
import json
from pathlib import Path
import tempfile
import unittest

from tcr_search.data_input import read_annotations, read_records, read_sequences


class DataInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_annotations(self, rows, compressed=True):
        path = self.root / ("annotations.jsonl.gz" if compressed else "annotations.jsonl")
        opener = gzip.open if compressed else open
        with opener(path, "wt", encoding="utf-8") as file:
            for row in rows:
                file.write(json.dumps(row, ensure_ascii=False) + "\n")
        return path

    def annotation(self, rid="R1", sequence="CASS", source_line=2):
        return {"reference_id": rid, "records": [{
            "source_line": source_line,
            "fields": {"cdr3": sequence, "antigen.epitope": "PEPTIDE", "empty": ""},
        }]}

    def test_tsv_preserves_extra_columns_and_sequence_api(self):
        path = self.root / "queries.tsv"
        path.write_text(
            "query_id\tsequence\tcategory\tsource_id\tnotes\n"
            "Q1\tCASS\texact_copy\tR1\t原始记录\n"
            "Q2\tCAS\tdeletion\tR1\t\n", encoding="utf-8-sig",
        )
        self.assertEqual(read_records(path, "query_id"), {
            "Q1": {"query_id": "Q1", "sequence": "CASS", "category": "exact_copy",
                   "source_id": "R1", "notes": "原始记录"},
            "Q2": {"query_id": "Q2", "sequence": "CAS", "category": "deletion",
                   "source_id": "R1", "notes": ""},
        })
        self.assertEqual(read_sequences(path, "query_id"), {"Q1": "CASS", "Q2": "CAS"})

    def test_tsv_rejects_duplicate_ids_with_line_context(self):
        path = self.root / "queries.tsv"
        path.write_text("query_id\tsequence\nQ1\tCASS\nQ1\tCAS\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, r"queries.tsv:3: duplicate ID Q1"):
            read_records(path, "query_id")

    def test_gzip_and_plain_keep_all_records_for_subset(self):
        selected = self.annotation()
        second = self.annotation(source_line=9)["records"][0]
        second["fields"].update({"antigen.epitope": "OTHER", "mhc.a": "HLA-A*02:01"})
        selected["records"].append(second)
        rows = [self.annotation("UNSELECTED", "CA"), selected]
        for compressed in (True, False):
            with self.subTest(compressed=compressed):
                path = self.write_annotations(rows, compressed)
                self.assertEqual(read_annotations(path, {"R1": "CASS"}),
                                 {"R1": selected["records"]})
                self.assertEqual(sorted(p.name for p in self.root.iterdir()),
                                 ["annotations.jsonl.gz"] if compressed else
                                 ["annotations.jsonl", "annotations.jsonl.gz"])

    def test_missing_reference_annotation_is_not_silently_empty(self):
        path = self.write_annotations([self.annotation()])
        with self.assertRaisesRegex(ValueError, r":EOF: missing annotations.*R2"):
            read_annotations(path, {"R1": "CASS", "R2": "CA"})

    def test_duplicate_selected_id_is_rejected_even_after_all_ids_found(self):
        path = self.write_annotations([self.annotation(), self.annotation()])
        with self.assertRaisesRegex(ValueError, r":2: duplicate annotations.*R1"):
            read_annotations(path, {"R1": "CASS"})

    def test_annotation_sequence_must_match_reference(self):
        path = self.write_annotations([self.annotation(sequence="CAS")])
        with self.assertRaisesRegex(ValueError, r":1: reference R1: annotation cdr3"):
            read_annotations(path, {"R1": "CASS"})

    def test_invalid_record_structures_are_rejected(self):
        cases = [None, [], ["not a record"], [{"source_line": True, "fields": {}}],
                 [{"source_line": 2, "fields": {"cdr3": "CASS", "score": 0}}]]
        for records in cases:
            with self.subTest(records=records):
                path = self.write_annotations([{"reference_id": "R1", "records": records}])
                with self.assertRaisesRegex(ValueError, r":1: reference R1:"):
                    read_annotations(path, {"R1": "CASS"})

    def test_invalid_json_and_damaged_gzip_have_path_context(self):
        path = self.root / "invalid.jsonl"
        path.write_text("{invalid json}\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, r"invalid.jsonl:1:"):
            read_annotations(path, {"R1": "CASS"})
        path = self.root / "damaged.jsonl.gz"
        path.write_bytes(b"not gzip data")
        with self.assertRaisesRegex(ValueError, r"damaged.jsonl.gz:"):
            read_annotations(path, {"R1": "CASS"})


if __name__ == "__main__":
    unittest.main()
