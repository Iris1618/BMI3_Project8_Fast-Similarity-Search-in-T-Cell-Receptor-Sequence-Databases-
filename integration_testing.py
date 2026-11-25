"""BMI3 Project 8 集成测试入口。

定位：这里不做单个函数的白盒测试，而是把 Input_Source、
Inverted_Index、两种编辑距离验证器和 Exhaustive_Benchmark 串起来，
验证它们组合后仍能得到与穷举基线一致的结果。各目录下的
test_*.py 仍然负责单元级测试。

运行方式（在仓库根目录）：

    python test.py -v
    python -m unittest test -v

前置条件：Python 3.10+；Input_Source/examples 下的示例 TSV 存在。
测试只依赖示例数据和临时生成的合成数据，不依赖真实 VDJdb 下载。
所有写操作都限制在 tempfile.TemporaryDirectory() 中。
"""

import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

WORKSPACE = Path(__file__).resolve().parent
# 代码目录名通常是固定的，但允许克隆目录带后缀，因此先按完整名查找，
# 找不到再用 glob 兜底，避免 test.py 因目录改名而直接导入失败。
PROJECT = (
    WORKSPACE
    / "BMI3_code"
    / "BMI3_Project8_Fast-Similarity-Search-in-T-Cell-Receptor-Sequence-Databases--main"
)
if not PROJECT.is_dir():
    matches = sorted((WORKSPACE / "BMI3_code").glob("BMI3_Project8*"))
    if not matches:
        raise RuntimeError("找不到 BMI3_Project8 代码目录")
    PROJECT = matches[0]
INDEX_DIR = PROJECT / "Inverted_Index"
# 必须同时加入两个路径：
# - PROJECT：让 Exhaustive_Benchmark / Input_Source 以包的形式导入；
# - INDEX_DIR：integration_example.py 内部使用扁平导入 inverted_index，
#   只有把该目录放进 sys.path 才能复用它的真实实现。
for module_path in (str(PROJECT), str(INDEX_DIR)):
    if module_path not in sys.path:
        sys.path.insert(0, module_path)

# 注意：索引侧使用扁平导入，评测侧使用包导入，这与各模块自身的导入方式
# 保持一致，避免为了测试去修改生产代码。
from Exhaustive_Benchmark.baseline import exhaustive_matches, search_one as baseline_search_one
from Exhaustive_Benchmark.compare import (
    benchmark,
    load_truth,
    main as compare_main,
    memory_comparison,
)
from Input_Source.prepare_real import prepare as prepare_real
from Input_Source.prepare_synthetic import prepare as prepare_synthetic
from integration_example import search_one as index_search_one
from inverted_index import InvertedIndex, read_sequences
from verification import verify
from verification_basic import verification_basic

EXAMPLES = PROJECT / "Input_Source" / "examples"
REFERENCE_TSV = EXAMPLES / "reference.tsv"
QUERIES_TSV = EXAMPLES / "queries.tsv"
STRATEGIES = ("spread", "prefix")


def basic_verifier(query, reference, max_distance):
    """把基础 DP 适配成与 verify 相同的返回约定。

    verification_basic 总是返回 int，而 verification.verify 在距离超过
    阈值时返回 None。统一成 None 后，两个验证器才能共用 run_indexed，
    测试比较的也就是“最终命中集合和距离”，而不是各自的返回值风格。
    """
    distance = verification_basic(query, reference)
    return distance if distance <= max_distance else None


def run_indexed(index, query, d, verifier, strategy="spread"):
    """复现线上检索链路：生成候选 -> 只对候选跑 DP -> 收集 ≤ d 的命中。

    返回 (CandidateResult, {reference_id: distance})。
    测试中把它与穷举结果比较，即可同时检查：
    - 候选集是否漏掉了真实命中（假阴性）；
    - DP 给出的距离是否正确。
    """
    candidate = index.generate_candidates(query, d, strategy)
    matches = {}
    for reference_id in candidate.candidate_ids:
        distance = verifier(query, index.references[reference_id], d)
        if distance is not None and distance <= d:
            matches[reference_id] = distance
    return candidate, matches


def read_rows(path):
    """读取带表头的 TSV，保留所有原始列，用于检查类别、来源等元数据。"""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ExamplePipelineTests(unittest.TestCase):
    """针对仓库示例数据集的端到端测试。

    示例集很小（14 条参考、8 条查询），但覆盖普通查询、短查询回退、
    零命中和多命中，适合做全参数组合的交叉验证。
    """

    @classmethod
    def setUpClass(cls):
        # 整个测试类共享同一份只读输入，避免每个用例重复读盘。
        cls.references = read_sequences(REFERENCE_TSV, "reference_id")
        cls.queries = read_sequences(QUERIES_TSV, "query_id")

    def test_all_thresholds_and_strategies_match_exhaustive(self):
        """核心正确性：所有 d/k/种子策略组合都必须与穷举基线一致。

        验证逻辑分为四层：
        1. 两个验证器给出的 {参考ID: 距离} 与穷举结果完全相同；
        2. 真实命中必须包含在候选集中（候选只允许多、不允许漏）；
        3. 候选必须满足长度差过滤；
        4. 索引模式必须符合鸽巢原理的启用条件。
        """
        for d in (0, 1, 2):
            for k in (2, 3):
                # k 会改变索引内容；d 只在查询时参与候选生成和验证。
                index = InvertedIndex(self.references, k)
                self.assertEqual(index.k, k)
                for strategy in STRATEGIES:
                    for query_id, query in self.queries.items():
                        # subTest 让失败信息精确到具体的 d/k/策略/查询。
                        with self.subTest(d=d, k=k, strategy=strategy, query=query_id):
                            # 独立真值：对全库做长度过滤 + 基础 DP，不经过倒排索引。
                            expected = exhaustive_matches(query, self.references, d)
                            candidate, got_basic = run_indexed(
                                index, query, d, basic_verifier, strategy
                            )
                            _, got_optimized = run_indexed(index, query, d, verify, strategy)

                            # 两种 DP 实现必须给出相同的命中集合和距离。
                            self.assertEqual(got_basic, expected)
                            self.assertEqual(got_optimized, expected)
                            # 候选只允许是真实命中的超集，不能漏召回。
                            self.assertLessEqual(set(expected), set(candidate.candidate_ids))
                            # 所有候选都必须通过长度差过滤。
                            self.assertTrue(
                                all(
                                    abs(len(self.references[rid]) - len(query)) <= d
                                    for rid in candidate.candidate_ids
                                )
                            )
                            # 查询长度不足以切出 d+1 段长度 k 的互不重叠种子时，
                            # 鸽巢保证失效，必须回退到全扫；否则应走索引。
                            expected_mode = (
                                "fallback" if len(query) < k * (d + 1) else "indexed"
                            )
                            self.assertEqual(candidate.mode, expected_mode)

    def test_final_result_helpers_agree(self):
        """最终输出契约：索引版和穷举版的 search_one 必须逐字段一致。

        这条用例不比较耗时，只确认两个对外结果构造函数在字段、注释、
        距离和排序上完全等价，避免集成阶段“换了入口就换结果”。
        """
        index = InvertedIndex(self.references, 2)
        # 故意只给一条参考加注释，检查注释是否原样带到最终命中里，
        # 而不是在索引路径中被丢掉。
        annotations = {"EDGE_R05": [{"source_line": 5, "fields": {}}]}
        for query_id, query in self.queries.items():
            with self.subTest(query=query_id):
                indexed = index_search_one(
                    index, query_id, query, 1, basic_verifier, annotations=annotations
                )
                baseline = baseline_search_one(
                    query_id, query, self.references, 1, annotations
                )
                # 比较整个 dict，包括 matches、query 元信息等所有字段。
                self.assertEqual(indexed, baseline)
                # 输出顺序也是契约的一部分，避免下游 diff 出现无意义抖动。
                self.assertEqual(
                    indexed["matches"],
                    sorted(indexed["matches"], key=lambda hit: hit["reference_id"]),
                )

    def test_candidate_cli_and_overwrite_guard(self):
        """真实 CLI 集成 + 输出安全。

        这里刻意用子进程而不是直接调用 main()，因为 CLI 自己会设置
        sys.path、解析参数并写 JSONL；只有子进程才能覆盖这条真实路径。
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "candidates.jsonl"
            # 参数与 README 中的推荐用法一致，输出写到临时目录。
            command = [
                sys.executable,
                str(PROJECT / "Inverted_Index" / "inverted_index.py"),
                "--reference", str(REFERENCE_TSV),
                "--queries", str(QUERIES_TSV),
                "--k", "2",
                "--d", "1",
                "--strategy", "spread",
                "--output", str(output),
            ]
            completed = subprocess.run(
                command, capture_output=True, text=True, check=True
            )
            rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
            self.assertEqual(len(rows), len(self.queries))
            # CLI 把汇总 JSON 写到 stderr，最后一行是机器可读的摘要。
            summary = json.loads(completed.stderr.strip().splitlines()[-1])
            self.assertEqual(summary["query_count"], len(self.queries))
            self.assertEqual(summary["reference_count"], len(self.references))
            # 每条查询要么走索引，要么回退，两者之和必须等于查询总数。
            self.assertEqual(
                summary["indexed_queries"] + summary["fallback_queries"],
                len(self.queries),
            )

            # 用与 CLI 相同的 k 重建索引，逐条复核候选和最终命中。
            cli_index = InvertedIndex(self.references, 2)
            for row in rows:
                with self.subTest(query=row["query_id"]):
                    query = self.queries[row["query_id"]]
                    expected = exhaustive_matches(query, self.references, row["d"])
                    # CLI 给出的候选必须覆盖所有真实命中。
                    self.assertLessEqual(set(expected), set(row["candidate_ids"]))
                    _, verified = run_indexed(
                        cli_index,
                        query,
                        row["d"],
                        verify,
                        row["strategy"],
                    )
                    self.assertEqual(verified, expected)

            # 覆盖保护：--output 与输入文件是同一路径时必须拒绝执行。
            # 先复制到临时目录，避免触碰仓库里的示例文件。
            guarded_reference = root / "guarded_reference.tsv"
            shutil.copyfile(REFERENCE_TSV, guarded_reference)
            original = guarded_reference.read_bytes()
            failed = subprocess.run(
                [
                    sys.executable,
                    str(PROJECT / "Inverted_Index" / "inverted_index.py"),
                    "--reference", str(guarded_reference),
                    "--queries", str(QUERIES_TSV),
                    "--k", "2",
                    "--d", "1",
                    "--output", str(guarded_reference),
                ],
                capture_output=True,
                text=True,
            )
            # 退出码 2 是该 CLI 约定的“输入/参数错误”。
            self.assertEqual(failed.returncode, 2)
            self.assertIn("must not overwrite an input file", failed.stderr)
            # 关键：拒绝执行后文件内容必须一字节不变。
            self.assertEqual(guarded_reference.read_bytes(), original)

    def test_compare_validate_with_generated_truth(self):
        """compare.py validate 模式 + 外部真值文件加载。

        真值完全由穷举基线独立生成，而不是从索引结果反推，
        因此这一步能发现“索引和验证器一起错”的系统性偏差。
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth_path = root / "truth_matches.tsv"
            summary_path = root / "truth_summary.tsv"
            output_path = root / "validation.csv"
            reference_set = "examples"

            # truth_matches.tsv：每个查询的每个命中一行，距离来自穷举 DP。
            with open(truth_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["reference_set", "threshold", "query_id", "reference_id", "distance"]
                )
                for query_id, query in self.queries.items():
                    for reference_id, distance in sorted(
                        exhaustive_matches(query, self.references, 1).items()
                    ):
                        writer.writerow(
                            [reference_set, 1, query_id, reference_id, distance]
                        )

            # truth_summary.tsv：每个查询一行，用命中数和参考库大小
            # 校验真值文件是否完整。零命中查询也必须写 n_matches=0。
            with open(summary_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["reference_set", "threshold", "query_id", "n_matches", "n_reference"]
                )
                for query_id, query in self.queries.items():
                    writer.writerow(
                        [
                            reference_set,
                            1,
                            query_id,
                            len(exhaustive_matches(query, self.references, 1)),
                            len(self.references),
                        ]
                    )

            # 先单独验证加载器：它应校验 reference_set、阈值、
            # 参考库大小和命中数，并返回 {query_id: {reference_id: distance}}。
            loaded = load_truth(
                truth_path, summary_path, reference_set,
                self.queries, 1, len(self.references),
            )
            for query_id, query in self.queries.items():
                self.assertEqual(
                    loaded[query_id], exhaustive_matches(query, self.references, 1)
                )

            # 再走完整的 validate CLI：索引两种 DP + 真值比对。
            return_code = compare_main(
                [
                    "validate",
                    "--reference", str(REFERENCE_TSV),
                    "--queries", str(QUERIES_TSV),
                    "--d", "1",
                    "--k", "2",
                    "--strategy", "spread",
                    "--truth", str(truth_path),
                    "--truth-summary", str(summary_path),
                    "--reference-set", reference_set,
                    "--output", str(output_path),
                ]
            )
            self.assertEqual(return_code, 0)

            with open(output_path, encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            # 每个查询 3 行：indexed_basic、indexed_optimized、真值比对。
            self.assertEqual(len(rows), 3 * len(self.queries))
            self.assertEqual(
                {row["comparison"] for row in rows},
                {
                    "indexed_basic",
                    "indexed_optimized",
                    "exhaustive_vs_supplied_truth",
                },
            )
            for row in rows:
                with self.subTest(query=row["query_id"], comparison=row["comparison"]):
                    # bool 写入 CSV 后是字符串 "True"。
                    self.assertEqual(row["match"], "True")
                    # 三类差异列表都应为空，表示无漏报、无多报、距离无误。
                    self.assertEqual(row["missing_ids"], "[]")
                    self.assertEqual(row["unexpected_ids"], "[]")
                    self.assertEqual(row["wrong_distances"], "[]")
                    self.assertEqual(
                        row["expected_matches"], row["observed_matches"]
                    )

    def test_benchmark_and_memory_modes(self):
        """计时与内存评测的接口契约（不断言具体耗时）。

        耗时受机器负载影响，断言绝对数值会让测试变得不稳定；
        这里只检查评测器是否产出了完整的结构、样本和方法集合。
        """
        # d=1, k=2, spread, repeats=2：repeats 至少为 2 才能计算标准误。
        metadata, rows = benchmark(self.references, self.queries, 1, 2, "spread", 2)
        # 元数据必须记录输入规模、查询顺序、建索引样本和中位数，
        # 这些是后续对比不同数据集时的可追溯信息。
        self.assertEqual(metadata["reference_count"], len(self.references))
        self.assertEqual(metadata["query_count"], len(self.queries))
        self.assertEqual(metadata["query_ids"], list(self.queries))
        self.assertEqual(len(metadata["index_build_ms_samples"]), 2)
        self.assertGreater(metadata["index_build_median_ms"], 0)
        # 5 类方法：3 个整批搜索 + 2 个在相同候选对上的验证器对比。
        self.assertEqual(
            {row["method"] for row in rows},
            {
                "exhaustive_basic",
                "indexed_basic",
                "indexed_optimized",
                "candidate_basic",
                "candidate_optimized",
            },
        )
        for row in rows:
            with self.subTest(method=row["method"]):
                # 每次重复都应留下一个原始样本，方便复核离散程度。
                self.assertEqual(len(json.loads(row["samples"])), 2)
                if row["method"].startswith("candidate_"):
                    # 候选对方法按“每对候选”计时，单位与整批搜索不同。
                    self.assertGreater(row["observations_per_repeat"], 0)
                    self.assertIsNotNone(row["median"])
                else:
                    self.assertEqual(row["unit"], "ms_per_query")
                    self.assertIsNotNone(row["mean"])
                    self.assertIsNotNone(row["median"])

        # 内存对比使用 tracemalloc；索引方法包含建索引成本。
        # 指标是 Python 增量峰值，不是进程 RSS，也不包含已预加载的数据。
        memory_rows = memory_comparison(
            self.references, self.queries, 1, 2, "spread"
        )
        self.assertEqual(
            {row["method"] for row in memory_rows},
            {"exhaustive_basic", "indexed_basic", "indexed_optimized"},
        )
        for row in memory_rows:
            self.assertGreater(row["incremental_peak_python_bytes"], 0)
            # scope 必须明确声明统计口径，避免被误当成进程内存。
            self.assertIn("not RSS", row["scope"])


class InputPreparationTests(unittest.TestCase):
    """输入准备模块的测试：合成数据生成 + 真实数据清洗。"""

    def test_synthetic_generation_and_search(self):
        """合成数据生成器 + 索引检索的联合验证。

        使用固定种子生成 3 个 GQE motif 比例的库，因此结果可复现；
        生成后不信任生成器自报的类别，而是重新用穷举基线核对。
        """
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "synthetic"
            # 只构造 prepare() 需要的参数，避免走 argparse 和文件发布逻辑。
            args = SimpleNamespace(
                out=output,
                size=40,
                lengths=[12],
                motif_percentages=[0, 50, 90],
                anchors=3,
                max_edits=2,
                probe_count=2,
                seed=20262111,
                probe_seed=20261113,
            )
            prepare_synthetic(args, output)

            # 顶层报告应列出 3 个数据集；每个数据集还会各有一份报告。
            root_report = json.loads(
                (output / "generation_report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(len(root_report["datasets"]), 3)

            for percentage in (0, 50, 90):
                name = f"L12_motif{percentage:03d}"
                library = output / name
                references = read_sequences(library / "reference.tsv", "reference_id")
                queries = read_sequences(library / "queries.tsv", "query_id")
                # rows 保留 category/source_id 等元数据，用于后面的分类检查。
                rows = {row["query_id"]: row for row in read_rows(library / "queries.tsv")}

                self.assertEqual(len(references), 40)
                self.assertEqual(len(queries), len(rows))

                # k=3，d=1/2；候选 + 带状 DP 必须与穷举基线完全一致。
                index = InvertedIndex(references, 3)
                for d in (1, 2):
                    for query_id, query in queries.items():
                        with self.subTest(library=name, d=d, query=query_id):
                            expected = exhaustive_matches(query, references, d)
                            candidate, got = run_indexed(index, query, d, verify)
                            self.assertEqual(got, expected)
                            self.assertLessEqual(set(expected), set(candidate.candidate_ids))

                for query_id, row in rows.items():
                    if row["category"] == "exact_copy":
                        # 精确拷贝：应与源序列相同，且在 d=0 时候选就包含源 ID。
                        self.assertEqual(row["sequence"], row["source_sequence"])
                        self.assertIn(row["source_id"], references)
                        self.assertEqual(references[row["source_id"]], row["sequence"])
                        candidate = index.generate_candidates(row["sequence"], 0)
                        self.assertIn(row["source_id"], candidate.candidate_ids)
                    if row["category"] == "short_query":
                        # 短查询（长度 1-3）不足以切出 d+1 段长度 3 的种子，
                        # 必须触发全库回退。
                        candidate = index.generate_candidates(row["sequence"], 1)
                        self.assertEqual(candidate.mode, "fallback")

            # motif probe 是人为构造的“接近但不等于任何参考”的查询。
            probes_path = output / "motif_probes" / "L12" / "queries.tsv"
            self.assertTrue(probes_path.is_file())
            probes = read_rows(probes_path)
            self.assertEqual(len(probes), 2)
            # 跨三个库收集序列，确保 probe 不是任何库的精确拷贝。
            all_reference_sequences = set()
            for percentage in (0, 50, 90):
                all_reference_sequences.update(
                    read_sequences(
                        output / f"L12_motif{percentage:03d}" / "reference.tsv",
                        "reference_id",
                    ).values()
                )
            for probe in probes:
                with self.subTest(probe=probe["query_id"]):
                    # probe 不是精确命中，但与源序列恰好差一次编辑。
                    self.assertNotIn(probe["sequence"], all_reference_sequences)
                    self.assertEqual(probe["source_distance"], "1")
                    self.assertEqual(
                        verification_basic(probe["source_sequence"], probe["sequence"]), 1
                    )

    def test_real_data_preparation(self):
        """真实数据清洗器的端到端检查（使用伪造的 VDJdb 小样本）。

        输入覆盖重复、非法氨基酸、空 CDR3、错误链、错误物种五类情况；
        不需要下载真实 VDJdb，也能验证清洗、去重、注释保留和子集生成。
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "vdjdb.slim.txt"
            # 7 行原始记录：2 行同序列、1 行非法字符、1 行空 CDR3、
            # 1 行错误链、1 行错误物种。
            source.write_text(
                "gene\tspecies\tcdr3\tweight\n"
                "TRB\tHomoSapiens\tCASS\t1\n"
                "TRB\tHomoSapiens\tCASS\t2\n"
                "TRB\tHomoSapiens\tCATS\t1\n"
                "TRB\tHomoSapiens\tCASSX\t1\n"
                "TRB\tHomoSapiens\t\t1\n"
                "TRA\tHomoSapiens\tCASS\t1\n"
                "TRB\tMusMusculus\tCASS\t1\n",
                encoding="utf-8",
            )
            output = root / "real"
            args = SimpleNamespace(
                input=source,
                out=output,
                sizes=[2],
                queries=2,
                seed=7,
                species="HomoSapiens",
                gene="TRB",
            )
            prepare_real(args, output)

            # 清洗报告必须准确区分“目标物种/链过滤后”和“通过 CDR3 检查”。
            report = json.loads(
                (output / "cleaning_report.json").read_text(encoding="utf-8")
            )
            counts = report["counts"]
            self.assertEqual(counts["raw_rows"], 7)
            self.assertEqual(counts["target_rows"], 5)
            self.assertEqual(counts["valid_rows"], 3)
            self.assertEqual(counts["unique_sequences"], 2)
            self.assertEqual(counts["duplicate_rows_merged"], 1)
            self.assertEqual(counts["missing_cdr3"], 1)
            self.assertEqual(counts["invalid_amino_acid"], 1)
            self.assertEqual(counts["excluded_species_or_chain"], 2)
            # 输入指纹用于实验可追溯，SHA-256 十六进制长度固定为 64。
            self.assertEqual(len(report["input_sha256"]), 64)

            # 去重后只保留 CASS 和 CATS；最小子集大小等于全量时两者一致。
            references = read_sequences(output / "reference_all.tsv", "reference_id")
            self.assertEqual(sorted(references.values()), ["CASS", "CATS"])
            subset = read_sequences(output / "reference_n2.tsv", "reference_id")
            self.assertEqual(subset, references)

            # 注释按参考 ID 聚合；重复的 CASS 应保留两条原始记录。
            with gzip.open(output / "annotations.jsonl.gz", "rt", encoding="utf-8") as handle:
                annotations = [json.loads(line) for line in handle]
            self.assertEqual(len(annotations), 2)
            records_by_id = {row["reference_id"]: row["records"] for row in annotations}
            cass_id = next(rid for rid, sequence in references.items() if sequence == "CASS")
            self.assertEqual(len(records_by_id[cass_id]), 2)

            # 查询应是参考序列的精确拷贝，距离为 0、来源 ID 可追溯。
            query_rows = read_rows(output / "queries.tsv")
            self.assertEqual(len(query_rows), 2)
            for row in query_rows:
                self.assertEqual(row["category"], "exact_copy")
                self.assertIn(row["source_id"], references)
                self.assertEqual(references[row["source_id"]], row["sequence"])
                self.assertEqual(row["source_distance"], "0")

            # 被排除的目标行只包含两类 CDR3 问题，不含物种/链过滤行。
            excluded = read_rows(output / "excluded_target_rows.tsv")
            self.assertEqual(
                sorted(row["reason"] for row in excluded),
                ["invalid_amino_acid", "missing_cdr3"],
            )
            # 两条长度 4 的序列都在 8-30 的推荐范围之外。
            self.assertEqual(len(read_rows(output / "length_review_flags.tsv")), 2)


if __name__ == "__main__":
    # 直接运行 test.py 时以详细模式输出每个用例名。
    unittest.main(verbosity=2)
