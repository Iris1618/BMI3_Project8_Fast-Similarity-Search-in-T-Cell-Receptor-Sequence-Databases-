# BMI3_Project8_Fast-Similarity-Search-in-T-Cell-Receptor-Sequence-Databases-
An indexed search algorithm that retrieves T-cell receptor  sequences within a specified edit distance of a query, while reducing the number of  full sequence comparisons.

## Input preparation

See [draft/Input_Source/README.md](draft/Input_Source/README.md) for VDJdb cleaning, configurable
reference subsets, synthetic datasets, common GQE queries and English pseudocode.

## Draft code

The initial module implementations are archived in [draft/](draft/):

```text
draft/
├── Verification/
│   ├── verification.py
│   └── verification_basic.py
├── Exhaustive_Benchmark/
├── Input_Source/
└── Inverted_Index/
```

Source files are preserved as originally written.
