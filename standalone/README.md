# SINCOPA — Standalone Command-Line Version

SINCOPA detects selective sweeps in bacteria based on **S**imilarity of **IN**con**CO**ngruent **PA**tterns. It identifies regions of a DNA multiple sequence alignment where allelic patterns are highly incongruent with the species phylogeny, consistent with a transfer-mediated selective sweep.

SINCOPA is also available as a [web server](https://sincopa.tau.ac.il/).

## Requirements

- Python ≥ 3.6
- conda
- git ≥ 2.25

## Installation

```bash
git clone --filter=blob:none --sparse https://github.com/orenavram/SINCOPA.git
cd SINCOPA
git sparse-checkout set standalone
cd standalone
conda env create -f environment.yml
conda activate sincopa
```

## Usage

```bash
python sincopa.py <alignment.fasta> <tree.newick> <output_dir> [--window_size 50]
```

| Argument | Description |
|---|---|
| `alignment.fasta` | DNA multiple sequence alignment in FASTA format |
| `tree.newick` | Species phylogeny in Newick format |
| `output_dir` | Directory where results will be written |
| `--window_size` | Sliding window size (default: 50) |

**Input requirements:**
- Alignment must be ≥ 300 bp after trimming
- The species tree must contain all taxa present in the alignment
- The tree should be reconstructed from independent data, not from the input alignment

## Output

| File | Description |
|---|---|
| `homoplasy.txt` | Parsimony-based homoplasy score per alignment column |
| `sweeps_scores.txt` | S* score per sliding window |
| `sweeps_summary.txt` | Summary statistics including peak S* score |
| `sweeps_scores.png` | S* score plot across the alignment |
| `done.txt` | Created upon successful completion |

## Quick test

```bash
python sincopa.py test/example.fasta test/example.newick test/out --window_size 50
```
