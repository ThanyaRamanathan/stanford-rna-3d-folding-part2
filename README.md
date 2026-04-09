# 🧬 Stanford RNA 3D Folding — Part 2

[![Open in Kaggle](https://kaggle.com/static/images/open-in-kaggle.svg)](https://www.kaggle.com/code/thanyaramanathan/stanford-rna-3d-folding-pt2-v3)
[![Competition](https://img.shields.io/badge/Kaggle-Competition-20BEFF?logo=kaggle)](https://www.kaggle.com/competitions/stanford-rna-3d-folding-2)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

> Predicting RNA 3D structure from sequence using Template-Based Modeling (TBM) + [RNAPro](https://github.com/NVIDIA-Digital-Bio/RNAPro) — NVIDIA's state-of-the-art RNA folding model built with winners of the Stanford RNA 3D Folding Kaggle competition.

---

## 🔬 Problem

RNA molecules fold into precise 3D shapes that determine their biological function — from gene regulation to drug targeting. Unlike proteins (where AlphaFold transformed the field), RNA 3D structure prediction remains largely unsolved due to limited experimental data and the complexity of RNA folding dynamics.

**Part 2** raises the bar significantly over Part 1:

| Challenge | Part 1 | Part 2 |
|---|---|---|
| Target types | Single-chain RNA | Multi-chain complexes, RNA machines, protein-bound RNA |
| Max sequence length | ~1,000 nt | ~5,500 nt |
| Templates available | Most targets had templates | Includes template-free novel folds |
| Evaluation metric | TM-score (sequence-agnostic) | TM-score with **matching residue numbering** (stricter) |

Submissions are scored using TM-score averaged over the **best of 5 predictions** per target. The public leaderboard `best_template_oracle` baseline is **0.554**.

---

## 🏗️ Architecture

The pipeline combines two stages: a **template search** that finds structurally similar known RNAs, and **RNAPro** — a 500M-parameter diffusion model that uses those templates as structural priors.

```
test_sequences.csv
       │
       ▼
┌──────────────────────────────┐
│   Template-Based Modeling    │  ← searches training + validation
│   (TBM / PDB template search)│    structures for similar sequences
└──────────────┬───────────────┘
               │  submission_tbm.csv  (x,y,z coords per residue)
               ▼
┌──────────────────────────────┐
│  convert_templates_to_pt.py  │  ← converts CSV → binary .pt format
└──────────────┬───────────────┘
               │  templates.pt
               ▼
┌──────────────────────────────────────────────────────────────┐
│                         RNAPro (500M)                        │
│                                                              │
│  RibonanzaNet2 encoder  →  Pairformer  →  Diffusion sampler │
│  (RNA sequence features)   (template     (N_sample=5,        │
│                             embedder,     N_step=200,         │
│                             MSA)          N_cycle=10)         │
└──────────────────────────────┬───────────────────────────────┘
                               │
                               ▼
                        submission.csv
               (5 predicted structures per target)
```

**Key components:**
- **RNAPro** — post-trained Protenix (open AlphaFold 3 reproduction) with RNA-specific template embedder and RibonanzaNet2 as a frozen sequence encoder. Paper: [bioRxiv 2025.12.30](https://pmc.ncbi.nlm.nih.gov/articles/PMC12776560/)
- **RibonanzaNet2** — RNA foundation model pretrained on ~23.7M RNA sequences, used to extract sequence and pairwise features
- **Template embedder** — 2-block Pairformer that injects structural template information into the diffusion process
- **TBM fallback** — for sequences > `max_len` that RNAPro skips due to GPU memory limits

---

## 📊 Score Progression

| Version | Public LB Score | Key Changes |
|---|---|---|
| **v2** (baseline) | **0.354** | Custom BioPython TBM (global alignment, train data only), `N_sample=1`, `max_len=1000` |
| v3 | 0.302 ❌ | Switched to local alignment + relaxed length filter + validation templates → degraded TBM quality fed into RNAPro |
| **v4** (in progress) | **TBD** | jaejohn PDB-searched templates for Part 2 targets, `N_sample=5`, `max_len=1500` |

**What v3 taught us:** Template quality fed into RNAPro's embedder is the single most important variable. Switching to local (Smith-Waterman) alignment and relaxing the sequence length filter caused the TBM to select weaker, noisier templates — and since those feed directly into RNAPro as structural priors, the entire pipeline degraded. Better-sounding algorithm design ≠ better results without empirical validation at each stage.

---

## 🔑 Key Findings

### 1. Template quality dominates everything
RNAPro's template embedder was trained on high-quality PDB-searched templates (via MMseqs2 or the `jaejohn` BLAST pipeline). Feeding it weak BioPython-aligned templates from training data only produces noisy structural priors and actively hurts prediction quality vs. using no templates at all.

**Lesson:** Use the right TBM tool. The winning Part 1 approach searched the entire PDB using proper sequence homology tools, not just the competition training set.

### 2. The best-of-5 metric rewards diversity — but only with meaningful base predictions
Setting `N_sample=5` generates 5 distinct diffusion trajectories and is theoretically free points — the scorer takes the best. However, if the template is bad, all 5 samples collapse around the same wrong structure. Sample diversity only helps when base prediction quality is already reasonable.

### 3. Part 2's stricter residue-number matching is punishing for naive TBM
Part 1's TM-score allowed flexible residue alignment. Part 2 requires matching residue numbers between prediction and ground truth. Template coordinate transplanting — which naturally shifts residue indices — must carefully preserve the original numbering scheme. A perfectly shaped prediction with wrong resids scores zero.

### 4. Multi-chain targets need a fundamentally different approach
The custom TBM in this notebook treats all sequences as single-chain. Multi-chain complexes (new in Part 2) require chain-aware template matching and coordinate assembly — a major open gap vs. the Part 1 baseline.

---

## 🗂️ Repository Structure

```
stanford-rna-3d-folding/
├── README.md
├── notebooks/
│   └── stanford-rna-3d-folding-pt2-v3.ipynb   ← main notebook (code only)
├── results/
│   └── score_progression.png                   ← leaderboard screenshots
└── requirements.txt
```

> **Note on data:** All datasets (competition data, model checkpoints, CCD cache) are hosted on Kaggle and referenced by the notebook. They are too large to include in this repository. See the [Data Sources](#-data-sources) section below for Kaggle links.

---

## ⚡ Quickstart

The recommended way to run this is directly on Kaggle, where all data sources are pre-linked:

**1. Open the notebook on Kaggle**

[![Open in Kaggle](https://kaggle.com/static/images/open-in-kaggle.svg)](https://www.kaggle.com/code/thanyaramanathan/stanford-rna-3d-folding-pt2-v3)

**2. Add the required datasets** (if forking):

| Dataset | Kaggle Path |
|---|---|
| Competition data | `stanford-rna-3d-folding-2` |
| RNAPro source + checkpoint | `theoviel/rnapro-src` |
| CCD cache | `jaejohn/rnapro-ccd-cache` |
| RibonanzaNet2 | `shujun717/ribonanzanet2` |
| Python packages (offline) | `thanyaramanathan/pypackage` |

**3. Enable GPU** (T4 or better) and run all cells.

---

## 🧪 Running Locally (Advanced)

To run RNAPro outside Kaggle, follow the [official NVIDIA setup](https://github.com/NVIDIA-Digital-Bio/RNAPro):

```bash
git clone https://github.com/NVIDIA-Digital-Bio/RNAPro
cd RNAPro
conda create -n rnapro python=3.12 -y
conda activate rnapro
pip install -r requirements.txt
pip install -e .
```

Then generate templates using either:
- [MMseqs2 template identification](https://www.kaggle.com/code/rhijudas/mmseqs2-3d-rna-template-identification) (fast, good for large-scale PDB search)
- [jaejohn TBM-only approach](https://www.kaggle.com/code/jaejohn/rna-3d-folds-tbm-only-approach) (generally stronger — 1st place Part 1 method)

Convert to RNAPro format:
```bash
python preprocess/convert_templates_to_pt_files.py \
  --input_csv path/to/submission.csv \
  --output_name templates.pt \
  --max_n 40
```

---

## 📦 Data Sources

| Source | Description | Link |
|---|---|---|
| Competition data | Test sequences, training labels, validation labels, MSAs | [Kaggle](https://www.kaggle.com/competitions/stanford-rna-3d-folding-2/data) |
| RNAPro | Model source code + 500M private-best checkpoint | [Kaggle Dataset](https://www.kaggle.com/datasets/theoviel/rnapro-src) |
| RibonanzaNet2 | Pretrained RNA foundation model (encoder) | [Kaggle Model](https://www.kaggle.com/models/shujun717/ribonanzanet2) |
| CCD Cache | Chemical Component Dictionary for RNA residues | [Kaggle Dataset](https://www.kaggle.com/datasets/jaejohn/rnapro-ccd-cache) |
| jaejohn TBM | 1st-place Part 1 PDB template search output | [Kaggle Notebook](https://www.kaggle.com/code/jaejohn/rna-3d-folds-tbm-only-approach) |
| MMseqs2 Templates | MMseqs2-based PDB template identification | [Kaggle Notebook](https://www.kaggle.com/code/rhijudas/mmseqs2-3d-rna-template-identification) |

---

## 📚 References

- **RNAPro paper:** Lee Y. et al. "Template-based RNA structure prediction advanced through a blind code competition." *bioRxiv* (2025). [PMC link](https://pmc.ncbi.nlm.nih.gov/articles/PMC12776560/)
- **RNAPro code:** [github.com/NVIDIA-Digital-Bio/RNAPro](https://github.com/NVIDIA-Digital-Bio/RNAPro)
- **RibonanzaNet2:** He S. et al. *Nature Methods* (2024). [Paper](https://www.nature.com/articles/s41592-024-02487-0)
- **Protenix (AlphaFold 3 reproduction):** [github.com/bytedance/Protenix](https://github.com/bytedance/Protenix)
- **Competition:** [kaggle.com/competitions/stanford-rna-3d-folding-2](https://www.kaggle.com/competitions/stanford-rna-3d-folding-2)

---

## 🙏 Acknowledgements

This work builds directly on:
- The **NVIDIA Digital Biology** team for open-sourcing RNAPro and its weights
- **jaejohn** (G John Rao) for the 1st-place Part 1 TBM pipeline
- **theoviel** for the reference RNAPro inference notebook
- **Stanford Das Lab** and **HHMI** for organizing the competition and releasing experimental structures
- The **RibonanzaNet2** team for the RNA foundation model

---

<p align="center">
  <sub>Built for the Stanford RNA 3D Folding Part 2 Kaggle Competition · 2026</sub>
</p>
