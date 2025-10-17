# diSBPred

*A machine-learning framework for predicting disulfide bonds directly from protein sequence.*

Authors: Avdesh Mishra, Md Wasi Ul Kabir, Md Tamjidul Hoque

---

## Overview

**diSBPred** predicts intra-protein disulfide bonds from primary sequence. It ships with a benchmark set and scripts to run natively or via Docker.

* **Input:** text file at `Input/input.txt`
* **Output:** result files written to `Output/`
* **Tested OS:** Ubuntu 20.04

If you use diSBPred in published work, please cite the paper listed in [Citation](#citation).

---

## Dataset

The `Dataset/` directory contains:

* `list_1859_Uniprot_Cutoff25_33_Intra_Cutoff25_Combined_Cutoff25_TotalProts1866.txt` — UniProt IDs used to benchmark diSBPred.

> Tip: use this list to reproduce reported benchmarks or to sanity-check your setup.

---

## Quick Start

```bash
# 1) Get the code
git clone https://github.com/wasicse/diSBPred.git
cd diSBPred

# 2) Install dependencies (pyenv, Python 3.7.4, Poetry 1.1.13)
./install_dependencies.sh

# 3) Run (reads Input/input.txt, writes to Output/)
./run_diSBPred.sh
```

---

## Installation

### Native (recommended for development)

**Prereqs**

* `pyenv` (latest)
* `Python 3.7.4` (installed via pyenv)
* `Poetry 1.1.13`

**One-liner:**

```bash
./install_dependencies.sh
```

This script installs/initializes pyenv, sets Python 3.7.4 locally, and installs project packages with Poetry.

### Docker (fastest to try)

Build locally:

```bash
docker build -t wasicse/disbpred https://github.com/wasicse/disbpred.git#master
```

—or pull the prebuilt image:

```bash
docker pull wasicse/disbpred:latest
```

Run (mounts the current repo’s `Input/` and `Output/`):

```bash
./run_diSBPred_Docker.sh "$(pwd)/Input/input.txt" Output
```

---

## Running diSBPred

### Inputs

* Default input path: `Input/input.txt`
* Format: one entry per line (see examples shipped in the repo).
* Ensure the file is readable and non-empty before running.

### Command

```bash
./run_diSBPred.sh
```

### Outputs

* Results are written to the `Output/` directory (created if missing).
* Check the console log for the exact output filenames produced by your run.

---

## Project Layout

```
diSBPred/
├─ Dataset/
│  └─ list_1859_Uniprot_..._TotalProts1866.txt
├─ Input/
│  └─ input.txt
├─ Output/                # created/populated after a run
├─ run_diSBPred.sh
├─ run_diSBPred_Docker.sh
├─ install_dependencies.sh
└─ (source and config files)
```

---

## Reproducibility Notes

* Validated on **Ubuntu 20.04**. Other Linux distros may work but are not guaranteed.
* For strict reproducibility, prefer the **Docker** workflow.

---

## Troubleshooting

**Poetry/pyenv not found**

* Re-open your shell to refresh PATH, or source your profile:

  ```bash
  source ~/.bashrc  # or ~/.zshrc
  ```

**Wrong Python version**

* Inside the repo:

  ```bash
  pyenv local 3.7.4
  poetry env use 3.7.4
  poetry install
  ```

**Permission errors on scripts**

```bash
chmod +x run_diSBPred.sh run_diSBPred_Docker.sh install_dependencies.sh
```

**Docker “permission denied” on mounts**

* On Linux, ensure your user can access the repo path and that `$(pwd)/Output` exists (or let the script create it).

---

## Support

Questions/bugs: Md Tamjidul Hoque — [thoque@uno.edu](mailto:thoque@uno.edu)

---

## Citation

Mishra, A., Kabir, M. W. U., & Hoque, M. T. (2021). *diSBPred: A Machine Learning Based Approach for Disulfide Bond Prediction.* **Computational Biology and Chemistry**, 91, 107436. [https://doi.org/10.1016/j.compbiolchem.2021.107436](https://doi.org/10.1016/j.compbiolchem.2021.107436)

**BibTeX**

```bibtex
@article{Mishra2021diSBPred,
  title   = {diSBPred: A Machine Learning Based Approach for Disulfide Bond Prediction},
  author  = {Mishra, Avdesh and Kabir, Md Wasi Ul and Hoque, Md Tamjidul},
  journal = {Computational Biology and Chemistry},
  volume  = {91},
  pages   = {107436},
  year    = {2021},
  doi     = {10.1016/j.compbiolchem.2021.107436}
}
```


