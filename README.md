# WSINDyStructuredPopulations

This repository contains the code used to generate the figures in the paper:

> **Learning Structured Population Models from Data with WSINDy**  
> R. Lyons, V. Dukic, and D. M. Bortz  
> [arXiv:2506.24101](https://arxiv.org/abs/2506.24101)

---

## 📂 Contents

- `Main_artificialData.m` – reproduces results for artificial examples **L.1–L.4** from the paper.
- `Main_RealData.m` – processes and analyzes real population data from Jackson *et al.*

---

## 🔬 Artificial Data (L.1–L.4)

These examples use synthetic data for structured population models and can be run directly using the provided scripts.

⚠️ The datasets for the nonlinear examples (e.g., NL.1–NL.2) are too large to be hosted on GitHub.  
If you are interested in running these examples, please contact the authors to request access.

---

## 📊 Real Population Data

The script `Main_RealData.m` uses empirical demographic data from the following study:

> **Changes in age-structure over four decades were a key determinant of population growth rate in a long-lived mammal**  
> J. Jackson, M. Coulson, L. Börger, et al.  
> *Journal of Animal Ecology*  
> [DOI: 10.1111/1365-2656.13290](https://doi.org/10.1111/1365-2656.13290)

We include only the subset of the dataset necessary to reproduce our results, in accordance with the public domain license.  
The full dataset is available via Dryad:  
➡️ [https://datadryad.org/dataset/doi:10.5061/dryad.m905qftwx](https://datadryad.org/dataset/doi:10.5061/dryad.m905qftwx)

---

## 📜 License and Data Use

All included example data are either synthetic or derived from public domain sources.  
The Jackson *et al.* dataset is shared under its original public domain designation.

---
License:
Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)

You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material

Under the following terms:
- Attribution — You must give appropriate credit
- NonCommercial — You may not use the material for commercial purposes

Full license: https://creativecommons.org/licenses/by-nc/4.0/
---
