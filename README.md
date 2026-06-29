# Meso2P

This repository hosts the core data analysis code and model files used in the study described in the paper:

> **Cortex-Wide, Cellular-Resolution Volumetric Imaging with a Modular Two-Photon Imaging Platform**  
> Jiahao Hu, Yanfeng Zhu, Shoupei Liu, Chengyu Li, Min Zhang, Xinyang Gu, Jingchuan Wu, Fang Xu, Ying Mao, Bo Li (2025)  
> *BioRxiv* 662899

## Contents

The code is organized into two main modules, corresponding to the two key experimental paradigms in the paper:

- **`Anesthesia-propagation`** – Analysis scripts for the anesthesia experiment, tracking cortical activity propagation under anesthetic conditions.
- **`Whisker-dynamics`** – Analysis scripts for the air-puff whisker stimulation experiment, characterizing whisker‑evoked dynamic responses.

Each module contains fully documented Jupyter notebooks/Python scripts, along with pre-trained models and intermediate data processing pipelines to reproduce the figures and statistical results reported in the publication.

## Getting Started

Please refer to the individual module READMEs for detailed usage instructions, dependency lists, and step‑by‑step reproduction guides.

## Citation

If you use this package in your research, please cite the original paper:

```bibtex
@article {Hu2025Meso2P,
  author = {Jiahao Hu and Yanfeng Zhu and Shoupei Liu and Chengyu Li and Min Zhang and Xinyang Gu and Jingchuan Wu and Fang Xu and Ying Mao and Bo Li},
  title = {Cortex-Wide, Cellular-Resolution Volumetric Imaging with a Modular Two-Photon Imaging Platform},
  year = {2025},
  journal = {BioRxiv},
  doi = {10.1101/662899}
}
