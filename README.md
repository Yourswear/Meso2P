# Meso2P

This repository hosts the core data analysis code and model files used in the study described in the paper:

> **Cortex-Wide, Cellular-Resolution Volumetric Imaging with a Modular Two-Photon Imaging Platform**  
> Jiahao Hu, Yanfeng Zhu, Shoupei Liu, Chengyu Li, Min Zhang, Xinyang Gu, Jingchuan Wu, Fang Xu, Ying Mao, Bo Li (2025)  
> *BioRxiv* 662899

## Contents

The code is organized into two main modules, corresponding to the two key experimental paradigms in the paper:

- **`Anesthesia-propagation`** – Analysis scripts for the anesthesia experiment, tracking cortical activity propagation under anesthetic conditions.
- **`Whisker-dynamics`** – Analysis scripts for the air-puff whisker stimulation experiment, characterizing whisker‑evoked dynamic responses.

Each module contains fully documented matlab scripts, along with intermediate data processing pipelines to reproduce the figures and statistical results reported in the publication.

## Getting Started

### Data Download
First, download the required test dataset from Zenodo: [10.5281/zenodo.20956973](https://doi.org/10.5281/zenodo.20956973).

### Whisker-dynamics (Air-puff stimulation)
1. Place the downloaded files `17-AH-all_result.mat` and `AP0,6mm_recode.mat` into the `Whisker-dynamics/` folder.
2. Run the script `Airpuff_analysis_for_submitt.m`.
3. The script will generate 16 figures:

- **Figure 1** – Atlas alignment calibration pattern  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%201.png" width="400px">
- **Figure 2** – Distribution of neuronal counts across different brain regions  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%202.png" width="400px">
- **Figure 3** – Calcium signal heatmap across 10 trials  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%203.png" width="400px"> 
- **Figure 4** – Spike rate statistics for control group and 10 trials  
   <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%204.png" width="400px">
- **Figure 5** – Normalized neuronal spike probabilities over time and position for trials 1 and 10. Bar plots show the normalized spike probability at each time point and anterior-posterior bin for the first (Trial 1) and last (Trial 10) trials.
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%205.png" width="400px">
- **Figure 6** – Spike probabilities across 10 trials in M region  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%206.png" width="400px">
- **Figure 7** – Spike probabilities across 10 trials in S region  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%207.png" width="400px">
- **Figure 8** – Spike probabilities across 10 trials in V region  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%208.png" width="400px">
- **Figure 9** – Spike probabilities across 10 trials in RSP region  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%209.png" width="400px">
- **Figure 10** – Cross‑correlation coefficient matrix between brain regions across 10 trials  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2010.png" height="150px">
- **Figure 11** – The projected activity of S region in the CCA dimensions identified across 10 trials
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2011.png" width="400px">
- **Figure 12** – The projected activity of M region in the CCA dimensions identified across 10 trials
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2012.png" width="400px">
- **Figure 13** – The projected activity of V region in the CCA dimensions identified across 10 trials 
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2013.png" width="400px">
- **Figure 14** – The projected activity of RSP region in the CCA dimensions identified across 10 trials
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2014.png" width="400px">
- **Figure 15** – Heatmap of neuronal loading weights across trials  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2015.png" width="400px">
- **Figure 16** – PC1 derived from PCA of neuronal loading weight matrices  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Whisker-dynamics/Figure%2016.png" width="400px">

### Anesthesia-propagation (Anesthesia experiment)
1. Place the downloaded file `17_induce_deep.mat` into the `Anesthesia-propagation/` folder (note: keep it in this module folder to maintain consistency with the repository structure).
2. Run the script `Anesthesia_analysis_for_submitt.m`.
3. The script will generate 4 figures:

- **Figure 1** – Mean spike probabilities during anesthesia induction  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Anethesia-propagation/Figure%201.png" height="150px">
- **Figure 2** – Normalized global response delay. Standard deviation (SD) of peak times in spike probabilities across 16 anterior-posterior (AP) bins was computed during anesthesia progression  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Anethesia-propagation/Figure%202.png" width="400px">
- **Figure 3** – Bar plot of normalized neuronal spike probability across anterior-posterior positions  
  <img src="https://github.com/Yourswear/Meso2P/blob/master/Anethesia-propagation/Figure%203.png" width="400px">

## Citation

If you use these codes in your research, please cite the original paper:

```bibtex
@article {Hu2025Meso2P,
  author = {Jiahao Hu and Yanfeng Zhu and Shoupei Liu and Chengyu Li and Min Zhang and Xinyang Gu and Jingchuan Wu and Fang Xu and Ying Mao and Bo Li},
  title = {Cortex-Wide, Cellular-Resolution Volumetric Imaging with a Modular Two-Photon Imaging Platform},
  year = {2025},
  journal = {BioRxiv},
  doi = {10.1101/662899}
}
