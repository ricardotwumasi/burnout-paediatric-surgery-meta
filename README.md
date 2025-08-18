# Burnout Prevalance in Paediatric Surgeons Meta-Analysis Code

[![DOI](https://zenodo.org/badge/DOI/[pending].svg)](https://doi.org/[pending])
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.1.0-blue.svg)](https://cran.r-project.org/)

## Version 1.0.2
**Updated:** August 2025
  
Full R code for reproducing our meta-analysis of the prevalance of burnout in Paediatric Surgeons. This repository contains the code used in our systematic review and meta-analysis (Kirdar-Smith et al., in preparation) 

## Registration  
**PROSPERO Registration:** [42025640570](https://www.crd.york.ac.uk/PROSPERO/view/CRD42025640570)

## Overview

This repository contains the complete code and data for a systematic review and meta-analysis examining burnout prevalence among paediatric surgeons. The analysis includes 15 studies capturing data from a total of 2757 paediatric surgeons.

## 💻 Requirements

- R (≥ 4.1.0)
- Required R packages:
  ```R
  tidyverse
  metafor
  ```

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/ricardotwumasi/burnout-paediatric-surgery-meta.git
   ```

2. Install required R packages:
   ```R
   required_packages <- c("tidyverse", "metafor")
   install.packages(required_packages)
   ```

## 🤖 AI Statement

This repo was vibe coded with the assistance of Claude Sonnet 4 (Anthropic, San Francisco: CA)

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citations

Key methodological references:

```bibtex
@book{harrer2021,
      title     = {Doing Meta-Analysis With {R}: A Hands-On Guide},
      author    = {Harrer, Mathias and Cuijpers, Pim and Furukawa Toshi A and Ebert, David D},
      year      = {2021},
      publisher = {Chapman & Hall/CRC Press},
      address   = {Boca Raton, FL and London},
      isbn      = {9780367610074},
      edition   = {1st}
    }

@book{borenstein2021,
	title = {Introduction to {Meta}-{Analysis}},
	isbn = {978-1-119-55835-4},
	url = {https://books.google.co.uk/books?id=2oYmEAAAQBAJ},
	publisher = {Wiley},
	author = {Borenstein, M. and Hedges, L.V. and Higgins, J.P.T. and Rothstein, H.R.},
	year = {2021},
}


@article{viechtbauer2010,
  title={Conducting meta-analyses in R with the metafor package},
  author={Viechtbauer, Wolfgang},
  journal={Journal of Statistical Software},
  volume={36},
  number={3},
  pages={1--48},
  year={2010}
}
```

[Citation](#citation) 
For citing this repository, please use:

<details>
<summary>BibTeX</summary>
<pre><code>@article{kirdar-smith2025,
  title={Burnout Prevalence in Paediatric Surgeons: A Systematic Review and Meta-Analysis},
  author={Kirdar-Smith, Sebastian; Twumasi, Ricardo; Capon, Charlotte; Pearse, Callum; Smychkovich, Vasilisa; Knight, Alec},
  journal={SSRN},
  year={2025},
  publisher={Elsevier},
  doi={10.2139/ssrn.5382634}
}
</code></pre>
</details>
<details>
<summary>APA</summary>
<pre><code>Kirdar-Smith, S., Twumasi, R. Capon, C., Pearse, C., Smychkovich, V. & Knight, A. (2025). Burnout Prevalence in Paediatric Surgeons: A Systematic Review and Meta-Analysis. SSRN. 10.2139/ssrn.5382634 </code></pre>
</details>
<details>
<summary>Vancouver</summary>
<pre><code>Kirdar-Smith, S, Twumasi, R Capon, C, Pearse, C, Smychkovich, V, Knight, A (2025). Burnout Prevalence in Paediatric Surgeons: A Systematic Review and Meta-Analysis.SSRN. 10.2139/ssrn.5382634 </code></pre>
</details>

---
Contributors: Ricardo Twumasi, Sebastian Kirdar-Smith
