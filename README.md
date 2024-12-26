# Estimating the Probability of CPVT from RYR2 Genetic Variant Properties

Here we propose a method to estimate the post-test probability of CPVT for variants in the cardiac calcium channel **RYR2**. All code and data used and referenced in the manuscript (not yet accepted) are included here.

An overview of the method and results is provided in **literature_process_RYR2.R**.  
The raw data are contained in **RYR2_20241129.xlsx**, which we curated from literature and gnomAD. The first half of **literature_process_RYR2.R** and **literature_process_RYR2_test.R** processes these data to produce **d_full.csv**.  

We then added three pathogenic scores (AlphaMissense, REVEL, and ClinVar) to each variant via ChatGPT, producing two additional files:  
- **d_full_Model1_AM_REVEL_ClinVar.csv**  
- **d_full_AM_REVEL_ClinVar_againstALL.csv**  

The raw distances between centroid atoms in the three-dimensional structure of Ryanodine receptor 2 (RYR2) can be found in the **RYR2_AF2_distances** folder. These distance data are used to generate the CPVT probability density via the **func_dist_seq.R** script.

---

## Folders

- **RYR2_AF2_distances/**  
  Contains two datasets that must be combined by the script. (They are split due to file size constraints on GitHub.)  

- **pathogenic_scores/**  
  Contains all pathogenic score data, including AlphaMissense and REVEL (used as predictive covariates). ClinVar data are used only for model evaluation.  

- **variant_browser/**  
  Holds the raw data used for the variant browser website.

---

Feel free to open an issue or contact us if you have questions about this dataset or the scripts. We will update this repository as our manuscript progresses through review.
