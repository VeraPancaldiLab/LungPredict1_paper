# LungPredict 

Profiling the Non-Small Cell Lung Cancer (NSCLC) microenvironment by integrating transcriptomics to uncover potential  phenotypic profiles associated to patterns in immune infiltration

## Summary
Here, we applied a computational immunology approach involving differential expression and pathway analyses, along with the quantification of immune cell proportions using computational deconvolution methods, inference of transcription factor activities, and specification of immune score types to better characterize bulk transcriptomics samples of lung adenocarcinoma. This comprehensive analysis enabled us to identify biomarkers of disease progression and immune infiltration patterns across disease stages and patient groups, thus allowing us to compute and evaluate specific hallmarks of immune response across patients.Through our methodology and feature integration pipeline, we identified groups of immune cells related to disease stage as well as potential immune response or evasion. More specifically, we reported a duality in the behavior of immune cells, notably natural killer (NK) cells, which could be relevant for prognosis and potentially, immune response or evasion. The dual profile of other immune cells, most notably T-cell populations, have been discussed in the context of disease such as cancer. Here, we report the duality of NK cells which should be taken into account in conjunction with other immune cell populations and behaviors in predicting prognosis, immune response or evasion. The potentially predictive value of these populations in determining prognosis at different stages of the disease will be further investigated in other cohorts. 

## Data used in the analysis
- LungPredict bulk RNAseq data: 62 patients - Adenocarcinoma - Stage I, II, III, IV
- Vanderbilt bulk RNAseq data: 70 patients - Adenocarcinoma - Early Stage (I, II)
- Vanderbilt single-cell RNAseq data: 15 patients - Adenocarcinoma - Early Stage (I,II)

![image](https://github.com/VeraPancaldiLab/LungPredict1/assets/37853385/2641fa06-91e4-46f5-bc6f-4f83baacb035)

## Project organization
- **Scripts**: Codes used for analysis.
  -  `Data_preprocessing.Rmd`:
  -  `LP_all_analysis.Rmd`: Analysis for LungPredict all stages samples. 
  -  `LP_early_analysis.Rmd`: Analysis for LungPredict early stages samples.
  -  `Vanderbilt_analysis.Rmd`: Analysis for Vanderbilt early stages samples.
  -  `Vanderbilt_scRNAseq_analysis.Rmd`: scRNAseq analysis of Vanderbilt cohort.
  -  `Vanderbilt_survival_analysis.Rmd`: Survival analysis for Vanderbilt early stage samples.
  -  `TCGA_data_preprocessing.Rmd`: Data retrieval and processing from TCGA.
  -  `TCGA_analysis.Rmd`: Analysis using TCGA data.
- **Figures**: Figures used in the paper. 
- **Results**: Results files from analysis.

## Computational methods
- Differential expression analysis.
- Immune cell type deconvolution.
- Transcription factor (TF) activity inference.
- Pathway activity calculation.
- Hallmarks of immune response.
- Feature selection algorithms.
- Enrichment analysis.
- scRNAseq analysis.
- Reference-based deconvolution.
- Survival analysis.
- TCGA analysis.

## Specifications
Analysis was done using R version 4.3.1 with the OS Ubuntu 22.04.3 LTS.

## Reproducibility
If you would like to reproduce the analysis done here, we invite you to use our provided r-environment. Setting it up will install all the neccessary packages, along with their specific versions in an isolated environment.  

For this, open the project `LungPredict.Rproj` inside the `scripts/` folder and in the R console run:

```r
# Download renv package (if not installed)
install.packages('renv')
# To activate the R environment 
renv::activate()
# To download and install all the require libraries and packages 
renv::restore() 
```

Note that this is an **once-step** only when running LungPredict for the **first time**. For the following times, you will only need to open the `LungPredict.Rproj` and you are ready to go!

Once all packages have been installed, you can start reproducing the analysis using the scripts inside the `scripts/` folder.

Make sure to run `renv::deactivate()` when finishing, to avoid conflicts whenever you start a different R project.

For more information about how R-environments work, visit the main page of the tool [renv](https://rstudio.github.io/renv/articles/renv.html).

## Contributing
If you are interested or have questions about the analysis done in this project, we invite you to open an issue in https://github.com/VeraPancaldiLab/LungPredict1_paper/issues or contact Marcelo Hurtado (marcelo.hurtado@inserm.fr) for more information.

## Authors
- [Marcelo Hurtado](https://github.com/mhurtado13)
- [Leila Khajavi](https://github.com/LeilaKhajavi)
- [Mouneem Essabar](https://github.com/mouneem)
- [Vera Pancaldi](https://github.com/VeraPancaldi)
