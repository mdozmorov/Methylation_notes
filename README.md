# Notes on DNA methylation analysis

Primarily array-oriented, human/mouse. These notes are not intended to be comprehensive. They include notes about methods, packages and tools I would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/Bioinformatics_notes) and [collections of links to various resources](https://github.com/mdozmorov/MDmisc_notes). Issues with suggestions and pull requests are welcome!

* [Pipelines](#pipelines)
* [Alignment](#alignment)
* [Preprocessing](#preprocessing)
* [Differential methylation](#differential-methylation)
* [CNV](#cnv)
* [Integrative](#integrative)
* [Misc](#misc)


## Pipelines

- `scRRBS` - single-cell RRBS notes and scripts by Ming Tang, https://gitlab.com/tangming2005/scRRBS, start with https://gitlab.com/tangming2005/scRRBS/blob/master/quality_control.md

- `MEAL` - Package to integrate methylation and expression data. It can also perform methylation or expression analysis alone. Wraps other packages

- `GLINT` - methylation array data analysis pipeline. Illumina 27K/450K/EPIC. ReFACTor to adjust for tissue heterogeneity, inferring population structure, imputation, association testing, basic visualization. Python implementation. https://github.com/cozygene/glint
    - Rahmani, Elior, Reut Yedidim, Liat Shenhav, Regev Schweiger, Omer Weissbrod, Noah Zaitlen, and Eran Halperin. “GLINT: A User-Friendly Toolset for the Analysis of High-Throughput DNA-Methylation Array Data.” Edited by John M Hancock. Bioinformatics 33, no. 12 (June 15, 2017): 1870–72. https://doi.org/10.1093/bioinformatics/btx059.


## Alignment

- `bwa-meth` - fast and accurate alignment of BS-Seq reads https://arxiv.org/abs/1401.1129. https://github.com/brentp/bwa-meth
- `GemBS` - alignment to converted and regular reference genome, calling methylated CpGs. http://statgen.cnag.cat/gemBS/


## Preprocessing

- `ENmix` - data preprocessing, batch, PCA. https://www.bioconductor.org/packages/release/bioc/html/ENmix.html


## Differential methylation

- Review of differential methylation methods and 22 tools. Categorized by approaches. [Pros and cons of each approach](https://academic.oup.com/view-large/figure/121909042/bbx013f3.tif), [Table 1. Summary of the important characteristics of the 22 surveyed approaches](https://academic.oup.com/view-large/121909044), [Table 2. Comparison of the available implementations of the 22 surveyed approaches](https://academic.oup.com/view-large/121909057)
    - Shafi, Adib, Cristina Mitrea, Tin Nguyen, and Sorin Draghici. “A Survey of the Approaches for Identifying Differential Methylation Using Bisulfite Sequencing Data.” Briefings in Bioinformatics, March 8, 2017. https://doi.org/10.1093/bib/bbx013.

- `DMRcate` - https://bioconductor.org/packages/release/bioc/html/DMRcate.html
- `metilene` - https://www.bioinf.uni-leipzig.de/Software/metilene/
MOABS, BSmooth, 
- `BiSeq` - https://bioconductor.org/packages/release/bioc/html/BiSeq.html
- `RADMeth` - Regression Analysis of Differential Methilation is a software for computing individual differentially methylated sites and genomic regions in data from whole genome bisulfite sequencing (WGBS) experiments. http://smithlabresearch.org/software/radmeth/. A part of `MethPipe` - a computational pipeline for analyzing bisulfite sequencing data. http://smithlabresearch.org/software/methpipe/
- `DSS` - Differential methylation analysis for general experimental design, based on a beta-binomial generalized linear model with arcsine link function. https://bioconductor.org/packages/release/bioc/html/DSS.html
- `DMRfinder` - Following Bismark, extracts CpG methylation, cluster CpGs into regions, tests for differential methylation using DSS package (Bayesian beta-binomial hierarchical modeling). https://github.com/jsh58/DMRfinder

- `dmrff` - differentially methylated regions based on inverse-variance weighted meta-analysis. Outperforms bumphunter, Comb-p, DMRcate, seqlm, has more power. https://github.com/perishky/dmrff
    - Suderman, Matthew, James R Staley, Robert French, Ryan Arathimos, Andrew Simpkin, and Kate Tilling. “Dmrff: Identifying Differentially Methylated Regions Efficiently with Power and Control.” BioRxiv, January 1, 2018, 508556. https://doi.org/10.1101/508556.


## Deconvolution

- https://github.com/stephaniehicks/methylCC - R/BioC package to estimate the cell composition of whole blood in DNA methylation samples in microarray or sequencing platforms


## CNV

- `conumee` - copy-number variation analysis using Illumina DNA methylation arrays


## Integrative

- `ELMER` - Inferring Regulatory Element Landscapes and Transcription Factor Networks Using Cancer Methylomes. https://bioconductor.org/packages/release/bioc/html/ELMER.html

- `omicade4` - Multiple co-inertia analysis of omics datasets. https://bioconductor.org/packages/release/bioc/html/omicade4.html

- `mixOmics` - mixOmics R Package. merging datasets by different measures on samples, or on genes. Various unsupervised and supervised analyses. Regularization, lasso for reature selection. http://mixomics.org/

- `MLExpResso` - Package for analyzing genes expression and CpG probes metylation. https://geneticsmining.github.io/MLGenSig/index.html


## Misc

- `CoRSIVs` - 9926 correlated regions of systemin interindividual variation of DNA methylation. GTEx data. Enriched in subtelomeric regions, transposable elements, depleted in TFBSs. Enriched in Quiescent regions, repressive polycomb marks, depleted in heterochromatin, active promoters and enhancers (bivalent). Likely genetically driven. [Supplementary matierial](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1708-1#Sec13): table S2 - all significant 39,424, S3 - filtered 9,926, S13 - 1659 450K probes overlapping 819 CoRSIVs.
    - Gunasekara, Chathura J., C. Anthony Scott, Eleonora Laritsky, Maria S. Baker, Harry MacKay, Jack D. Duryea, Noah J. Kessler, et al. “A Genomic Atlas of Systemic Interindividual Epigenetic Variation in Humans.” Genome Biology 20, no. 1 (December 2019): 105. https://doi.org/10.1186/s13059-019-1708-1.

- Methylation inhibits TFBSs, but some factors, like homeodomain, POU, NFAT, prefer binding to methylated DNA. These TFs play a role in embryonic and organismal development. Yin, Yimeng, Ekaterina Morgunova, Arttu Jolma, Eevi Kaasinen, Biswajyoti Sahu, Syed Khund-Sayeed, Pratyush K. Das, et al. “Impact of Cytosine Methylation on DNA Binding Specificities of Human Transcription Factors.” Science (New York, N.Y.) 356, no. 6337 (May 5, 2017). doi:10.1126/science.aaj2239.

- The 450K array measures the methylation status of 485,512 methylcytosine sites in the human genome at a single nucleotide resolution, representing approximately 1.5% of total genomic CpG sites [22126295, 21593595]. While the assayed CpG sites are concentrated around promoter regions and gene bodies, approximately 25% are located in intergenic regions [22126295].
