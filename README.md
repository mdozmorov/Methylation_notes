# DNA methylation notes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

Methylation-related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Pipelines](#pipelines)
  - [Bisulfite sequening](#bisulfite-sequening)
- [Preprocessing](#preprocessing)
- [Differential methylation](#differential-methylation)
- [Deconvolution](#deconvolution)
- [CNV](#cnv)
- [Integrative](#integrative)
- [Data](#data)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Pipelines

- [**RnBeads 2.0**](https://rnbeads.org/) - Full analysis pipeline for Illumina 27K/450K/EPIC arrays, WGBS/RRBS/target bisulfite seq. Figure 1 - all steps, from QC, preprocessing, exploratory analysis, to differential methylation/region analysis. New analyses include sex prediction, genetic purity, CNV analysis, cell type inference, age prediction, and more. 
    - Müller, Fabian, Michael Scherer, Yassen Assenov, Pavlo Lutsik, Jörn Walter, Thomas Lengauer, and Christoph Bock. “[**RnBeads 2.0: Comprehensive Analysis of DNA Methylation Data.**](https://doi.org/10.1186/s13059-019-1664-9)” Genome Biology 20, no. 1 (December 2019): 55. 

- [**MEAL**](https://www.bioconductor.org/packages/release/bioc/html/MEAL.html) - R/Bioconductor Package to integrate methylation and expression data. It can also perform methylation or expression analysis alone. Wraps other packages

- [meffil](https://github.com/perishky/meffil) - Efficient algorithms for analyzing DNA methylation data. R implementation.

- [**GLINT**](https://github.com/cozygene/glint) - methylation array data analysis pipeline. Illumina 27K/450K/EPIC. ReFACTor to adjust for tissue heterogeneity, inferring population structure, imputation, association testing, basic visualization. Python implementation. 
    - Rahmani, Elior, Reut Yedidim, Liat Shenhav, Regev Schweiger, Omer Weissbrod, Noah Zaitlen, and Eran Halperin. “[**GLINT: A User-Friendly Toolset for the Analysis of High-Throughput DNA-Methylation Array Data.**](https://doi.org/10.1093/bioinformatics/btx059)” Edited by John M Hancock. Bioinformatics 33, no. 12 (June 15, 2017): 1870–72.

### Bisulfite sequening

- Benchmarking of whole-genome bisulfite converted sequencing data (WGS methylation). Nine tools, tested on simulated data. BWA-mem and BSMAP perform best, Bismark and others are close, often faster and less memory demanding
    - Nunn, Adam, Cristian Otto, Peter F Stadler, and David Langenberger. “[Comprehensive Benchmarking of Software for Mapping Whole Genome Bisulfite Data: From Read Alignment to DNA Methylation Analysis](https://doi.org/10.1101/2020.08.28.271585).” BioRxiv, 2020.

- [Meth10X](https://github.com/luuloi/Meth10X) - Whole Genome Bisulfite Sequencing Alignment pipeline. Includes SNP calling, region calling (partial methylation, low, high etc.), basic downstream analyses. Designed for SGE/PBS clusters. Built based on in-house bash/python/perl/R script, Bpipe and common alignment tools.

- [bwa-meth](https://github.com/brentp/bwa-meth) - fast and accurate alignment of BS-Seq reads https://arxiv.org/abs/1401.1129. 

- [GemBS](http://statgen.cnag.cat/gemBS/) - alignment to converted and regular reference genome, calling methylated CpGs. 

- [MethylDackel](https://github.com/dpryan79/MethylDackel) - A (mostly) universal methylation extractor for BS-seq experiments. Will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq alignments and extract per-base methylation metrics from them.

- [methylpy](https://github.com/yupenghe/methylpy) - processing bisulfite sequencing and NOME-seq data, supports single- and paired-end data, read trimming, QC, duplicate removal, outputs methylation states and/or open chromatin, differentially methylated region calling at single cytosine level, multi-group comparison. 
    - Schultz, Matthew D., Yupeng He, John W. Whitaker, Manoj Hariharan, Eran A. Mukamel, Danny Leung, Nisha Rajagopal, et al. “[**Human Body Epigenome Maps Reveal Noncanonical DNA Methylation Variation.**](https://doi.org/10.1038/nature14465)” Nature 523, no. 7559 (July 9, 2015): 212–16.

- [MOABS](https://code.google.com/archive/p/moabs/) - bisulfite sequencing data processing and differential analysis pipeline. Beta-binomial hierarchical model, priors estimated using an Empirical Bayes approach. Credible methylation differences (CDIFs) better capture the biology of methylation differences. Outperforms BSmooth. Can account for CNVs. C++ implementation 
    - Sun, Deqiang, Yuanxin Xi, Benjamin Rodriguez, Hyun Jung Park, Pan Tong, Mira Meong, Margaret A Goodell, and Wei Li. “[MOABS: Model Based Analysis of Bisulfite Sequencing Data](https://doi.org/10.1186/gb-2014-15-2-r38),” 2014, 12.

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

- [deconvR](https://github.com/BIMSBbioinfo/deconvR) - an R package designed for analyzing deconvolution of the bulk sample(s) using an atlas of reference omic signature profiles and a user-selected model. Works with RNA-seq, array- and whole genome DNA methylation.

- https://github.com/stephaniehicks/methylCC - R/BioC package to estimate the cell composition of whole blood in DNA methylation samples in microarray or sequencing platforms


## CNV

- `conumee` - copy-number variation analysis using Illumina DNA methylation arrays


## Integrative

- `ELMER` - Inferring Regulatory Element Landscapes and Transcription Factor Networks Using Cancer Methylomes. https://bioconductor.org/packages/release/bioc/html/ELMER.html

- `omicade4` - Multiple co-inertia analysis of omics datasets. https://bioconductor.org/packages/release/bioc/html/omicade4.html

- `mixOmics` - mixOmics R Package. merging datasets by different measures on samples, or on genes. Various unsupervised and supervised analyses. Regularization, lasso for reature selection. http://mixomics.org/

- `MLExpResso` - Package for analyzing genes expression and CpG probes metylation. https://geneticsmining.github.io/MLGenSig/index.html

## Data

- DNA methylation atlas of normal human cell types (39 sorted cell types). Whole-genome bisulfite sequencing, fragment-level analysis. 99.5% reproducibility within the same cell type. Reconstructs lineages. Cell type-specific markers as putative enhancers, associated with various epigenetic marks. A role for DNA methylation to attenuate CTCF binding and modulate cell-specific 3D structure. hg19, hg38 genome assembly, BED, BigWig files at [GSE186458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458). Supplementary tables with most variable methylation blocks, cell type-specific unmethylated markers, markers enriched for functional terms, gene associations, cell type-specific markers on methylation arrays. [wgbs_tools](https://github.com/nloyfer/wgbs_tools) for alignment/visualization, [UXM_deconv](https://github.com/nloyfer/UXM_deconv) for deconvolution. <details>
    <summary>Paper</summary>
    Loyfer, N., Magenheim, J., Peretz, A. et al. A DNA methylation atlas of normal human cell types. Nature 613, 355–364 (2023). https://doi.org/10.1038/s41586-022-05580-6
</details>

- [DNA methylation atlas](http://neomorph.salk.edu/omb/) of the mouse brain, single-nucleus sequencing (snmC-seq2, [pipeline](https://hq-1.gitbook.io/mc/), [ALLCools: ALL methyl-Cytosine tools](https://github.com/lhqing/ALLCools)). 103,982 nuclei (including 95,815 neurons and 8,167 non-neuronal cells) from 45 regions of the mouse cortex, hippocampus, striatum, pallidum and olfactory areas. 161 cell clusters with distinct spatial locations and projection targets. Integration with scATAC-seq confirms cell types and identifies enhancer-gene interactions. HDF5 data on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132489), and links to individual experiments and data in TSV format. Data on [Nemoarchive](https://assets.nemoarchive.org/dat-vmivr5x). [Merged methylome viewer](http://neomorph.salk.edu/mouse_brain.php). [Brain Cell Methylation Viewer](http://neomorph.salk.edu/omb/) to view gene, brain region, cell type. <details>
    <summary>Paper</summary>
    Liu, Hanqing, Jingtian Zhou, Wei Tian, Chongyuan Luo, Anna Bartlett, Andrew Aldridge, Jacinta Lucero, et al. “DNA Methylation Atlas of the Mouse Brain at Single-Cell Resolution.” Nature 598, no. 7879 (October 7, 2021): 120–28. https://doi.org/10.1038/s41586-020-03182-8.
</details>

- [EWASdb](http://www.bioapp.org/ewasdb/index.php/Index/index) - EWASdb is the first epigenome-wide association database (first online at 2015, and first published on Nucleic Acids Res. 2018 Oct 13) which stores the results of 1319 EWAS studies associated with 302 diseases/phenotypes (p<1e-7). Three types of EWAS results were stored in EWASdb: EWAS for single epi-marker; EWAS for KEGG pathway and EWAS for GO (Gene Ontology) categories. http://www.bioapp.org/ewasdb/

- [EWAS Atlas](https://bigd.big.ac.cn/ewas/index) - A knowledgebase of epigenome-wide association studies. https://bigd.big.ac.cn/ewas/index

- `CoRSIVs` - 9926 correlated regions of systemin interindividual variation of DNA methylation. GTEx data. Enriched in subtelomeric regions, transposable elements, depleted in TFBSs. Enriched in Quiescent regions, repressive polycomb marks, depleted in heterochromatin, active promoters and enhancers (bivalent). Likely genetically driven. [Supplementary matierial](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1708-1#Sec13): table S2 - all significant 39,424, S3 - filtered 9,926, S13 - 1659 450K probes overlapping 819 CoRSIVs.
    - Gunasekara, Chathura J., C. Anthony Scott, Eleonora Laritsky, Maria S. Baker, Harry MacKay, Jack D. Duryea, Noah J. Kessler, et al. “A Genomic Atlas of Systemic Interindividual Epigenetic Variation in Humans.” Genome Biology 20, no. 1 (December 2019): 105. https://doi.org/10.1186/s13059-019-1708-1.

- The 450K array measures the methylation status of 485,512 methylcytosine sites in the human genome at a single nucleotide resolution, representing approximately 1.5% of total genomic CpG sites [22126295, 21593595]. While the assayed CpG sites are concentrated around promoter regions and gene bodies, approximately 25% are located in intergenic regions [22126295].

## Misc

- Methylation inhibits TFBSs, but some factors, like homeodomain, POU, NFAT, prefer binding to methylated DNA. These TFs play a role in embryonic and organismal development. Yin, Yimeng, Ekaterina Morgunova, Arttu Jolma, Eevi Kaasinen, Biswajyoti Sahu, Syed Khund-Sayeed, Pratyush K. Das, et al. “Impact of Cytosine Methylation on DNA Binding Specificities of Human Transcription Factors.” Science (New York, N.Y.) 356, no. 6337 (May 5, 2017). doi:10.1126/science.aaj2239.

