
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intratumormeth

The goal of this package is to facilitate analysis of DNA methylation
data obtained from Illumina MethylationEPIC platform for multiple
samples of the same tumor.

## Installation

You can install the development version of intratumormeth from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fynnwi/intratumormeth")
```

## Workflow

Every stage of the workflow is described in a separate vignette:

1.  `preprocessing`: read IDAT files, extract methylation values
2.  `quality_control`: investigate detection p-values, tumor purity,
    find poor quality samples and probes
3.  `copynumber`: run conumee and generate CNV data
4.  `classification`: evaluate glioma subtype and MGMT promoter
    methylation
5.  `dimensionality_reduction`: perform PCA and t-SNE to find clusters
6.  `phylogeny`: run TuMult separate CNV data

<!-- ## To be implemented: -->
<!-- Batch correction -->
<!-- - require user input of technical variables -->
<!--   - bisulfite conversion batch -->
<!--   - chip, row, column id -->
<!-- - run PCA and t-SNE and color by potential confounder variables -->
<!-- Probe filtering -->
<!-- - implement option to apply various filters to methylation data: -->
<!--   - sex chromosomes -->
<!-- Dimensionality reduction -->
<!-- - implement tsne function -->
<!-- Probe-wise heterogeneity analysis -->
<!-- - calculate intra-patient standard devation of beta values for each probe -->
<!-- - average these sds across all patients -->
<!-- - find probes that are commonly "methylation-instable"? -->
<!--   - what is their genomic context? distribution across chromosomes, regulatory features etc. -->
<!-- Gene-wise CNV analysis -->
<!-- - calculate intra-tumor mean and standard deviation for each gene (i.e. conumee detail region) -->
<!-- Sample table -->
<!-- | sample_id | patient_id | idat_location                                       | project | ceccarelli_subtype | capper_subtype  | mgmt_stp27 | -->
<!-- |-----------|------------|-----------------------------------------------------|---------|--------------------|-----------------|------------| -->
<!-- | sfb01_1   | sfb01      | ~/data/methylation_data/Pat01.1/202212330261_R07C01 | SFB824  | Classic-like       | RTK1            | M          | -->
<!-- | sfb01_2   | sfb01      | ~/data/methylation_data/Pat01.2/202242420132_R07C01 | SFB824  | PA-like            | GBM Mesenchymal | U          | -->
<!-- | ...       | ...        | ...                                                 | ...     | ...                | ...             | ...        | -->
<!-- Patient table -->
<!-- | patient_id | sex | age | n_samples | beta_mean | beta_sd | ... | -->
<!-- |------------|-----|-----|-----------|-----------|---------|-----| -->
<!-- | sfb01      | m   | 57  | 4         | 0.55      | 0.300   | ... | -->
<!-- | sfb02      | m   | 76  | 3         | 0.53      | 0.222   | ... | -->
<!-- | ...        | ... | ... | ...       | ...       | ...     | ... | -->
<!-- Detail region table -->
<!-- Collect cohort-level information about detail regions: -->
<!-- | detail_region | chr | start    | end      | n_probes | cnv_mean | cnv_sd | beta_mean | beta_sd | surroundings_cnv_mean | ... | -->
<!-- |---------------|-----|----------|----------|----------|----------|--------|-----------|---------|-----------------------|-----| -->
<!-- | PTEN          | 10  | 89621740 | 89725510 | 20       | -0.351   | 0.032  | 0.388     | 0.1     | 0.0234                | ... | -->
<!-- | PDGFRA        | ... | ...      | ...      | ...      | ...      | ...    | ...       | ...     | ...                   | ... | -->
<!-- | CDK2A/B       | ... | ...      | ...      | ...      | ...      | ...    | ...       | ...     | ...                   | ... | -->
<!-- | ...           | ... | ...      | ...      | ...      | ...      | ...    | ...       | ...     | ...                   | ... | -->
<!-- CpG table -->
<!-- Idea: summarize all cohort-level information related to single CpG sites in a table like this: -->
<!-- | probe_id     | mean | sd  | min | max | mean_intratumor_sd | sd_intratumor_mean | sd_intratumor_sd | neighborhood_correlation | -->
<!-- | ------------ | ---- | --- | --- | --- | ------------------ | ------------------ | ---------------- | ---- | -->
<!-- | cg12412432   | 1    |  1  | 0.1 | 0.9 | 0.002              | 0.75               | 0.0013           | ...  | -->
<!-- | cg12412432   | 1    |  1  | 0.1 | 0.9 | 0.002              | 0.75               | 0.0013           | ...  | -->
<!-- | ...          | ...  | ... | ... | ... | ...                | ...                | ...              | ...  | -->
<!-- ## Tumor purity -->
<!-- ### InfiniumPurify -->
<!-- - [paper](https://www.sciencedirect.com/science/article/pii/S2352304218300163) from 2018, 26 citations -->
<!-- - Wenger et al. use this tool -->
<!--   - using 51 adult brain tissue samples from GSE36278 (Sturm 2012, six normal brain samples), GSE40360 (Casaccia 2013, 19 control samples), GSE50798 (Dracheva 2013, postmortem 24 samples) and GSE52556 (Kleinman 2013, 34 control samples) as the normal data set -->
<!-- - normal tissue not necessarily required -->
<!--   - from vignette: "If either the number of tumor samples or number of normal smaples is less than 20, the tumor.type argument should be specified according to CancerTypeAbbr. If the numbers of tumor and normal samples are both more than 20, tumor.type could be null." -->
<!-- ### PAMES -->
<!-- - Purity Assesment from DNA MEthylation Sites -->
<!-- - [publication Benelli et al.](https://academic.oup.com/bioinformatics/article/34/10/1642/4792963?login=true) 2018 -->
<!--   - contains comparison to InfiniumPurify, they find high correlation -->
<!-- - [GitHub repo] (https://github.com/cgplab/PAMES) -->
<!-- - used by Verburg et al.  -->
<!--   - validation of purity estimates by looking at correlation of purity with chr7 amplification and chr10 deletion -->
<!-- - methylation data from normal tissue required! -->
<!-- Workflow: -->
<!-- 1. find CpG sites that are predictive for tumor purity by comparing tumor methylation against control group methylation data - these "informative CpG sites" are cancer type specific -->
<!-- 2. derive tumor purity (ratio of tumor to normal cells) from these informative sites -->
<!-- ## Continuous phylogeny tools -->
<!-- ### [BEAST](https://beast.community/index.html) -->
<!-- - there is also a [book from the author](https://www.amazon.com/dp/1107019656/) -->
<!-- - there exists [BEAST 2](http://www.beast2.org/) -->
<!--   - University of Auckland -->
<!-- ### [RevBayes](https://revbayes.github.io/) -->
<!-- - [publication](https://academic.oup.com/sysbio/article/65/4/726/1753608) with 360 citations -->
<!-- - developed by Sebastian Höhna (LMU Paleontology) -->
<!-- ### [Mesquite](http://www.mesquiteproject.org/) -->
<!-- - "modular system for evolutionary analysis" -->
<!-- - updated Aug 2021 -->
<!-- - support for continuous characters -->
<!-- ### [BayesPhylogenies](http://www.evolution.rdg.ac.uk/BayesPhy.html) -->
<!-- - University of Reading, UK -->
<!-- - from 2004 -->
<!-- ## Neuro-Oncology feedback -->
<!-- ### Reviewer 1 -->
<!-- - Ceccarelli classifier only discrete bins -->
<!-- - ~~correction for tumor purity~~ -->
<!-- - conclusions on copy number variations of individual genes are over interpreted -->
<!--   - many genes are on same chromosome -> LOH -->
<!-- - which CpGs make up most of the differences between biopsies?  -->
<!--   - Friedman test -->
<!-- - samples snap frozen or FFPE or mixed? -->
<!-- - evolutionary trajectories for all samples in supplementary information -->
<!-- ### Reviewer 2 -->
<!-- - clinical implications on heterogeneity? -->
<!-- - ~~tumor purity assessment based on methylation profiles~~ -->
<!-- - insufficient figure legends -->
<!-- - data availability missing -->
<!-- ### Reviewer 3 -->
<!-- - methylation heterogeneity missing -->
<!-- - ~~tumor purity~~ -->
<!--   - purity distribution as assessed by pathologist -->
<!-- - probe filtering (SNPs, sex chr etc.) -->
<!-- - settings for CNV analysis -->
<!-- - better description of TuMult -->
<!-- - classification probability scores -->
<!-- - percentage class-mes, mes-pa-like -->
<!-- - clustering:  -->
<!--   - color by Ceccarelli classification -->
<!--   - description of how patients clustered in relation to each other -->
<!-- - drop IDH-mut and use 914 CpGs described by Ceccarelli et al. to visualize samples and patients -->
<!-- - ~~tumor purity vs MGMT heterogeneity~~ -->
<!-- - compare MGMT heterogeneous and subtype heterogeneous patients -->
<!-- - Fig 4A y-axis label -->
<!-- - neutral vs. gain/loss heterogeneity clinically more interesting than medium vs. high gain/loss -->
<!-- - percentages how many samples (and patients) had specific CNV event -->
<!-- - phylogenetic distances between samples -->
<!-- - methylation vs. CNV heterogeneity -->
<!--   - do subtypes have their own branches in phylogeny? -->
<!-- - phylotrees of other patients in supplementary materials -->
<!-- - figure legends too brief -->
<!-- - compare CNV event frequencies with literature (Sottoriva, Taylor) -->
