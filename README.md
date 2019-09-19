# readme
A Tool.config file is necessary to run this tool!

This tool is mostly wrote in R3.5 and Python3.7, please make sure it can run on your enviroment.

Before analysis, you can run a script to install dependence by running:
```
Rscript dependence.R
```

The Tool.config file is actually a *json* with **one** object “config", which include two object "transcriptome" and "genotyping". In both transcriptome and genotyping, the first key is **"usage"** which represent if this part is used or not. For example, if the usage in transcriptome is "FALSE" and in genotyping is "TRUE" then the tool will only run the genotyping part, all parameter in transcriptome is ignored.

The details of each parameters are listed below.
## Transcriptome

### Create TPM
If **usage** in transcriptome is TRUE *(or true, T, t)*, you need to select input data type including **fastq, bam, readcounts and TPM matrix** by assigning the corresponding values. Note that if the input is a readcounts file, a gene length file is also needed.

If input files is fastq or bam, one(or two if fastq is paried-end)file(s) only represent one individual, so you have to put all files in one directory and sent the directory to the parameter.

*e.g. "fastq":"/home/Document/experiment/" means all fastq files are in this directory. So does the "bam"*

After determind fastq or bam, the tool will count TPM(Transcripts per million reads) with <a href="https://github.com/deweylab/RSEM" target="_blank">**RSEM**</a> tool, you may need to install it first and confirm the value of each key in object "RSRM", the details of them are listed:

| parameter for RSEM| necessity(Y/N)|  remark  |
| --------   | -----:  | :----:    |
|rsem_path|  Y  |The directory of RSEM|
|pnum|N|number of process for RSEM|
|paired-end|Y|TRUE or FALSE|
|reffile|necessary|The name of the reference used|
|mapping software|necessary|Choose a tool to align reads including bowtie/bowtie2/star|
|bowtie software path|necessary|Directory of mapping software chosen|
|output_genome_bam|optional|alignments in genomic coordinates instead in transcript|
|estimate-rspd|optional|Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD|
|append_names|optional|tells RSEM to append gene_name/transcript_name to the result files|

After RSEM, each individuals will have a .result file as result, which contain a colume "TPM", then the tool will read all the .result file and create a "TPM_matrix.csv" for later analysis.

If input is readcount and gene length file, the tool will call readcounts2TPM.R and create a "TPM_matrix.csv".

If input is a TPM matrix, or we have already create it from other format, we can do the next step.

### transcriptomeQC
All the following process is base on a TPM file.
In this part, we will do quality control and some normalizatation work base on TPM data.

#### Three statistics to identify outliers

##### RLE

RLE(Relative Log Expression)is assume that most expressions in a sample should be near a mean value, only a few genes have differential expression, which means higher or lower than other genes. To calculate it, for each gene *i* , calculate it median in *N* sample as *Medi* , for each sample *j* , and expression value *eij*, count the difference between *eij* and *Medi* : *deij = eij-Medi* , than boxplot each sample base on *deij* and sort by IQR(interquartile range), the sample with lager IQR is more likely to be an outlier. 

##### Hierarchical clustering

We use(1 - spearman's correlation) as distances between samples and assume samples should be homogeneous, so all samples should have short distances between others and clustered homogeneously. So samples far away others may be outliers. So we first use distances which is 1-spearmen correlation to hierarchically clustered, than use top 100 gene sort by variance to calculate Mahalanobis Distance, then a chi2 p-value will be calculated based on mahalanobis distance. Then clusters with $ \geq 60\%$ samples with Bonferroni-corrected p-values < 0.05 were marked as outliers.

##### D-statistic

For each sample, D-statistic represent the average correlation between its expression and other samples. A sample with lower D-statistic is more likely to be an outlier.

After this three steps, we will do quantile normalization, if users have provided coveriance information, we can do gender test if have sex information; do combat to remove batch effect if have batch information. Then use sva to adjust other and hidden coveriance.

| parameter        | necessity(Y/N)   |  remark  |
| --------   | -----:  | :----:  |
|RLEFilterPercent|Y|sample with higher RLE may be outliers|
|DSFilterPercent|Y|D-statistic|
|pvalues.cut|Y|cutoff p-value|
|cluster_level|Y|The max level of hcluster when finding outlier|
|coveriance|Y|TRUE or FALSE|
|coveriance_file|N|coveriance information file, necessary when coveriance is TRUE|

## genotyping

When "usage":"TRUE", the genotyping part is functional. The tools is use plink to process genotype data, if your data is in VCF fromat you can also use plink to change it into plink format. However, only plink 1.9+ can change VCF to plink so make sure you have install plink1.9+ or you may have to change VCF to plink yourself. The tool also use *smartpca* to run PCA analysis, please install it first.

| parameter        | necessity(Y/N)   |  remark  |
| --------   | -----:  | :----:  |
|plink|N|Rawdata in plink format|
|VCF|N|Rawdata in VCF format(plink and VCF should have one and only one value)|
|workdir|N|Set the work directory|
|genotyping_rate|Y|the cutoff of genotyping rate|
|call_rate|Y|the cutoff of call rate|
|HWE|Y|HWE Test's cutoff p-value|
|MAF|Y|Minor Allele Frequency|
|F_outlier_n_sd|Y|determind F outlier with fmean $\pm$F_outlier_n_sd\*fsd|
|population|N|Population stratification| 

#### PCA
| parameter for smartpca| remark  |
| ------ | :----:  |
|m|default(5)maximum number of outlier removal iterations.To turn off outlier removal, set -m 0.|
|k|default(10)number of principal components to output|
|t|default(10)number of principal components along which to remove outliers during each outlier removal iteration|

#### Imputation
If the usage of Imputation is TRUE, the tool will use <a href="https://www.well.ox.ac.uk/~wrayner/tools/" target="_blank">**McCarthy's tool**</a>McCarthy's tool to do prepreparing checking. The tool can use three kinds of reference panel: HRC, 1000G and CAAPA, use the parameter in config file to seletc and determind the reference panel.

