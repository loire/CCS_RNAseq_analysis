# CCS_RNAseq_analysis

Dataset and analysis of *Aedes aegypti* expression profiles following infection with dengue or RVF viruses

## Reproductibility

Analysis and figures can be generated using the following command on an unix terminal

```
git clone https://github.com/loire/CCS_RNAseq_analysis
cd CCS_RNAseq_analysis
Open analyse.rmd with Rstudio to knit the html or pdf report and generate tables and figures  
```

## Dependencies

R packages: [tidyverse](https://www.tidyverse.org/), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html),[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) 


## Material and methods

* 30 samples
	* 12 biological samples (3 replicates each)
		* RVF 24H 
		* RVF 6 days 
		* Dengue 24H 
		* Dengue 6 days
	* 18 controls (mock community; 3 replicates);
		* "water" (mockA),  inactivated viruses (mockB), stress inducing media mockC)) 
		* 24h and 6 days  

Samples are described in file Data/sample.csv

Library construction and sequencing have been performed by Montpellier Genomix

## Bioinformatic pipeline

MGX contribution
* fastq files were mapped with bowtie2 on *Aedes aegypti* [genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000004015.4/)  and transcripts raw abundance inferred with Tophat2 pipeline.

Raw transcripts abundance are reported in in file Data/Raw_Counts_RNA-Seq_CetreSossah.txt



![MDS plot of filtered dataset](Figures/MDS_GOOD_DATA.pdf)



