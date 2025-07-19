# Data used in the preprint/publication

## FASTQ and BLOW5

The ENA project [PRJEB86853](https://www.ebi.ac.uk/ena/browser/view/PRJEB86853) contains raw data with file name prefixes as in this tsv file [here](data-fastq.tsv). In the tsv table:

- For PacBio HiFi data there should be one fastq.gz file.
- For ONT data there should ba a fastq.gz containing basecalled reads and a tar.gz containing raw signal BLOW5 files.											
- Note that FASTQ files contain what was given to the assembler, that are already filtered (e.g., pass reads, qscore cut off).
- FASTQ for ONT duplex samples only contain duplex reads.
- For cornetto rounds for ONT simplex, FASTQ files contain reads >30kbase.											
- BLOW5 files contain everything that came out of the sequencing run without any filtering.
- Those with T2TC inside brackets are where we used publicly available data from T2T Consortium's human pangenomics AWS bucket.	The links we used are:
  - [m84005_220827_014912_s1 (T2TC)](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84005_220827_014912_s1.hifi_reads.fastq.gz)
  - [HG002_3 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_3_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz), [HG002_6 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_6_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz), [HG002_2 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/1_3_23_R1041_Duplex_HG002_2_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz), [HG002_9 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_9_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz), [HG002_4 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_4_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz), [HG002_8 (T2TC)](https://human-pangenomics.s3.amazonaws.com/submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex/11_15_22_R1041_Duplex_HG002_8_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz)      				
- *hg002-Cornetto-x.1 is made of fastq0+fastq1; hg002-Cornetto-x.2 is made of fastq0+fastq1+fastq2; and so on.											


## FASTA assembly

### HG002

The files in the table [here](data-asm-hg002.tsv) are found inside the cornetto-hg002-asm directory created when you download and extract the cornetto-hg002-asm.tar.gz file from [https://doi.org/10.5061/dryad.kkwh70sfr](https://doi.org/10.5061/dryad.kkwh70sfr)

### Animals

The files in the table [here](data-asm-animals.tsv) are found inside the cornetto-animal-asm directory created when you download and extract the cornetto-animal-asm.tar.gz file from [https://doi.org/10.5061/dryad.kkwh70sfr](https://doi.org/10.5061/dryad.kkwh70sfr)
