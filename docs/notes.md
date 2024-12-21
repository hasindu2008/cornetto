# Simplex manual

slow5tools merge slow5/ -o D_1_PGXXXX240596.blow5


/usr/bin/time -v /install/slow5-dorado-0.8.3/bin/slow5-dorado basecaller -x cuda:all /install/dorado-0.8.3/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0/ D_1_PGXXXX240596.blow5 --emit-fastq --min-qscore 10  > D_1_PGXXXX240596.fastq

/install/seqkit-v2.3.0/seqkit stats --threads 16 -o D_1_PGXXXX240596.seqkit_stats.tsv --all --tabular ${IN} D_1_PGXXXX240596.fastq

/install/seqkit-v2.3.0/seqkit seq -m 2000 D_1_PGXXXX240596.fastq -o D_1_PGXXXX240596.filt.fastq