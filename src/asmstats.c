/**
MIT License

Copyright (c) 2019 Hasindu Gamaarachchi (hasindu@unsw.edu.au)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/

#define _XOPEN_SOURCE 700
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include "error.h"
#include "pafrec.h"
#include "khash.h"
#include "ksort.h"

//just to get the statistics, not written for efficiency

char *human_chr[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
        "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
uint32_t n_human_chr = 24;

// paf recs and telostat each assembly contig
typedef struct {

    //from paf
    paf_rec_t **paf_recs;
    size_t n_paf_recs;
    size_t c_paf_recs;

    uint32_t len; //length of contig

    uint32_t ntelo; //updated using telo.bed
    char *mapped_chr; //best matched (updated using report.tsv)

} as_ctg_t;

// paf recs per each reference chromosome
typedef struct {
    //paf_rec_t **paf_recs; //not yet used
    //size_t n_paf_recs;  //not yet used
    uint32_t len; //length of chr
} as_chr_t;

// assembly contig name to pafrec hash
KHASH_MAP_INIT_STR(as_map_ctgs, as_ctg_t*)
// reference chromosome name to pafrec hash
KHASH_MAP_INIT_STR(as_map_chr, as_chr_t*)

KSORT_INIT_GENERIC(uint32_t)

static as_ctg_t *init_as_ctg(){
    as_ctg_t *ctg = (as_ctg_t *)malloc(sizeof(as_ctg_t));
    MALLOC_CHK(ctg);
    ctg->c_paf_recs = 5;
    ctg->paf_recs = (paf_rec_t **)malloc(sizeof(paf_rec_t*) * ctg->c_paf_recs);
    MALLOC_CHK(ctg->paf_recs);
    ctg->n_paf_recs = 0;
    ctg->len = 0; //length of the contig, will be set later
    ctg->ntelo = 0;
    ctg->mapped_chr = NULL; // initially no mapped chromosome
    return ctg;
}

static void free_as_ctg(as_ctg_t *ctg) {
    if (ctg) {
        for(int i=0; i<ctg->n_paf_recs; i++){
            paf_rec_t *rec = ctg->paf_recs[i];

            free(rec->rid);
            free(rec->tid);
            free(rec);

        }
        free(ctg->paf_recs);
        free(ctg->mapped_chr);
        free(ctg);
    }
}

static as_chr_t *init_as_chr() {
    as_chr_t *chr = (as_chr_t *)malloc(sizeof(as_chr_t));
    MALLOC_CHK(chr);
    //chr->paf_recs = NULL;
    //chr->n_paf_recs = 0;
    chr->len = 0;
    return chr;
}

static void free_as_chr(as_chr_t *chr) {
    if (chr) {
        //free(chr->paf_recs);
        free(chr);
    }
}

static void destroy_hashmap_ctgs(khash_t(as_map_ctgs) *h){
    for(khiter_t k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free_as_ctg(kh_value(h, k));
            free((char*)kh_key(h, k)); // free the key (contig name)
        }
    }
    kh_destroy(as_map_ctgs, h);
}

static void destroy_hashmap_chr(khash_t(as_map_chr) *h_chr){
    for (khint_t k = 0; k < kh_end(h_chr); ++k){
        if (kh_exist(h_chr, k)){
            free_as_chr(kh_value(h_chr, k));
            free((char*)kh_key(h_chr, k));
        }
    }
    kh_destroy(as_map_chr, h_chr);
}

void trim_mat_pat(char *chr){

    char *pch=strstr(chr, "_PATERNAL");
    if(pch != NULL) *pch='\0';

    pch=strstr(chr, "_MATERNAL");
    if(pch !=NULL) *pch='\0';

    NULL_CHK(chr);
    assert(strlen(chr) > 0);

}

static void load_paf(const char *paffile, khash_t(as_map_ctgs) *h_ctg, khash_t(as_map_chr) *h_chr, uint8_t trim) {

    // Initialize buffers
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char) * bufferSize);
    MALLOC_CHK(buffer);

    // Read PAF file
    FILE *fp = fopen(paffile, "r");
    F_CHK(fp, paffile);

    while (getline(&buffer, &bufferSize, fp) != -1) {

        paf_rec_t *rec = parse_paf_rec(buffer);
        if(trim){
            trim_mat_pat(rec->tid);
        }
        //get relevant contig
        khiter_t k = kh_get(as_map_ctgs, h_ctg, rec->rid);
        if (k != kh_end(h_ctg)) {
            as_ctg_t *ctg = kh_value(h_ctg, k);

            //update the contig length
            if (ctg->len == 0) {
                ctg->len = rec->qlen;
            } else if (ctg->len != rec->qlen) {
                ERROR("Contig '%s' has inconsistent lengths in PAF file", rec->rid);
                exit(EXIT_FAILURE);
            }

            //store the paf record
            if (ctg->n_paf_recs == ctg->c_paf_recs) {
                ctg->c_paf_recs *= 2;
                ctg->paf_recs = (paf_rec_t **)realloc(ctg->paf_recs, sizeof(paf_rec_t*) * ctg->c_paf_recs);
                MALLOC_CHK(ctg->paf_recs);
            }
            ctg->paf_recs[ctg->n_paf_recs++] = rec;

            //update chr len
            khiter_t k_chr = kh_get(as_map_chr, h_chr, rec->tid);
            if (k_chr != kh_end(h_chr)) {
                as_chr_t *chr = kh_value(h_chr, k_chr);
                if (chr->len == 0) {
                    chr->len = rec->tlen;
                } else if (chr->len != rec->tlen) {
                    ERROR("Chromosome '%s' has inconsistent lengths in PAF file", rec->tid);
                    exit(EXIT_FAILURE);
                }
            } else {
                WARNING("Chromosome '%s' in PAF file was not there in the tsv report or the telomere bed", rec->tid);
            }


        } else {
            WARNING("Contig '%s' in PAF file was not there in the tsv report or the telomere bed", rec->rid);
            free(rec->rid);
            free(rec->tid);
            free(rec);
        }




    }
    fclose(fp);

    free(buffer);
}


static void load_telobed(khash_t(as_map_ctgs) *h, const char *bedfile){

    FILE *bedfp = fopen(bedfile,"r");
    F_CHK(bedfp,bedfile);

    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);

    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;
    int64_t line_no = 0;

    while ((readlinebytes = getline(&buffer, &bufferSize, bedfp)) != -1) {

        char *ref = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(ref);
        int64_t beg=-1;
        int64_t end=-1;

        int ret=sscanf(buffer,"%s\t%ld\t%ld",ref,&beg, &end);
        if(ret!=3 || end<beg){
            ERROR("Malformed bed entry at line %ld",line_no);
            exit(EXIT_FAILURE);
        }

        if(beg<0 || end<0){
            ERROR("Malformed bed entry at %s:%ld. Coordinates cannot be negative",bedfile,line_no);
            exit(EXIT_FAILURE);
        }
        if(beg>=end){
            ERROR("Malformed bed entry at %s:%ld. start must be smaller than end coordinate",bedfile,line_no);
            exit(EXIT_FAILURE);
        }

        khiter_t k = kh_get(as_map_ctgs, h, ref);
        if (k == kh_end(h)){ //add
            int absent;
            khint_t k = kh_put(as_map_ctgs, h, ref, &absent);
            if(absent == 1){
                kh_key(h, k) = strdup(ref);
                as_ctg_t *ctg = init_as_ctg();
                ctg->ntelo++;
                kh_value(h, k) = ctg;

            }
            else{
                ERROR("Contig '%s' insertion failed", ref);
                exit(EXIT_FAILURE);
            }
        } else { //update
            as_ctg_t *ctg = kh_value(h, k);
            ctg->ntelo++;
        }

        free(ref);
        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    //*count = reg_i;
    return;
}

static void load_fixasm_report(khash_t(as_map_ctgs) *h_ctg, khash_t(as_map_chr) *h_chr, const char *reportfile) {

    FILE *reportfp = fopen(reportfile, "r");
    F_CHK(reportfp, reportfile);

    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);

    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;
    int64_t line_no = 0;

    while ((readlinebytes = getline(&buffer, &bufferSize, reportfp)) != -1) {

        char *chr = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(chr);
        char *ctg = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(ctg);

        int ret=sscanf(buffer,"%s\t%s",ctg,chr);
        if(ret!=2){
            ERROR("Malformed report entry at line %ld. Expected format: <ctg>\t<chr>",line_no);
            exit(EXIT_FAILURE);
        }

        // contig hash table
        khiter_t k = kh_get(as_map_ctgs, h_ctg, ctg);
        if (k == kh_end(h_ctg)){ //add
            int absent;
            khint_t k = kh_put(as_map_ctgs, h_ctg, ctg, &absent);
            if(absent == 1){
                kh_key(h_ctg, k) = strdup(ctg);
                as_ctg_t *ctg = init_as_ctg();
                ctg->mapped_chr = strdup(chr);
                kh_value(h_ctg, k) = ctg;

            }
            else{
                ERROR("Contig '%s' insertion failed", ctg);
                exit(EXIT_FAILURE);
            }
        } else { //update
            as_ctg_t *ctg = kh_value(h_ctg, k);
            ctg->mapped_chr = strdup(chr);
        }


        //chr hash table
        khiter_t k_chr = kh_get(as_map_chr, h_chr, chr);
        if (k_chr == kh_end(h_chr)){ //add
            int absent;
            khint_t k = kh_put(as_map_chr, h_chr, chr, &absent);
            if(absent == 1){
                kh_key(h_chr, k) = strdup(chr);
                as_chr_t *chr = init_as_chr();

                kh_value(h_chr, k) = chr;

            }
            else{
                ERROR("Chromosome '%s' insertion failed", chr);
                exit(EXIT_FAILURE);
            }
        } else { //update
            //as_chr_t *chr = kh_value(h_chr, k_chr);
            //chr->len = 0;
        }

        free(chr);
        free(ctg);
        line_no++;
    }

    fclose(reportfp);
    free(buffer);
}

static struct option long_options[] = {
    {"report", required_argument, 0, 'r'},          //0
    {"human-chr", no_argument, 0, 0},             //1
    {"trim-pat-mat", no_argument, 0, 0},                  //2
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {0, 0, 0, 0}
};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto asmstats <asm2ref.paf> <telomere.bed> -r <fixasm.report.tsv>\n");
    fprintf(fp_help,"   -r FILE                    report from fixasm\n");
    fprintf(fp_help,"   --human-chr                use inbuilt human chromosome names and order when printing report\n");
    //fprintf(fp_help,"   --trim-pat-mat                     trim chr name suffixes _MATERNAL and _PATERNAL from PAF\n");
    fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -h                         help\n");
}

static void telo_table(khash_t(as_map_chr) *h_chr, khash_t(as_map_ctgs) *h_ctg, char const **chr_list, size_t chr_list_size){

    fprintf(stdout,"chr\tT2T?\tNTelo\tTelocontiglen\n");

    //go through each chromosome
    for (int i = 0; i < chr_list_size; i++){

        int32_t total_telo = 0; // Total number of telomeres for this chr

        size_t n_contigs = 0;
        size_t capacity = 10;
        char *t2t = malloc(sizeof(char) * capacity);
        uint32_t *len = malloc(sizeof(uint32_t) * capacity);
        MALLOC_CHK(t2t);
        MALLOC_CHK(len);

        //iterate through the contig hash table as performance is not a key goal
        for (khiter_t k_ctg = kh_begin(h_ctg); k_ctg != kh_end(h_ctg); ++k_ctg) {
            if (kh_exist(h_ctg, k_ctg)) {
                as_ctg_t *ctg = kh_value(h_ctg, k_ctg);
                NULL_CHK(ctg);

                if (ctg->mapped_chr && strcmp(ctg->mapped_chr, chr_list[i]) == 0) {

                    if(ctg->ntelo>0){
                        char this_t2t = (ctg->ntelo == 2)? 'y' : 'n';
                        uint32_t this_len = ctg->len;

                        if (n_contigs >= capacity) {
                            capacity *= 2;
                            t2t = realloc(t2t, sizeof(char) * capacity);
                            len = realloc(len, sizeof(uint32_t) * capacity);
                            MALLOC_CHK(t2t);
                            MALLOC_CHK(len);
                        }
                        t2t[n_contigs] = this_t2t;
                        len[n_contigs] = this_len;
                        n_contigs++;

                        total_telo += ctg->ntelo;
                    }

                }
            }
        }

        fprintf(stdout,"%s\t", chr_list[i]);

        if(n_contigs > 0){

            assert(total_telo>0 && n_contigs<=total_telo);
            for(int j=0;j<n_contigs;j++){
                fprintf(stdout,"%c,", t2t[j]);
            }
            fprintf(stdout,"\t%d\t", total_telo);
            for(int j=0;j<n_contigs;j++){
                fprintf(stdout,"%d,", len[j]);
            }
        } else {
            fprintf(stdout,"\t\t");
        }

        fprintf(stdout,"\n");

        free(t2t);
        free(len);

    }

    return;

}

static void process_chr_and_print(khash_t(as_map_ctgs) *h_ctg, char const *chr, uint32_t len, uint8_t invert){

    uint32_t c0=0,c01=0,c1=0,c5=0,c10=0; //counts
    uint64_t s0=0,s01=0,s1=0,s5=0,s10=0; //sums

    //iterate through the contig hash table as performance is not a key goal
    for (khiter_t k_ctg = kh_begin(h_ctg); k_ctg != kh_end(h_ctg); ++k_ctg) {
        if (kh_exist(h_ctg, k_ctg)) {
            as_ctg_t *ctg = kh_value(h_ctg, k_ctg);
            NULL_CHK(ctg);

            // if contig is mapped to any chromome [ctg->mapped_chr exists] and
            // if invert = 0 then we should process contigs whose majority aligned to the currently considered chromosome
            // if invert = 1 then we should process contigs whose majority not aligned to the currently considered chromosome
            if (ctg->mapped_chr && (invert == 0 ? (strcmp(ctg->mapped_chr, chr) == 0) : (strcmp(ctg->mapped_chr, chr) != 0))) {

                uint64_t ta = 0; //total aligned length

                //get the paf records for this contig
                paf_rec_t **paf_recs = ctg->paf_recs;
                size_t n_paf_rec = ctg->n_paf_recs;
                if(n_paf_rec > 0){
                    for(size_t j = 0; j < n_paf_rec; j++){
                        paf_rec_t *rec = paf_recs[j];
                        NULL_CHK(rec);
                        if(strcmp(rec->tid,chr)==0){
                            //contig is aligned to the chromosome
                            uint32_t aligned_len = rec->target_end - rec->target_start; // length of the alignment
                            ta += aligned_len;
                        }
                    }
                    if(ta>0)        { c0++; s0+=ta;   }
                    if(ta>=100000)  { c01++; s01+=ta; }
                    if(ta>=1000000) { c1++; s1+=ta;   }
                    if(ta>=5000000) { c5++; s5+=ta;   }
                    if(ta>=10000000){ c10++; s10+=ta; }
                }
            }
        }
    }

    fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t", chr, c0,c01,c1,c5,c10);
    fprintf(stdout, "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", (double)s0/len*100,(double)s01/len*100,(double)s1/len*100,(double)s5/len*100,(double)s10/len*100);

    return;
}


static void process_lx_chr_and_print(khash_t(as_map_ctgs) *h_ctg, char const *chr, uint32_t len){

    size_t cap_aln_lens = 100;
    uint32_t *aln_lens = (uint32_t *)malloc(cap_aln_lens * sizeof(uint32_t));
    MALLOC_CHK(aln_lens);
    size_t n_aln_lens  = 0;

    //iterate through the contig hash table as performance is not a key goal
    for (khiter_t k_ctg = kh_begin(h_ctg); k_ctg != kh_end(h_ctg); ++k_ctg) {
        if (kh_exist(h_ctg, k_ctg)) {
            as_ctg_t *ctg = kh_value(h_ctg, k_ctg);
            NULL_CHK(ctg);

            //contig is majority aligned to the chromosome
            if (ctg->mapped_chr && strcmp(ctg->mapped_chr, chr) == 0) {
                //get the paf records for this contig
                paf_rec_t **paf_recs = ctg->paf_recs;
                size_t n_paf_rec = ctg->n_paf_recs;
                uint64_t ta = 0; //total aligned length

                if(n_paf_rec > 0){
                    for(size_t j = 0; j < n_paf_rec; j++){
                        paf_rec_t *rec = paf_recs[j];
                        NULL_CHK(rec);
                        if(strcmp(rec->tid,chr)==0){
                            //contig is aligned to the chromosome
                            uint32_t aligned_len = rec->target_end - rec->target_start; // length of the alignment
                            ta += aligned_len;
                        }
                    }
                    if(n_aln_lens >= cap_aln_lens){
                        cap_aln_lens *= 2;
                        aln_lens = (uint32_t *)realloc(aln_lens, cap_aln_lens * sizeof(uint32_t));
                        MALLOC_CHK(aln_lens);
                    }
                    aln_lens[n_aln_lens++] = ta;

                }
            }
        }
    }

    ks_mergesort(uint32_t, n_aln_lens, aln_lens, 0);

    uint32_t l50=0,l90=0,l95=0,l99=0;
    uint64_t CumCov1=0,CumCov2=0,CumCov3=0,CumCov4=0,CumCov5=0;
    uint64_t sum = 0;
    for(int i=0; i<n_aln_lens; i++){
        int j=n_aln_lens-i-1;
        sum += aln_lens[j];
        if(sum>=len*0.50 && l50==0){l50=i+1;}
        if(sum>=len*0.90 && l90==0){l90=i+1;}
        if(sum>=len*0.95 && l95==0){l95=i+1;}
        if(sum>=len*0.99 && l99==0){l99=i+1;}
        if(i<1) {CumCov1+=aln_lens[j];}
        if(i<2) {CumCov2+=aln_lens[j];}
        if(i<3) {CumCov3+=aln_lens[j];}
        if(i<4) {CumCov4+=aln_lens[j];}
        if(i<5) {CumCov5+=aln_lens[j];}
    }

    fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t",chr,l50,l90,l95,l99);
    fprintf(stdout, "%.3f,%.3f,%.3f,%.3f,%.3f\n", (double)CumCov1/len*100,(double)CumCov2/len*100,(double)CumCov3/len*100,(double)CumCov4/len*100,(double)CumCov5/len*100);

    free(aln_lens);
    return;
}


void contig_majority_common(khash_t(as_map_chr) *h_chr, khash_t(as_map_ctgs) *h_ctg, char const **chr_list, size_t chr_list_size, uint8_t invert, uint8_t lx){

    if(lx==1 && invert==1){
        ERROR("%s", "Not supported lx==1 && invert==1");
        exit(EXIT_FAILURE);
    }

    //go through each chromosome
    for (int i = 0; i < chr_list_size; i++){

        khiter_t k = kh_get(as_map_chr, h_chr, chr_list[i]);

        if (k != kh_end(h_chr)) {
            as_chr_t *chr = kh_value(h_chr, k);
            uint32_t chrlen = chr->len;
            if(chrlen == 0){
                ERROR("Failed to get chromosome %s length from hash table. Check your input files.", chr_list[i]);
                exit(EXIT_FAILURE);
            }
            if(lx){
                process_lx_chr_and_print(h_ctg, chr_list[i], chrlen);
            } else {
                process_chr_and_print(h_ctg, chr_list[i], chrlen, invert);
            }
        } else {
            WARNING("Failed to get chromosome %s from hash table. Ignoring.", chr_list[i]);
            fprintf(stdout, "%s\n", chr_list[i]);
        }
    }
}


static void contig_majority_correct_table(khash_t(as_map_chr) *h_chr, khash_t(as_map_ctgs) *h_ctg, char const **chr_list, size_t chr_list_size){

    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Contigs whose majority is mapped to the corresponding chromosome\n");
    fprintf(stdout, "\tNcontigsofsize>=KMbasealignedtochr\t\t\t\t\t%%ofchrsequencecoveredbycontigsofsize>=KMbase\n");
    fprintf(stdout, "chr\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\n");

    contig_majority_common(h_chr, h_ctg, chr_list, chr_list_size, 0, 0);

}


static void lx_contig_majority_correct_table(khash_t(as_map_chr) *h_chr, khash_t(as_map_ctgs) *h_ctg, char const **chr_list, size_t chr_list_size){

    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "LX of Contigs whose majority is mapped to the corresponding chromosome\n");
    fprintf(stdout, "\tL50\tL90\tL95\tL99\tCumCovN5\n");

    contig_majority_common(h_chr, h_ctg, chr_list, chr_list_size, 0, 1);

}


void contig_majority_wrong_table(khash_t(as_map_chr) *h_chr, khash_t(as_map_ctgs) *h_ctg, char const **chr_list, size_t chr_list_size){

    fprintf(stdout, "\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Contigs whose majority is mapped to another chromosome\n");
    fprintf(stdout, "\tNcontigsofsize>=KMbasealignedtochr\t\t\t\t\t%%ofchrsequencecoveredbycontigsofsize>=KMbase\n");
    fprintf(stdout, "chr\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\n");

    contig_majority_common(h_chr, h_ctg, chr_list, chr_list_size, 1, 0);

}

char const **get_chr_list(khash_t(as_map_chr) *h_chr, size_t *chr_list_size) {
    size_t c = 100;
    char const **chr_list = malloc(sizeof(char *) * c);
    MALLOC_CHK(chr_list);
    size_t n = 0;

    for (khiter_t k = kh_begin(h_chr); k != kh_end(h_chr); ++k) {
        if (kh_exist(h_chr, k)) {
            if (n >= c) {
                c *= 2;
                chr_list = realloc(chr_list, sizeof(char *) * c);
                MALLOC_CHK(chr_list);
            }
            chr_list[n] = kh_key(h_chr, k);
            n++;
        }
    }

    *chr_list_size = n;
    return chr_list;
}

// static void print_hashmap_ctgs(khash_t(as_map_ctgs) *h_ctg){
//     for (khiter_t k = kh_begin(h_ctg); k != kh_end(h_ctg); ++k) {
//         if (kh_exist(h_ctg, k)) {
//             as_ctg_t *ctg = kh_value(h_ctg, k);
//             fprintf(stdout, "Contig: %s, Mapped Chr: %s; telos: %d\n", kh_key(h_ctg, k), ctg->mapped_chr, ctg->ntelo);
//         }
//     }
// }

int asmstats_main(int argc, char* argv[]) {

    const char* optstring = "r:h";

    int longindex = 0;
    int32_t c = -1;

    const char *paf = NULL;
    const char *bed = NULL;
    const char *report = NULL;

    uint8_t use_human_chr = 0;
    uint8_t trim = 0;

    FILE *fp_help = stderr;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        }
        else if (c == 'r') {
            report = optarg;
        } else if (c == 0 && longindex == 1){ // --human-chr{
            use_human_chr = 1;
        } else if (c == 0 && longindex == 2){ // --trim
            trim = 1;
        }


    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    paf = argv[optind];
    bed = argv[optind + 1];


    if (paf == NULL || bed == NULL || report == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    // Initialise hash table for assembly contigs
    khash_t(as_map_ctgs) *h_ctg = kh_init(as_map_ctgs);

    // Initialise hash table for reference chrs
    khash_t(as_map_chr) *h_chr = kh_init(as_map_chr);

    //load telobed
    load_telobed(h_ctg, bed);

    //load the contig to chr mapping (report from fixasm)
    load_fixasm_report(h_ctg, h_chr, report);

    //load the asm to ref paf
    load_paf(paf, h_ctg, h_chr, trim);

    size_t chr_list_size = 0;
    char const **chr_list = NULL;
    if(use_human_chr){
        chr_list = (char const **)human_chr;
        chr_list_size = n_human_chr;
    } else {
        chr_list = get_chr_list(h_chr, &chr_list_size);
    }

    fprintf(stdout, "%s\n\n", paf);

    telo_table(h_chr, h_ctg, chr_list, chr_list_size);

    contig_majority_correct_table(h_chr, h_ctg, chr_list, chr_list_size);

    lx_contig_majority_correct_table(h_chr, h_ctg, chr_list, chr_list_size);

    contig_majority_wrong_table(h_chr, h_ctg, chr_list, chr_list_size);

    if(!use_human_chr) free(chr_list);

    destroy_hashmap_ctgs(h_ctg);
    destroy_hashmap_chr(h_chr);

    return 0;

}