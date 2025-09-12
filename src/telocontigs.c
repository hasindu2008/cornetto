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
#include <zlib.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include "error.h"
#include "kseq.h"
#include "ksort.h"
#include "khash.h"
#include "misc.h"

KSEQ_INIT(gzFile, gzread)


static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto telocontigs <assembly.fasta> <telomere.bed>\n");
    //fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -h                         help\n");

}

typedef struct {
    char *name;
    uint64_t len; //length of contig
    uint32_t ntelo; //updated using telo.bed
} telo_ctg_t;

KHASH_MAP_INIT_STR(telo_map_ctgs, telo_ctg_t*)

static void load_telobed(telo_ctg_t *contigs, khash_t(telo_map_ctgs) *h, const char *bedfile){

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

        khiter_t k = kh_get(telo_map_ctgs, h, ref);
        if (k == kh_end(h)){ //add
            ERROR("Contig '%s' in bed file not found in fasta", ref);
            exit(EXIT_FAILURE);
        } else { //update
            telo_ctg_t *ctg = kh_value(h, k);
            ctg->ntelo++;
        }

        free(ref);
        line_no++;
    }

    VERBOSE("%ld bed entries loaded from %s", line_no, bedfile);

    fclose(bedfp);
    free(buffer);
    //*count = reg_i;
    return;
}

int compare_ctg_len_desc(const void *a, const void *b){
    const telo_ctg_t *ctg_a = (const telo_ctg_t *)a;
    const telo_ctg_t *ctg_b = (const telo_ctg_t *)b;
    if(ctg_a->len < ctg_b->len) return 1;
    else if(ctg_a->len > ctg_b->len) return -1;
    else return 0;
}


int telocontigs_main(int argc, char* argv[]) {

    const char* optstring = "h";

    int longindex = 0;
    int32_t c = -1;

    const char *fasta = NULL;
    const char *bed = NULL;

    FILE *fp_help = stderr;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
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
    fasta = argv[optind];

    if (fasta == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    bed = argv[optind+1];
    if (bed == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    khash_t(telo_map_ctgs) *h = kh_init(telo_map_ctgs);

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta, "r");
    F_CHK(fp,fasta);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    uint64_t n = 0;
    uint64_t cap = 100;
    telo_ctg_t *contigs = malloc(cap * sizeof(telo_ctg_t));
    MALLOC_CHK(contigs);

    //uint64_t sum = 0;

    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));
        if(n==cap){
            cap*=2;
            contigs = realloc(contigs, cap * sizeof(telo_ctg_t));
            MALLOC_CHK(contigs);
        }
        telo_ctg_t tmp;
        tmp.name = strdup(seq->name.s);
        tmp.len = l;
        tmp.ntelo = 0;
        contigs[n] = tmp;

        int absent;
        khint_t k = kh_put(telo_map_ctgs, h, tmp.name, &absent);
        if(absent == 1){
            kh_key(h, k) = tmp.name;
            kh_value(h, k) = &contigs[n];
        }
        else{
            ERROR("Duplicate contig '%s' found in fasta", tmp.name);
            exit(EXIT_FAILURE);
        }
        n++;

    }

    kseq_destroy(seq);
    gzclose(fp);

    //load the telobed
    load_telobed(contigs, h, bed);

    //sort by len
    qsort(contigs, n, sizeof(telo_ctg_t), compare_ctg_len_desc);

    //print the report
    fprintf(stdout, "Contig\tLength\tNTelomeres\n");
    for(uint64_t i=0; i<n; i++){
        fprintf(stdout, "%s\t%lu\t%u\n", contigs[i].name, contigs[i].len, contigs[i].ntelo);
    }

    for(int i=0; i<n; i++){
        free(contigs[i].name);
    }
    free(contigs);

    kh_destroy(telo_map_ctgs, h);

    return 0;

}