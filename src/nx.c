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
#include "misc.h"

KSEQ_INIT(gzFile, gzread)
KSORT_INIT_GENERIC(uint64_t)

static struct option long_options[] = {
    {"genome-size", required_argument, 0, 'g'},     //0 genome size
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto nx <assembly.fasta> \n");
    //fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -g STR                     genome size (e.g. 3.1G). if unspecified, will use total contig length\n");
    fprintf(fp_help,"   -h                         help\n");

}

void sort_contigs(uint64_t *length, uint64_t n){
    ks_mergesort(uint64_t, n, length, 0);
}

int nx_main(int argc, char* argv[]) {

    const char* optstring = "g:h";

    int longindex = 0;
    int32_t c = -1;

    const char *fasta = NULL;

    FILE *fp_help = stderr;
    int64_t genome_size = -1;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if (c == 'g') {
            genome_size=mm_parse_num(optarg);
            if(genome_size<=0){
                ERROR("%s","Genome size should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        }

    }

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
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

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fasta, "r");
    F_CHK(fp,fasta);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    uint64_t n = 0;
    uint64_t cap = 100;
    uint64_t *length = malloc(cap * sizeof(uint64_t));
    MALLOC_CHK(length);

    uint64_t sum = 0;

    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));
        if(n==cap){
            cap*=2;
            length = realloc(length, cap * sizeof(uint64_t));
            MALLOC_CHK(length);
        }
        length[n++] = l;
        sum += l;
    }

    kseq_destroy(seq);
    gzclose(fp);

    sort_contigs(length, n);

    fprintf(stdout, "#x\tcontig_len\n");

    uint64_t cumsum = 0;
    double percent = 0;
    for (uint64_t i = 0; i < n; ++i) {

        fprintf(stdout, "%f\t%lu\n", percent, length[n-i-1]);
        cumsum += length[n-i-1];
        if(genome_size>0){
            percent = (double)cumsum/genome_size*100;
        } else {
            percent = (double)cumsum/sum*100;
        }
        fprintf(stdout, "%f\t%lu\n", percent, length[n-i-1]);

    }

    free(length);

    return 0;

}