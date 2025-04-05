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


//#define _XOPEN_SOURCE 700
#include <zlib.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include "error.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"min-len", required_argument, 0, 'm'},       //1 min length
    {"help", no_argument, 0, 'h'},                 //2
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto seq <reads.fastq> \n");
    //fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -m INT                     min length [%d]\n",30000);
    fprintf(fp_help,"   -h                         help\n");
}

int seq_main(int argc, char* argv[]) {

    const char* optstring = "hm:";

    int longindex = 0;
    int32_t c = -1;

    const char *fastq = NULL;

    FILE *fp_help = stderr;
    int min_len = 30000;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if (c == 'm') {
            min_len = atoi(optarg);
            if (min_len < 0) {
                fprintf(stderr, "Error: min-len must be a positive integer\n");
                print_help_msg(fp_help);
                exit(EXIT_FAILURE);
            }
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[optind - 1]);
            print_help_msg(fp_help);
            exit(EXIT_FAILURE);
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
    fastq = argv[optind];

    if (fastq == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fastq, "r");
    F_CHK(fp,fastq);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    uint64_t before = 0;
    uint64_t after = 0;
    uint64_t before_n = 0;
    uint64_t after_n = 0;

    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));
        before += seq->seq.l;
        before_n ++;
        if(seq->seq.l >= min_len) {
            after += seq->seq.l;
            after_n ++;
            printf("@%s", seq->name.s);
            if (seq->comment.l){
                printf("\t%s", seq->comment.s);
            }
            printf("\n");
            printf("%s\n+\n%s\n", seq->seq.s, seq->qual.s);
        }
    }
    fprintf(stderr, "total reads: %lu\t%lu bases\t%.2f Gbases\n", before_n, before, before / 1e9);
    fprintf(stderr, "reads >= %d: %lu\t%lu bases\t%.2f Gbases\n", min_len, after_n, after, after / 1e9);

    kseq_destroy(seq);
    gzclose(fp);

    return 0;

}