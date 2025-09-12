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
// KSORT_INIT_GENERIC(uint64_t)

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto report <assembly.fasta> ... \n");
    //fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -h                         help\n");

}

//in nx.c
void sort_contigs(uint64_t *length, uint64_t n);

int report_main(int argc, char* argv[]) {

    const char* optstring = "h";

    int longindex = 0;
    int32_t c = -1;

    const char *fasta = NULL;

    FILE *fp_help = stderr;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        }
    }

    // No arguments given
    if (argc - optind < 1 || fp_help == stdout) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    fprintf(stdout, "#asm\tNcontigs\tLargestcontig(Mbase)\tN50(Mbase)\tN90(Mbase)\n");


    while(optind < argc){
        fasta = argv[optind];
        optind++;
        if (fasta == NULL) {
            print_help_msg(fp_help);
            if(fp_help == stdout){
                exit(EXIT_SUCCESS);
            }
            exit(EXIT_FAILURE);
        }
        fprintf(stdout, "%s\t", fasta);

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

        uint64_t cumsum = 0;
        uint64_t n50 = 0;
        uint64_t n90 = 0;

        for (uint64_t i = 0; i < n; ++i) {

            uint64_t len = length[n-i-1];
            cumsum += len;
            if(cumsum >= sum * 0.5 && n50 == 0){
                n50 = len;
            }
            if(cumsum >= sum * 0.9 && n90 == 0){
                n90 = len;
            }

        }

        fprintf(stdout, "%ld\t%.3f\t%.3f\t%.3f\n", n, length[n-1]/1e6, n50/1e6, n90/1e6);

        free(length);

    }







    return 0;

}