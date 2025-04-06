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

// Type definitions
typedef struct {
    char *rid;
    int32_t qlen;
    int32_t query_start;
    int32_t query_end;
    int8_t strand;
    char *tid;
    int32_t tlen;
    int32_t target_start;
    int32_t target_end;
    uint8_t mapq;
    char tp;
} paf_rec_t; //todo remove duplicate of paf_rec_t in src/fixdir_main.c


static paf_rec_t *parse_paf_rec(char *buffer) {
    char *pch = NULL;
    paf_rec_t *paf = (paf_rec_t *)malloc(sizeof(paf_rec_t));
    MALLOC_CHK(paf);

    // Read fields from buffer
    pch = strtok(buffer, "\t\r\n"); assert(pch != NULL);
    paf->rid = strdup(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->qlen = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->query_start = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->query_end = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->strand = (strcmp(pch, "+") == 0) ? 0 : 1;

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->tid = strdup(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->tlen = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->target_start = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->target_end = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);

    pch = strtok(NULL, "\t\r\n"); assert(pch != NULL);
    paf->mapq = atoi(pch);

    paf->tp = 'P';
    while ((pch = strtok(NULL, "\t\r\n"))) {
        if (strcmp("tp:A:P", pch) == 0) {
            paf->tp = 'P';
        } else if (strcmp("tp:A:S", pch) == 0) {
            paf->tp = 'S';
        }
    }

    return paf;
}


static void load_paf(const char *paffile) {

    // Initialize buffers
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char) * bufferSize);
    MALLOC_CHK(buffer);

    // Read PAF file
    FILE *fp = fopen(paffile, "r");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    while (getline(&buffer, &bufferSize, fp) != -1) {
        paf_rec_t *rec = parse_paf_rec(buffer);



        free(rec->rid);
        free(rec->tid);
        free(rec);
    }
    fclose(fp);

    free(buffer);
}

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto asmbed <assembly.fasta> \n");
    //fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -h                         help\n");
}

int asmstats_main(int argc, char* argv[]) {

    const char* optstring = "h";

    int longindex = 0;
    int32_t c = -1;

    const char *paf = NULL;
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
    paf = argv[optind];
    bed = argv[optind + 1];

    if (paf == NULL || bed == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    load_paf(paf);

    return 0;

}