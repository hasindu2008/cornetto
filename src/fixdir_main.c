/**
 * @file fixdir_main.c
 * @brief entry point to fixdir 
 * @author Kavindu Jayasooriya (k.jayasooriya@unsw.edu.au)

MIT License

Copyright (c) 2019 Kavindu Jayasooriya (k.jayasooriya@unsw.edu.au)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY kIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/
#define _XOPEN_SOURCE 700
#include "cornetto.h"
#include "khash.h"
#include "error.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include "kseq.h"

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
} paf_rec_t;

typedef struct {
    char* id;
    int32_t sump;
    int32_t sumn;
} ctg_t;

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(map_ctgs, ctg_t)

// Function declarations
static inline void print_help_msg(FILE *fp_help);
char *strdup(const char *src);
paf_rec_t *parse_paf_rec(char *buffer);
void reverse_complement(kseq_t *seq);
int fixdir_main(int argc, char* argv[]);

// Function implementations
static inline void print_help_msg(FILE *fp_help) {
    fprintf(fp_help, "Usage: cornetto fixdir <incorrect_assembly.fa> <a.paf>\n");
}

paf_rec_t *parse_paf_rec(char *buffer) {
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

void reverse_complement(kseq_t *seq) {
    int len = seq->seq.l;
    for (int i = 0; i < len / 2; ++i) {
        char tmp = seq->seq.s[i];
        seq->seq.s[i] = seq->seq.s[len - 1 - i];
        seq->seq.s[len - 1 - i] = tmp;
    }
    for (int i = 0; i < len; ++i) {
        switch (seq->seq.s[i]) {
            case 'A': seq->seq.s[i] = 'T'; break;
            case 'T': seq->seq.s[i] = 'A'; break;
            case 'G': seq->seq.s[i] = 'C'; break;
            case 'C': seq->seq.s[i] = 'G'; break;
            default: break;
        }
    }
}

int fixdir_main(int argc, char* argv[]) {
    FILE *fp_help = stderr;

    // Check arguments
    if (argc - optind != 2) {
        print_help_msg(fp_help);
        exit(EXIT_FAILURE);
    }

    const char *fastafile = argv[optind];
    const char *paffile = argv[optind + 1];

    if (paffile == NULL || fastafile == NULL) {
        print_help_msg(fp_help);
        exit(EXIT_FAILURE);
    }

    // Initialize buffers
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char) * bufferSize);
    MALLOC_CHK(buffer);

    // Initialize hash table
    khash_t(map_ctgs) *h = kh_init(map_ctgs);
    int absent;

    // Read PAF file
    FILE *fp = fopen(paffile, "r");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    while (getline(&buffer, &bufferSize, fp) != -1) {
        paf_rec_t *rec = parse_paf_rec(buffer);
        khiter_t k = kh_put(map_ctgs, h, rec->rid, &absent);
        if (absent) {
            ctg_t new_ctg = { .id = strdup(rec->rid), .sump = 0, .sumn = 0 };
            kh_value(h, k) = new_ctg;
        }

        ctg_t *cur_ctg = &kh_value(h, k);
        int32_t length = rec->target_end - rec->target_start;
        if (rec->strand == 0) {
            cur_ctg->sump += length;
        } else {
            cur_ctg->sumn += length;
        }

        free(rec->rid);
        free(rec->tid);
        free(rec);
    }
    fclose(fp);

    // Read FASTA file and create corrected FASTA
    gzFile fp_fasta = gzopen(fastafile, "r");
    F_CHK(fp_fasta, fastafile);

    kseq_t *seq = kseq_init(fp_fasta);
    MALLOC_CHK(seq);

    FILE *fp_corrected_fasta = fopen("corrected_contigs.fasta", "w");
    if (!fp_corrected_fasta) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    int total = 0, neg = 0;
    while (kseq_read(seq) >= 0) {
        khiter_t k = kh_get(map_ctgs, h, seq->name.s);
        if (k != kh_end(h)) {
            ctg_t *cur_ctg = &kh_value(h, k);
            if (cur_ctg->sump < cur_ctg->sumn) { // relative strand is negative
                reverse_complement(seq);
                neg++;
            }
            fprintf(fp_corrected_fasta, ">%s\n%s\n", seq->name.s, seq->seq.s);
            total++;
        }
    }
    printf("count: %d\nneg: %d\n", total, neg);

    // Clean up
    kseq_destroy(seq);
    gzclose(fp_fasta);
    fclose(fp_corrected_fasta);
    kh_destroy(map_ctgs, h);
    free(buffer);

    return 0;
}
