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

static inline void print_help_msg(FILE *fp_help){
    fprintf(fp_help,"Usage: cornetto fixdir a.paf\n");
}

char *strDup(const char *src) {
    char *dst = malloc(strlen (src) + 1);  // Space for length plus nul
    if (dst == NULL) return NULL;          // No memory
    strcpy(dst, src);                      // Copy the characters
    return dst;                            // Return the new string
}

typedef struct{
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
}paf_rec_t;

typedef struct{
    char* id;
    int64_t sump;
    int64_t sumn;
}ctg_t;

KHASH_MAP_INIT_STR(map_ctgs, ctg_t)

paf_rec_t *parse_paf_rec(char *buffer){

    char *pch=NULL;

    paf_rec_t *paf = (paf_rec_t *)malloc(sizeof(paf_rec_t));
    MALLOC_CHK(paf);

    //read name
    pch = strtok (buffer,"\t\r\n"); assert(pch!=NULL);
    paf->rid = strDup(pch);

    //readlen
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->qlen = atoi(pch);

    //query start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_start = atoi(pch);

    //query end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->query_end= atoi(pch);

    //relative strand
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    if(strcmp(pch,"+")==0){
        paf->strand=0;
    }
    else if(strcmp(pch,"-")==0){
        paf->strand=1;
    }
    else{
        assert(0);
    }

    //targetname
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tid = strDup(pch);

    //target len
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->tlen = atoi(pch);

    //target start
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_start = atoi(pch);

    //target end
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->target_end= atoi(pch);

    //num residue
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //num block
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);

    //mapq
    pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
    paf->mapq = atoi(pch);

    paf->tp = 'P';
    while((pch = strtok(NULL,"\t\r\n"))){
        if(strcmp("tp:A:P",pch)==0){
            paf->tp = 'P';
        }
        else if (strcmp("tp:A:S",pch)==0){
            paf->tp = 'S';
        }
    }

    return paf;
}

int fixdir_main(int argc, char* argv[]) {

    const char *paffile = NULL;

    FILE *fp_help = stderr;

    // more arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    paffile = argv[optind];

    if (paffile == NULL) {
        print_help_msg(fp_help);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //buffers for getline
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char)*(bufferSize));
    MALLOC_CHK(buffer);
    int readlinebytes=1;

    int absent, is_missing;

    khash_t(map_ctgs) *h = kh_init(map_ctgs);  // allocate a hash table

    FILE *fp = fopen(paffile, "r");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    while ((readlinebytes = getline(&buffer, &bufferSize, fp)) != -1) {
        paf_rec_t *rec = parse_paf_rec(buffer);

        khiter_t k = kh_get(map_ctgs, h, rec->tid);
        if (k == kh_end(h)) {
            ctg_t new_ctg;
            new_ctg.id = strDup(rec->rid);
            new_ctg.sump = 0;
            new_ctg.sumn = 0;
            k = kh_put(map_ctgs, h, new_ctg.id, &absent);
            kh_value(h, k) = new_ctg;
        }

        ctg_t *cur_ctg = &kh_value(h, k);
        int64_t length = rec->target_end - rec->target_start;
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
    free(buffer);

    FILE *fp_plus = fopen("ctg_plus.txt", "w");
    FILE *fp_minus = fopen("ctg_minus.txt", "w");
    if (!fp_plus || !fp_minus) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    const char *key;
    ctg_t val;
    kh_foreach(h, key, val, {
        if (val.sump >= val.sumn) {
            fprintf(fp_plus, "%s\n", val.id);
        } else {
            fprintf(fp_minus, "%s\n", val.id);
        }
    });

    fclose(fp_plus);
    fclose(fp_minus);
    kh_destroy(map_ctgs, h);

    return 0;
}
