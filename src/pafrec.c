/*
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "error.h"
#include "pafrec.h"

char *strdup(const char *src);

void malformed_paf_rec(char *pch){
    if(pch == NULL){
        ERROR("%s", "Malformed PAF record. Exiting.");
        exit(EXIT_FAILURE);
    }
}

paf_rec_t *parse_paf_rec(char *buffer) {
    char *pch = NULL;
    paf_rec_t *paf = (paf_rec_t *)malloc(sizeof(paf_rec_t));
    MALLOC_CHK(paf);

    // Read fields from buffer
    pch = strtok(buffer, "\t\r\n"); malformed_paf_rec(pch);
    paf->rid = strdup(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->qlen = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->query_start = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->query_end = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->strand = (strcmp(pch, "+") == 0) ? 0 : 1;

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->tid = strdup(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->tlen = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->target_start = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->target_end = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->match_len = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
    paf->block_len = atoi(pch);

    pch = strtok(NULL, "\t\r\n"); malformed_paf_rec(pch);
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