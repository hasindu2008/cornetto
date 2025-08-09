#ifndef PAF_REC_H
#define PAF_REC_H

#include <stdint.h>

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
    int32_t match_len;
    int32_t block_len;
    uint8_t mapq;
    char tp;
} paf_rec_t;

paf_rec_t *parse_paf_rec(char *buffer);

#endif // PAF_H