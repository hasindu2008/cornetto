/**
 * @file fixasm_main.c
 * @brief entry point to fixasm
 * @author Kavindu Jayasooriya (k.jayasooriya@unsw.edu.au)

MIT License

Copyright (c) 2025 Kavindu Jayasooriya
Copyright (c) 2025 Hasindu Gamaarachchi

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
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"
#include "error.h"
#include "pafrec.h"

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"missing", required_argument, 0, 'm'},       //2 write missing contig names to file
    {"report", required_argument, 0, 'r'},        //3 write report to file
    {"trim-pat-mat", no_argument, 0, 'T'},       //4 trim paternal and maternal suffixes
    {0, 0, 0, 0}
};



static char* cleanup_str(char *str, uint8_t trim_suffixes){
    char *cleaned_str = (char *)malloc(strlen(str) + 1);
    strcpy(cleaned_str, str);
    MALLOC_CHK(cleaned_str);

    if(trim_suffixes) {
        char *pch=strstr(cleaned_str, "_PATERNAL");
        if(pch != NULL) *pch='\0';

        pch=strstr(cleaned_str, "_MATERNAL");
        if(pch !=NULL) *pch='\0';
    }

    return cleaned_str;
}

typedef struct {
    char* id;
    int64_t sump;
    int64_t sumn;

    // tally array: tally[i] gives how many times the chr pointed by chr_list_t.chr_list[i] and chr_list_t.chr_index[i] is in the paf file
    int tally_size; // size of the tally array
    int tally_capacity; // capacity of the tally array
    int *tally;

    //later used to save the deduced new contig name
    char *new_name; // this will be set to the chr name with the highest tally

} ctg_t;

typedef struct{
    int size;
    int capacity;
    char **names;
    int *counters;
} chr_list_t; //index to name and current printing counter mapping

KSEQ_INIT(gzFile, gzread)
// Hash table for assembly contigs
KHASH_MAP_INIT_STR(map_ctgs, ctg_t*)
// reference chr name to index mapping
KHASH_MAP_INIT_STR(map_chr, int)

// Function declarations
static inline void print_help_msg(FILE *fp_help);
void reverse_complement(kseq_t *seq);
int fixasm_main(int argc, char* argv[]);

// Function implementations
static inline void print_help_msg(FILE *fp_help) {
    fprintf(fp_help, "Usage: cornetto fixasm <assembly.fa> <asm_to_ref.paf>\n");
    fprintf(fp_help, "   -m FILE                    write missing contig names to FILE\n");
    fprintf(fp_help, "   -r FILE                    write report to FILE\n");
    fprintf(fp_help, "   -w FILE                    write fixed PAF to FILE\n");
    fprintf(fp_help, "   -T                         trim chr name suffixes _MATERNAL and _PATERNAL in outputs\n");
    fprintf(fp_help, "   -v INT                     verbosity level [%d]\n", (int)get_log_level());
    fprintf(fp_help, "   -h                         help\n");
}

static ctg_t *init_ctg(char *rid){
    ctg_t *new_ctg = (ctg_t *)malloc(sizeof(ctg_t));
    MALLOC_CHK(new_ctg);
    new_ctg->id = rid;
    new_ctg->sump = 0;
    new_ctg->sumn = 0;
    new_ctg->tally_size = 0;
    new_ctg->tally_capacity = 100;
    new_ctg->tally = (int *)malloc(sizeof(int) * new_ctg->tally_capacity);
    MALLOC_CHK(new_ctg->tally);
    for (int i = 0; i < new_ctg->tally_capacity; ++i) {
        new_ctg->tally[i] = 0;
    }
    new_ctg->new_name = NULL; // initially no new name is set
    return new_ctg;
}

static chr_list_t *init_chr_list(){
    chr_list_t *chr_list = (chr_list_t *)malloc(sizeof(chr_list_t));
    MALLOC_CHK(chr_list);
    chr_list->size = 0;
    chr_list->capacity = 100;

    chr_list->names = (char **)malloc(sizeof(char*) * chr_list->capacity);
    MALLOC_CHK(chr_list->names);
    chr_list->counters = (int *)malloc(sizeof(int) * chr_list->capacity);
    MALLOC_CHK(chr_list->counters);
    memset(chr_list->counters, 0, sizeof(int) * chr_list->capacity);

    return chr_list;
}

void insert_to_chr_list(chr_list_t *chr_list, char *tid){
    if (chr_list->size >= chr_list->capacity) {
        chr_list->capacity *= 2;
        chr_list->names = (char **)realloc(chr_list->names, sizeof(char*) * chr_list->capacity);
        MALLOC_CHK(chr_list);
        chr_list->counters = (int *)realloc(chr_list->counters, sizeof(int) * chr_list->capacity);
        MALLOC_CHK(chr_list->counters);
        memset(chr_list->counters + chr_list->size, 0, sizeof(int) * (chr_list->capacity - chr_list->size));
    }
    chr_list->names[chr_list->size] = tid;
    chr_list->size++;
}

void insert_to_ctg(chr_list_t *chr_list, ctg_t *cur_ctg, int this_chr_index){
    if (chr_list->size >= cur_ctg->tally_capacity) {
        cur_ctg->tally_capacity *= 2;
        cur_ctg->tally = (int *)realloc(cur_ctg->tally, sizeof(int) * cur_ctg->tally_capacity);
        MALLOC_CHK(cur_ctg->tally);
        for (int i = cur_ctg->tally_size; i < cur_ctg->tally_capacity; ++i) {
            cur_ctg->tally[i] = 0;
        }
    }
    assert(this_chr_index<chr_list->size);
    cur_ctg->tally_size = chr_list->size;
    cur_ctg->tally[this_chr_index]++;
}

void free_chr_list(chr_list_t *chr_list){
    free(chr_list->names);
    free(chr_list->counters);
    free(chr_list);
}

static void free_ctg(ctg_t *ctg) {
    free(ctg->id);
    free(ctg->tally);
    if (ctg->new_name) {
        free(ctg->new_name);
    }
    free(ctg);
}

static void destroy_hashmap_ctgs(khash_t(map_ctgs) *h){
    for(khiter_t k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free_ctg(kh_value(h, k));
        }
    }
    kh_destroy(map_ctgs, h);
}

static void destroy_hashmap_chr(khash_t(map_chr) *h_chr){
    for (khint_t k = 0; k < kh_end(h_chr); ++k){
        if (kh_exist(h_chr, k)){
            free((char*)kh_key(h_chr, k));
        }
    }
    kh_destroy(map_chr, h_chr);
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

void load_paf(const char *paffile, khash_t(map_ctgs) *h, khash_t(map_chr) *h_chr, chr_list_t *chr_list) {

    // Initialize buffers
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char) * bufferSize);
    MALLOC_CHK(buffer);


    // Read PAF file
    FILE *fp = fopen(paffile, "r");
    F_CHK(fp, paffile);

    int absent;
    int absent_chr;

    while (getline(&buffer, &bufferSize, fp) != -1) {
        paf_rec_t *rec = parse_paf_rec(buffer);

        // Check if the assembly contig is already in the hash table

        khiter_t k = kh_put(map_ctgs, h, rec->rid, &absent);
        if (absent == 1) {
            ctg_t *new_ctg = init_ctg(rec->rid);
            kh_value(h, k) = new_ctg;
        } else if (absent == -1) {
            fprintf(stderr, "Error: failed to insert key\n");
            exit(EXIT_FAILURE);
        } else {
            free(rec->rid);
        }

        // Check if the chr is already in the chrname->index hash table
        khiter_t k_chr = kh_put(map_chr, h_chr, rec->tid, &absent_chr);
        if (absent_chr == 1) {
            kh_value(h_chr, k_chr) = chr_list->size;
            insert_to_chr_list(chr_list, rec->tid);
        } else if (absent_chr == -1) {
            fprintf(stderr, "Error: failed to insert key\n");
            exit(EXIT_FAILURE);
        } else {
            free(rec->tid);
        }

        ctg_t *cur_ctg = kh_value(h, k);
        int32_t length = rec->target_end - rec->target_start;
        if (rec->strand == 0) {
            cur_ctg->sump += length;
        } else {
            cur_ctg->sumn += length;
        }
        int this_chr_index = kh_value(h_chr, k_chr);
        insert_to_ctg(chr_list, cur_ctg, this_chr_index);

        free(rec);
    }
    fclose(fp);

    free(buffer);
}


void write_corrected_paf(const char *out_paf, const char *paffile, khash_t(map_ctgs) *h) {

    // Initialize buffers
    size_t bufferSize = 4096;
    char *buffer = (char *)malloc(sizeof(char) * bufferSize);
    MALLOC_CHK(buffer);


    // Read PAF file
    FILE *fp = fopen(paffile, "r");
    F_CHK(fp, paffile);

    FILE *fw = fopen(out_paf, "w");
    F_CHK(fw, out_paf);

    while (getline(&buffer, &bufferSize, fp) != -1) {
        paf_rec_t *rec = parse_paf_rec(buffer);
        khiter_t k = kh_get(map_ctgs, h, rec->rid);
        if (k != kh_end(h)) {
            ctg_t *cur_ctg = kh_value(h, k);
            int8_t newdir = rec->strand;
            int32_t new_query_start = rec->query_start;
            int32_t new_query_end = rec->query_end;
            if (cur_ctg->sump < cur_ctg->sumn) { // relative strand is negative
                newdir = !newdir;
                new_query_start = rec->qlen - rec->query_end;
                new_query_end = rec->qlen - rec->query_start;
            }
            char *newname = cur_ctg->new_name;

            //print paf record
            fprintf(fw, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\ttp:A:%c\n",
                    newname, rec->qlen, new_query_start, new_query_end, newdir == 0 ? '+' : '-',
                    rec->tid, rec->tlen, rec->target_start, rec->target_end, rec->match_len, rec->block_len, rec->mapq, rec->tp);


        } else{
            fprintf(stderr, "Error: contig %s not found in hash table\n", rec->rid);
            exit(EXIT_FAILURE);
        }

        free(rec->rid);
        free(rec->tid);
        free(rec);
    }
    fclose(fp);
    fclose(fw);

    free(buffer);
}




void fix_the_assembly(const char *fastafile, khash_t(map_ctgs) *h, chr_list_t *chr_list, const char *missing_fn, const char *report_fn, uint8_t trim_suffixes) {
    // Read FASTA file and write corrected FASTA to stdout
    gzFile fp_fasta = gzopen(fastafile, "r");
    F_CHK(fp_fasta, fastafile);

    kseq_t *seq = kseq_init(fp_fasta);
    MALLOC_CHK(seq);

    FILE *fp_report = NULL;
    if(report_fn != NULL) {
        fp_report = fopen(report_fn, "w");
        F_CHK(fp_report, report_fn);
    }
    FILE *fp_missing = NULL;
    if(missing_fn != NULL) {
        fp_missing = fopen(missing_fn, "w");
        F_CHK(fp_missing, missing_fn);
    }

    int missing=0, total = 0, neg = 0;
    while (kseq_read(seq) >= 0) {
        khiter_t k = kh_get(map_ctgs, h, seq->name.s);
        if (k != kh_end(h)) {
            ctg_t *cur_ctg = kh_value(h, k);
            char dir = '+';
            if (cur_ctg->sump < cur_ctg->sumn) { // relative strand is negative
                reverse_complement(seq);
                dir = '-';
                neg++;
            }
            //get the chr tally array index max for the given contig
            int max_chr_index = -1;
            int max_chr_index_value = -1;
            for(int i = 0; i < cur_ctg->tally_size; ++i) {
                if (cur_ctg->tally[i] >= max_chr_index_value) {
                    max_chr_index_value = cur_ctg->tally[i];
                    max_chr_index = i;
                }
                //fprintf(stderr, "FFF:%s\t%d\t%d\n", seq->name.s, i, cur_ctg->tally[i]);
            }
            assert(max_chr_index >= 0);
            assert(max_chr_index < chr_list->size);

            char *cleaned_name = cleanup_str(chr_list->names[max_chr_index], trim_suffixes);

            //get the counter
            int cleaned_name_counter = chr_list->counters[max_chr_index];
            //save the new name
            cur_ctg->new_name = (char *)malloc(strlen(cleaned_name) + 100);
            MALLOC_CHK(cur_ctg->new_name);
            sprintf(cur_ctg->new_name, "%s_%d", cleaned_name, cleaned_name_counter);

            if (fp_report) fprintf(fp_report, "%s\t%s\t%c\t%s_%d\n", seq->name.s, cleaned_name, dir, cleaned_name, cleaned_name_counter);
            //change the name accordingly and write to the stat tsv as well
            fprintf(stdout, ">%s_%d\n%s\n", cleaned_name, cleaned_name_counter, seq->seq.s);
            free(cleaned_name);
            total++;
            chr_list->counters[max_chr_index]++; //increment the current index for the chr
        } else {
            if (fp_missing) fprintf(fp_missing, "%s\n", seq->name.s);
            missing++;
        }
    }
    fprintf(stderr, "total: %d\nnegative: %d\nmissing: %d\n", total, neg, missing);

    // Clean up
    kseq_destroy(seq);
    gzclose(fp_fasta);

    if(fp_report != NULL) {
        fclose(fp_report);
    }
    if(fp_missing != NULL) {
        fclose(fp_missing);
    }
}

int fixasm_main(int argc, char* argv[]) {

    const char* optstring = "v:r:m:w:Th";

    int longindex = 0;
    int32_t c = -1;

    const char *missing = NULL;
    const char *report = NULL;
    const char *out_paf = NULL;

    uint8_t trim_suffixes = 0;

    FILE *fp_help = stderr;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='m'){
            missing = optarg;
        } else if (c == 'r'){
            report = optarg;
        } else if (c == 'w'){
            out_paf = optarg;
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum log_level_opt)v);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        } else if (c == 'T') {
            trim_suffixes = 1;
        }

    }

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

    // Initialize hash table for assembly contigs
    khash_t(map_ctgs) *h = kh_init(map_ctgs);

    // Initialize hash table for chr name to index mapping
    khash_t(map_chr) *h_chr = kh_init(map_chr);

    // Initialize chr list and chr counter
    chr_list_t *chr_list = init_chr_list();

    // Load PAF file into hash tables
    load_paf(paffile, h, h_chr, chr_list);

    // Fix the assembly
    fix_the_assembly(fastafile, h, chr_list, missing, report, trim_suffixes);

    //todo write the corrected PAF file if requested
    if (out_paf) write_corrected_paf(out_paf, paffile, h);

    free_chr_list(chr_list);
    destroy_hashmap_ctgs(h);
    destroy_hashmap_chr(h_chr);

    return 0;
}
