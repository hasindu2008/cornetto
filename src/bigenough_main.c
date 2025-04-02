/**
 * @file bigenough_main.c
 * @brief entry point to bigenough_main
 * @author Hasindu Gamaarachchi (hasindu@unsw.edu.au)

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <stdint.h>
#include "cornetto.h"
#include "error.h"
#include "khash.h"

typedef struct {
    int start;
    int end;
    int covlen; //covered length by boring bits
} reg_t;

KHASH_MAP_INIT_STR(chr_map, reg_t)

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {"version", no_argument, 0, 'V'},              //5
    {"threshold", required_argument, 0, 'T'},      //6 threshold
    {0, 0, 0, 0}};

typedef struct {
    int threshold;
} optp_t;


static inline void print_help_msg(FILE *fp_help, optp_t opt){
    fprintf(fp_help,"Usage: cornetto bigenough [options] <assembly.bed> <boring.bed>\n");
    fprintf(fp_help,"   -T INT                     threshold for bigenough [0-100]\n");
    fprintf(fp_help,"   -v INT                     verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   --version                  print version\n");

}

static void init_optp(optp_t *opt){
    opt->threshold = 50;
}


void update_covlen(khash_t(chr_map) *h, const char *bedfile){

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

        //TODO can optimised though strtok etc later
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


        khint_t k = kh_get(chr_map, h, ref);
        if (k == kh_end(h)) {
            ERROR("Contig '%s' in %s is not found in assembly bed file", ref, bedfile);
            exit(EXIT_FAILURE);
        }

        reg_t *r = &kh_value(h, k);

        r->covlen += (end - beg);

        free(ref);
        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    //*count = reg_i;
    return;
}


void print_bigenough_bits(khash_t(chr_map) *h, const char *bedfile, int threshold){

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

        //TODO can optimised though strtok etc later
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


        khint_t k = kh_get(chr_map, h, ref);
        if (k == kh_end(h)) {
            ERROR("Contig '%s' in %s is not found in assembly bed file", ref, bedfile);
            exit(EXIT_FAILURE);
        }

        reg_t *r = &kh_value(h, k);

        if(r->covlen > (r->end - r->start) * threshold / 100){
            printf("%s\t%ld\t%ld\n", ref, beg, end);
        }

        free(ref);
        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    //*count = reg_i;
    return;
}

void read_bed_to_hashmap(khash_t(chr_map) *h, const char *bedfile){

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

        //TODO can optimised though strtok etc later
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
        if(beg != 0 ){
            ERROR("start coordinate should be 0 in the assembly chromosome bed. Not so at %s:%ld. ",bedfile,line_no);
            exit(EXIT_FAILURE);
        }

        int absent;

        khint_t k = kh_put(chr_map, h, ref, &absent);
        if (absent == -1 || absent == 0) {
            // error if read_id duplicated?
            ERROR("Contig '%s' is duplicated in %s", ref, bedfile);
            exit(EXIT_FAILURE);
        }

        kh_key(h, k) = strdup(ref);
        reg_t *r = &kh_value(h, k);

        r->start = beg;
        r->end = end;
        r->covlen = 0;

        free(ref);
        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    //*count = reg_i;
    return;
}

void destroy_hashmap(khash_t(chr_map) *h){
    for (khint_t k = 0; k < kh_end(h); ++k){
        if (kh_exist(h, k)){
            //reg_t *r = &kh_value(h, k);
            //fprintf(stderr,"%s\t%d\t%d\t%d\n", kh_key(h, k), r->start, r->end, r->covlen);
            free((char*)kh_key(h, k));
        }
    }
    kh_destroy(chr_map, h);
}

void bigenough_boringbits(const char *assbed, const char *boringbed, optp_t *opt){


    khash_t(chr_map) *h = kh_init(chr_map);

    read_bed_to_hashmap(h, assbed);

    update_covlen(h, boringbed);

    print_bigenough_bits(h, boringbed, opt->threshold);

    destroy_hashmap(h);


}

int bigenough_main(int argc, char* argv[]) {

    const char* optstring = "T:v:hV";

    int longindex = 0;
    int32_t c = -1;

    const char *assbed = NULL;
    const char *boringbed = NULL;

    FILE *fp_help = stderr;

    optp_t opt;
    init_optp(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c=='T'){
            int threshold = atoi(optarg);
            if(threshold<0 || threshold>100){
                ERROR("Threshold should be between 0 and 100. You entered %d",threshold);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum log_level_opt)v);
        } else if (c=='V'){
            fprintf(stdout,"cornetto %s\n",CORNETTO_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
            fp_help = stdout;
        }

    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    assbed = argv[optind];
    boringbed = argv[optind+1];

    if (assbed == NULL || boringbed == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    bigenough_boringbits(assbed, boringbed, &opt);



    return 0;

}