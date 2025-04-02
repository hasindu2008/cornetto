/**
 * @file boringbits_main.c
 * @brief entry point to boringbits_main
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
#include "cornetto.h"
#include "error.h"
#include "misc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


static struct option long_options[] = {
    {"threads", required_argument, 0, 't'},        //0 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //1 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //2 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //3 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //4
    {"version", no_argument, 0, 'V'},              //5
    {"output",required_argument, 0, 'o'},          //6 output to a file [stdout]
    {"debug-break",required_argument, 0, 0},       //7 break after processing the first batch (used for debugging)
    {"profile-cpu",required_argument, 0, 0},       //8 perform section by section (used for profiling - for CPU only)
    {"accel",required_argument, 0, 0},             //9 accelerator
    {"qual",required_argument, 0, 'q'},             //10 cov-mq20.bg
    {"window-size",required_argument, 0, 'w'},     //11 window size
    {"window-inc",required_argument, 0, 'i'},       //12 window increment
    {"low-thresh",required_argument, 0, 'L'},       //13 lowt hreshold for depth
    {"high-thresh",required_argument, 0, 'H'},      //14 high threshold for depth
    {"low-mq-thresh",required_argument, 0, 'Q'},    //15 lowthreshold for mapq depth
    {"min-ctg-len",required_argument, 0, 'm'},      //16 min contig length
    {"edge-len",required_argument, 0, 'e'},         //17 edge length to ignore
    {0, 0, 0, 0}};


typedef struct {
    int window_size;
    int window_inc;

    float low_cov_thresh;
    float high_cov_thresh;
    float low_mq_cov_thresh;
    int min_ctg_len;
    int edge_len;

    uint64_t flag;              //flags
    int32_t batch_size;         //max reads loaded at once: K
    int64_t batch_size_bytes;   //max bytes loaded at once: B

    int32_t num_thread; //t
    int32_t debug_break;

    int8_t boring;

} optp_t;

static inline void print_help_msg(FILE *fp_help, optp_t opt){
    fprintf(fp_help,"Usage: cornetto boringbits cov-total.bg -q cov-mq20.bg\n");
    fprintf(fp_help,"\nbasic options:\n");
    //fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    //fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    //fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bytes loaded at once [%.1fM]\n",opt.batch_size_bytes/(float)(1000*1000));
    fprintf(fp_help,"   -q FILE                    depth file with high mapq read coverage\n");
    fprintf(fp_help,"   -w INT                     window size [%d]\n",opt.window_size);
    fprintf(fp_help,"   -i INT                     window increment [%d]\n",opt.window_inc);
    fprintf(fp_help,"   -L FLOAT                   low coverage threshold factor [%.1f]\n",opt.low_cov_thresh);
    fprintf(fp_help,"   -H FLOAT                   high coverage threshold factor [%.1f]\n",opt.high_cov_thresh);
    fprintf(fp_help,"   -Q FLOAT                   mapq low coverage threshold factor [%.1f]\n",opt.low_mq_cov_thresh);
    fprintf(fp_help,"   -m INT                     minimum contig length [%d]\n",opt.min_ctg_len);
    fprintf(fp_help,"   -e INT                     edge length to ignore [%d]\n",opt.edge_len);

    fprintf(fp_help,"   -h                         help\n");
    //fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    //fprintf(fp_help,"\nadvanced options:\n");
    //fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
    //fprintf(fp_help,"   --profile-cpu=yes|no       process section by section (used for profiling on CPU)\n");
#ifdef HAVE_ACC
    fprintf(fp_help,"   --accel=yes|no             Running on accelerator [%s]\n",(opt.flag&CORNETTO_ACC?"yes":"no"));
#endif

}

typedef struct {
    char *ctg_name;
    int ctg_length;
    int c_depth;
    uint16_t *depth;
    uint16_t *mq_depth;
} ctg_depth_t;

typedef struct {
    int num_ctg;
    int c_ctg;
    ctg_depth_t *ctg_depth;
    int mean_depth;
    int mean_mq_depth;
} asm_depth_t;


typedef struct {
    int st;
    int end;
    int depth;
    int mq_depth;
} reg_t;

typedef struct {
    char *ctg_name;
    int ctg_length;
    int n_reg;
    reg_t *reg;
} ctg_reg_t;

typedef struct {
    int num_ctg;
    ctg_reg_t *ctg_reg;
    int mean_depth;
    int mean_mq_depth;
} asm_reg_t;

asm_depth_t *init_asm_depth(){
    asm_depth_t *asm_depth = malloc(sizeof(asm_depth_t));
    MALLOC_CHK(asm_depth);

    asm_depth->num_ctg = 0;
    asm_depth->c_ctg = 1;
    asm_depth->ctg_depth = (ctg_depth_t*)malloc(asm_depth->c_ctg*sizeof(ctg_depth_t));
    MALLOC_CHK(asm_depth->ctg_depth);


    return asm_depth;
}


void free_asm_depth(asm_depth_t *asm_depth){
    for(int i=0;i<asm_depth->num_ctg;i++){
        free(asm_depth->ctg_depth[i].depth);
        free(asm_depth->ctg_depth[i].mq_depth);
        free(asm_depth->ctg_depth[i].ctg_name);
    }
    free(asm_depth->ctg_depth);
    free(asm_depth);
}


//lazily load all the shit memory
asm_depth_t *get_depths(const char* covtotalfile, const char *covmqfile){

    asm_depth_t  *asm_depth = init_asm_depth();

    FILE *fp_covtotal = fopen(covtotalfile, "r");
    F_CHK(fp_covtotal, covtotalfile);

    FILE *fp_covmq = fopen(covmqfile, "r");
    F_CHK(fp_covmq, covmqfile);

    char buff1[10000];
    char buff2[10000];
    int st1,st2;
    int end1,end2;
    int depth1,depth2;

    int ret = 0;
    char prev_ctg[10000] = "";
    int prev_pos = 0;

    double tot_ctg_len = 0;
    double tot_depth = 0;
    double tot_mq_depth = 0;

    while (1){
        if ((ret = fscanf(fp_covtotal, "%s\t%d\t%d\t%d\n", buff1, &st1, &end1, &depth1)) == EOF){
            break;
        }

        if (ret != 4){
            ERROR("The depth files should have 4 columns. Had %d.",ret);
            exit(EXIT_FAILURE);
        }

        if ((ret = fscanf(fp_covmq, "%s\t%d\t%d\t%d\n", buff2, &st2, &end2, &depth2)) == EOF){
            ERROR("%s","The two files are not in the same order");
            exit(EXIT_FAILURE);
        }

        if (ret != 4){
            ERROR("The depth files should have 4 columns. Had %d.",ret);
            exit(EXIT_FAILURE);
        }

        if (strcmp(buff1, buff2) != 0 || st1 != st2 || end1 != end2){
            ERROR("%s","The two files are not in the same order");
            exit(EXIT_FAILURE);
        }

        if(strcmp(buff1,prev_ctg)!=0){
            strcpy(prev_ctg,buff1);

            if(asm_depth->num_ctg==asm_depth->c_ctg){
                asm_depth->c_ctg = asm_depth->c_ctg*2;
                asm_depth->ctg_depth = (ctg_depth_t*)realloc(asm_depth->ctg_depth, asm_depth->c_ctg*sizeof(ctg_depth_t));
                MALLOC_CHK(asm_depth->ctg_depth);
            }

            asm_depth->ctg_depth[asm_depth->num_ctg].ctg_name = strdup(buff1);
            asm_depth->ctg_depth[asm_depth->num_ctg].ctg_length = 0;
            asm_depth->ctg_depth[asm_depth->num_ctg].c_depth = 100;
            asm_depth->ctg_depth[asm_depth->num_ctg].depth = (uint16_t*)calloc(asm_depth->ctg_depth[asm_depth->num_ctg].c_depth, sizeof(uint16_t));
            MALLOC_CHK(asm_depth->ctg_depth[asm_depth->num_ctg].depth);
            asm_depth->ctg_depth[asm_depth->num_ctg].mq_depth = (uint16_t*)calloc(asm_depth->ctg_depth[asm_depth->num_ctg].c_depth, sizeof(uint16_t));
            MALLOC_CHK(asm_depth->ctg_depth[asm_depth->num_ctg].mq_depth);
            asm_depth->num_ctg++;
            prev_pos = 0;

        } else {
            if(prev_pos+1 != st1){
                ERROR("The depth files should be incremantal at one base resolution. Found %d to %d",prev_pos,st1);
                exit(EXIT_FAILURE);
            }
            prev_pos++;
        }

        if(st1+1 != end1){
            ERROR("The depth files should have end=start+1. Found %d to %d",st1,end1);
            exit(EXIT_FAILURE);
        }

        if(depth1>65535){
            WARNING("The depth at %s:%d-%d was truncated to 65535. Found %d",buff1, st1, end1, depth1);
            depth1 = 65535;
        }
        if(depth2>65535){
            WARNING("The depth at %s:%d-%d was truncated to 65535. Found %d",buff2, st2, end2, depth2);
            depth2 = 65535;
        }

        ctg_depth_t *ctg_depth = &asm_depth->ctg_depth[asm_depth->num_ctg-1];
        if(ctg_depth->ctg_length==ctg_depth->c_depth){
            ctg_depth->c_depth = ctg_depth->c_depth*2;
            ctg_depth->depth = (uint16_t*)realloc(ctg_depth->depth, ctg_depth->c_depth*sizeof(uint16_t));
            MALLOC_CHK(ctg_depth->depth);
            ctg_depth->mq_depth = (uint16_t*)realloc(ctg_depth->mq_depth, ctg_depth->c_depth*sizeof(uint16_t));
            MALLOC_CHK(ctg_depth->mq_depth);
        }

        ctg_depth->depth[ctg_depth->ctg_length] = depth1;
        ctg_depth->mq_depth[ctg_depth->ctg_length] = depth2;
        ctg_depth->ctg_length++;

        tot_depth += depth1;
        tot_mq_depth += depth2;
        tot_ctg_len ++;

    }


    fclose(fp_covtotal);
    fclose(fp_covmq);

    asm_depth->mean_depth = round(tot_depth/tot_ctg_len);
    asm_depth->mean_mq_depth = round(tot_mq_depth/tot_ctg_len);

    //fprintf(stderr, "Mean depth: %f\n", tot_depth/tot_ctg_len);
    //fprintf(stderr, "Mean mq depth: %f\n", tot_mq_depth/tot_ctg_len);

    return asm_depth;

}

void print_depth(asm_depth_t *asm_depth){
    for(int i=0;i<asm_depth->num_ctg;i++){
        ctg_depth_t *ctg_depth = &asm_depth->ctg_depth[i];
        for(int j=0;j<ctg_depth->ctg_length;j++){
            printf("%s\t%d\t%d\t%d\n",ctg_depth->ctg_name, j, j+1, ctg_depth->depth[j]);
        }
    }
}

void print_depth_mq(asm_depth_t *asm_depth){
    for(int i=0;i<asm_depth->num_ctg;i++){
        ctg_depth_t *ctg_depth = &asm_depth->ctg_depth[i];
        for(int j=0;j<ctg_depth->ctg_length;j++){
            printf("%s\t%d\t%d\t%d\n",ctg_depth->ctg_name, j, j+1, ctg_depth->mq_depth[j]);
        }
    }
}


asm_reg_t *get_regs(asm_depth_t *asm_depth, int window_size, int window_inc){

    asm_reg_t *asm_reg = malloc(sizeof(asm_reg_t));
    MALLOC_CHK(asm_reg);

    asm_reg->num_ctg = asm_depth->num_ctg;
    asm_reg->ctg_reg = (ctg_reg_t*)malloc(asm_reg->num_ctg*sizeof(ctg_reg_t));
    MALLOC_CHK(asm_reg->ctg_reg);

    for(int i=0; i<asm_reg->num_ctg; i++){
        ctg_depth_t *ctg_depth = &asm_depth->ctg_depth[i];
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        ctg_reg->ctg_name = strdup(ctg_depth->ctg_name);
        int length = ctg_depth->ctg_length;
        ctg_reg->ctg_length = length;

        ctg_reg->n_reg = (length - window_size + window_inc -1) / window_inc + 1;
        ctg_reg->n_reg = (ctg_reg->n_reg < 1) ? 1 : ctg_reg->n_reg; //this could be dodgy, need to think about it
        VERBOSE("Number of windows for %s: %d, length %d, total size %ld",ctg_reg->ctg_name, ctg_reg->n_reg, length, ctg_reg->n_reg*sizeof(reg_t));
        ctg_reg->reg = (reg_t*)malloc(ctg_reg->n_reg*sizeof(reg_t));
        MALLOC_CHK(ctg_reg->reg);

        int st = 0;
        int end = 0;
        for(int j=0;j<ctg_reg->n_reg;j++){
            st = j*window_inc;
            end = st + window_size;
            if(end>length){
                end = length;
            }
            //fprintf(stderr,"%s %d %d\n",ctg_reg->ctg_name,st,end);
            assert(st<end);
            int depth = 0;
            int mq_depth = 0;
            for(int k=st;k<end;k++){
                depth += ctg_depth->depth[k];
                mq_depth += ctg_depth->mq_depth[k];
            }
            depth = depth/(end-st);
            mq_depth = mq_depth/(end-st);
            ctg_reg->reg[j].st = st;
            ctg_reg->reg[j].end = end;
            ctg_reg->reg[j].depth = depth;
            ctg_reg->reg[j].mq_depth = mq_depth;
        }

        assert(end == length);
        assert(st < end);

    }

    asm_reg->mean_depth = asm_depth->mean_depth;
    asm_reg->mean_mq_depth = asm_depth->mean_mq_depth;

    return asm_reg;

}

void free_asm_reg(asm_reg_t *asm_reg){
    for(int i=0;i<asm_reg->num_ctg;i++){
        free(asm_reg->ctg_reg[i].ctg_name);
        free(asm_reg->ctg_reg[i].reg);
    }
    free(asm_reg->ctg_reg);
    free(asm_reg);
}

void print_regs(asm_reg_t *asm_reg){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        for(int j=0;j<ctg_reg->n_reg;j++){
            printf("%s\t%d\t%d\t%d\t%d\n",ctg_reg->ctg_name, ctg_reg->reg[j].st, ctg_reg->reg[j].end, ctg_reg->reg[j].depth, ctg_reg->reg[j].mq_depth);
        }
    }
}

void print_reg_low_cov(asm_reg_t *asm_reg, int thresh){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        for(int j=0;j<ctg_reg->n_reg;j++){
            if(ctg_reg->reg[j].depth<thresh) printf("depth<%d\t%s\t%d\t%d\t%d\t%d\n",thresh,ctg_reg->ctg_name, ctg_reg->reg[j].st, ctg_reg->reg[j].end, ctg_reg->reg[j].depth, ctg_reg->reg[j].mq_depth);
        }
    }
}

void print_reg_high_cov(asm_reg_t *asm_reg, int thresh){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        for(int j=0;j<ctg_reg->n_reg;j++){
            if(ctg_reg->reg[j].depth>thresh) printf("depth>%d\t%s\t%d\t%d\t%d\t%d\n",thresh,ctg_reg->ctg_name, ctg_reg->reg[j].st, ctg_reg->reg[j].end, ctg_reg->reg[j].depth, ctg_reg->reg[j].mq_depth);
        }
    }
}

void print_reg_low_mq_cov(asm_reg_t *asm_reg, float thresh_fac){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        for(int j=0;j<ctg_reg->n_reg;j++){
            if(ctg_reg->reg[j].mq_depth / (double)ctg_reg->reg[j].depth < thresh_fac) printf("depth_mq/depth<%.1f\t%s\t%d\t%d\t%d\t%d\n",thresh_fac,ctg_reg->ctg_name, ctg_reg->reg[j].st, ctg_reg->reg[j].end, ctg_reg->reg[j].depth, ctg_reg->reg[j].mq_depth);
        }
    }
}

void print_fun_bits(asm_reg_t *asm_reg, int thresh_low_depth, int thresh_high_depth, float low_mq_cov_thresh, int edge_len, int min_ctg_len){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        int ctg_len = ctg_reg->ctg_length;
        if(ctg_len < min_ctg_len){ // small contigs are always fun
            printf("%s\t%d\t%d\t.\t.\n",ctg_reg->ctg_name, 0, min_ctg_len);
        } else {
            //print edges
            printf("%s\t%d\t%d\t.\t.\n",ctg_reg->ctg_name, 0, edge_len);
            printf("%s\t%d\t%d\t.\t.\n",ctg_reg->ctg_name, ctg_len-edge_len, ctg_len);

            for(int j=0;j<ctg_reg->n_reg;j++){
                int depth = ctg_reg->reg[j].depth;
                int mq_depth = ctg_reg->reg[j].mq_depth;
                if(depth<thresh_low_depth || depth>thresh_high_depth || (mq_depth /(double)depth) < low_mq_cov_thresh ){
                    printf("%s\t%d\t%d\t%d\t%d\n",ctg_reg->ctg_name, ctg_reg->reg[j].st, ctg_reg->reg[j].end, ctg_reg->reg[j].depth, ctg_reg->reg[j].mq_depth);
                }
            }
        }
    }
}


/*
boringbits subtool: print_boring_bits is deprecated - it was never used to begin with

It prints windows that meet 1, 2 and 2
1. contigs must be > [1MBase] in size
2. excluding [100kb] edge regions at each
3. Does not fall into any of the below category
   - windows with low coverage: [<0.6x] genome average
   - windows with high coverage: [>1.6x] genome average
   - windows with low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]

To merge the windows that are overlapping or adjacent (within 500bp), you may use bedtools as follows:
./cornetto boringbits test/cov-total.bg -q test/cov-mq20.bg  | cut -f 1,2,3 | bedtools sort | bedtools merge -d 500 > boringbits.bed
*/

void print_boring_bits(asm_reg_t *asm_reg, int thresh_low_depth, int thresh_high_depth, float low_mq_cov_thresh, int edge_len, int min_ctg_len){
    for(int i=0;i<asm_reg->num_ctg;i++){
        ctg_reg_t *ctg_reg = &asm_reg->ctg_reg[i];
        int ctg_len = ctg_reg->ctg_length;
        if(ctg_len > min_ctg_len){ //only large contigs can have boring bits
            for(int j=0;j<ctg_reg->n_reg;j++){
                int depth = ctg_reg->reg[j].depth;
                int mq_depth = ctg_reg->reg[j].mq_depth;
                int st = ctg_reg->reg[j].st;
                int end = ctg_reg->reg[j].end;
                if ( st > edge_len && end < ctg_len - edge_len) { // exclude edges
                    if(!(depth<thresh_low_depth || depth>thresh_high_depth || (mq_depth /(double)depth) < low_mq_cov_thresh)){
                        printf("%s\t%d\t%d\t%d\t%d\n",ctg_reg->ctg_name, st, end, depth, mq_depth);
                    }
                }
            }
        }
    }
}

void the_boring_bits(const char* covtotalfile, const char *covmqfile, optp_t *opt){

    double realtime0 = realtime();
    asm_depth_t *asm_depth = get_depths(covtotalfile, covmqfile);
    VERBOSE("Loaded depth files in %.2f seconds", realtime() - realtime0);

    int window_size = opt->window_size;
    int window_inc = opt->window_inc;
    float low_cov_thresh = opt->low_cov_thresh;
    float high_cov_thresh = opt->high_cov_thresh;
    float low_mq_cov_thresh = opt->low_mq_cov_thresh;
    int min_ctg_len = opt->min_ctg_len;
    int edge_len = opt->edge_len;

    fprintf(stderr,"Number of contigs: %d\n",asm_depth->num_ctg);
    fprintf(stderr,"Average depth: %d\n",asm_depth->mean_depth);
    fprintf(stderr,"Average mq depth: %d\n",asm_depth->mean_mq_depth);
    fprintf(stderr, "Window size: %d\n", window_size);
    fprintf(stderr, "Window increment: %d\n", window_inc);
    fprintf(stderr, "Low coverage threshold: %.1fx%d\n", low_cov_thresh, asm_depth->mean_depth);
    fprintf(stderr, "High coverage threshold: %.1fx%d\n", high_cov_thresh, asm_depth->mean_depth);
    fprintf(stderr, "Low mapq coverage threshold: %.1f\n", low_mq_cov_thresh);
    fprintf(stderr, "Min contig length: %d\n", min_ctg_len);
    fprintf(stderr, "Edge length: %d\n", edge_len);

    //print_depth_mq(asm_depth);

    realtime0 = realtime();
    asm_reg_t *asm_reg = get_regs(asm_depth, window_size, window_inc);
    VERBOSE("Found regions in %.2f seconds", realtime() - realtime0);

    free_asm_depth(asm_depth);


    //print_regs(asm_reg);
    int thresh_low_depth = round(low_cov_thresh * asm_reg->mean_depth);
    int thresh_high_depth = round(high_cov_thresh * asm_reg->mean_depth);

    //print_reg_low_cov(asm_reg, thresh_low_depth);

    //print_reg_high_cov(asm_reg, thresh_high_depth);

    //print_reg_low_mq_cov(asm_reg, low_mq_cov_thresh);

    realtime0 = realtime();
    if(opt->boring)
        print_boring_bits(asm_reg, thresh_low_depth, thresh_high_depth, low_mq_cov_thresh, edge_len, min_ctg_len);
    else
        print_fun_bits(asm_reg, thresh_low_depth, thresh_high_depth, low_mq_cov_thresh, edge_len, min_ctg_len);
    VERBOSE("Printed the bits in %.2f seconds", realtime() - realtime0);

    free_asm_reg(asm_reg);

}



static void init_optp(optp_t *opt){
    opt->window_size = 2500;
    opt->window_inc = 50;
    opt->low_cov_thresh = 0.4;
    opt->high_cov_thresh = 2.5;
    opt->low_mq_cov_thresh = 0.4;
    opt->min_ctg_len = 1000000;
    opt->edge_len = 100000;

    opt->batch_size = 512;
    opt->batch_size_bytes = 1000*1000*1000;
    opt->num_thread = 8;
    opt->debug_break = 0;
    opt->flag = 0;

    opt->boring = 1;
}

int boringbits_main(int argc, char* argv[], int8_t boring) {

    //double realtime0 = realtime();

    const char* optstring = "t:B:K:v:o:q:Q:H:L:w:i:e:m:hV";

    int longindex = 0;
    int32_t c = -1;

    const char *covtotalfile = NULL;
    const char *covmqfile = NULL;

    FILE *fp_help = stderr;

    optp_t opt;
    init_optp(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'B') {
            opt.batch_size_bytes = mm_parse_num(optarg);
            if(opt.batch_size_bytes<=0){
                ERROR("%s","Maximum number of bytes should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
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
        } else if (c=='q'){
            covmqfile = optarg;
        } else if (c=='w'){
            opt.window_size = atoi(optarg);
        } else if (c=='i'){
            opt.window_inc = atoi(optarg);
        } else if (c=='L'){
            opt.low_cov_thresh = atof(optarg);
        } else if (c=='H'){
            opt.high_cov_thresh = atof(optarg);
        } else if (c=='Q'){
            opt.low_mq_cov_thresh = atof(optarg);
        } else if (c=='m'){
            opt.min_ctg_len = atoi(optarg);
        } else if (c=='e'){
            opt.edge_len = atoi(optarg);
        } else if(c == 0 && longindex == 6){ //output
            //opt.output = optarg;
        } else if(c == 0 && longindex == 7){ //debug break
            opt.debug_break = atoi(optarg);
        } else if(c == 0 && longindex == 8){ //sectional benchmark todo : warning for gpu mode
            yes_or_no(&opt.flag, CORNETTO_PRF, long_options[longindex].name, optarg, 1);
        } else if(c == 0 && longindex == 9){ //accel
        #ifdef HAVE_ACC
            yes_or_no(&opt.flag, CORNETTO_ACC, long_options[longindex].name, optarg, 1);
        #else
            WARNING("%s", "--accel has no effect when compiled for the CPU");
        #endif
        }
    }
    opt.boring = boring;

    // No arguments given
    if (argc - optind != 1 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }
    covtotalfile = argv[optind];

    if (covtotalfile == NULL || covmqfile == NULL) {
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }


    the_boring_bits(covtotalfile, covmqfile, &opt);


    return 0;
}
