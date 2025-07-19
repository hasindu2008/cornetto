/**
 * @file main.c
 * @brief entry point
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "misc.h"
#include "cornetto.h"

int depth_main(int argc, char* argv[]);
int fixasm_main(int argc, char* argv[]);
int minidot_main(int argc, char* argv[]);
int boringbits_main(int argc, char* argv[], int8_t boring);
int bigenough_main(int argc, char* argv[]);
int find_telomere_main(int argc, char* argv[]);
int telomere_windows_main(int argc, char* argv[]);
int telomere_breaks_main(int argc, char* argv[]);
int sdust_main(int argc, char *argv[]);
int assbed_main(int argc, char* argv[]);
int seq_main(int argc, char* argv[]);
int asmstats_main(int argc, char* argv[]);

int print_usage(FILE *fp_help){

    fprintf(fp_help,"Usage: cornetto <command> [options]\n\n");
    fprintf(fp_help,"commands:\n");
    fprintf(fp_help,"   create panel:\n");
    //fprintf(fp_help,"         boringbits      print boring bits in an assembly (deprecated)\n");
    fprintf(fp_help,"       noboringbits    print no boring bits in an assembly\n");
    fprintf(fp_help,"       bigenough       find contigs that have sufficient boring bits\n");
    //fprintf(fp_help,"         subtool2      do something\n");
    fprintf(fp_help,"   dotplot:\n");
    fprintf(fp_help,"       fixasm          fix the direction of contigs in an assembly\n");
    fprintf(fp_help,"       minidot         create dot plot (from https://github.com/lh3/miniasm)\n");
    fprintf(fp_help,"   telomere eval:\n");
    fprintf(fp_help,"       telowin         analyse telomere windows in a fasta file\n");
    fprintf(fp_help,"       telobreaks      find telomere breaks in a fasta file\n");
    fprintf(fp_help,"       telofind        find telomere sequences in a fasta file\n");
    fprintf(fp_help,"       sdust           symmetric DUST (https://github.com/lh3/sdust)\n");
    fprintf(fp_help,"       asmstats        calculate assembly statistics\n");
    fprintf(fp_help,"   misc:\n");
    fprintf(fp_help,"       fa2bed          create a bed file with assembly contig lengths\n");
    fprintf(fp_help,"       seq             extract reads equal to larger than a threshold from a fastq\n");
    fprintf(fp_help,"       --help, -h      print this help message\n");
    fprintf(fp_help,"       --version, -V   print version information\n");

    if(fp_help==stderr){
        return(EXIT_FAILURE);
    } else if(fp_help==stdout){
        return(EXIT_SUCCESS);
    } else {
        return(EXIT_FAILURE);
    }

}

int main(int argc, char* argv[]){

    double realtime0 = realtime();

    int ret=1;

    if(argc<2){
        return print_usage(stderr);
    } else if (strcmp(argv[1],"depth")==0){
        ret=depth_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"fixasm")==0){
        ret=fixasm_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"boringbits")==0){ //deprecated (this was never used)
        ret=boringbits_main(argc-1, argv+1, 1);
    } else if (strcmp(argv[1],"noboringbits")==0){
        ret=boringbits_main(argc-1, argv+1, 0);
    } else if (strcmp(argv[1],"telowin")==0) {
        ret=telomere_windows_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"telobreaks")==0) {
        ret=telomere_breaks_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"telofind")==0) {
        ret=find_telomere_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"minidot")==0){
        ret=minidot_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"bigenough")==0){
        ret=bigenough_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"sdust")==0){
        ret=sdust_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"fa2bed")==0){
        ret=assbed_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"seq")==0){
        ret=seq_main(argc-1, argv+1);
    } else if (strcmp(argv[1],"asmstats")==0){
        ret=asmstats_main(argc-1, argv+1);
    } else if(strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-V")==0){
        fprintf(stdout,"cornetto %s\n",CORNETTO_VERSION);
        exit(EXIT_SUCCESS);
    } else if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0){
        return print_usage(stdout);
    } else{
        fprintf(stderr,"[cornetto] Unrecognised command %s\n",argv[1]);
        return print_usage(stderr);
    }

    fprintf(stderr,"[%s] Version: %s\n", __func__, CORNETTO_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec; Peak RAM: %.3f GB\n\n",
            __func__, realtime() - realtime0, cputime(),peakrss() / 1024.0 / 1024.0 / 1024.0);

    return ret;
}
